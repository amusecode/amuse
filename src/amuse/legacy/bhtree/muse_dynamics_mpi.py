import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system

class BHTree(object):
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    extra_content = """
    int _add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
        dynamics_state state;
        state.id = id;
        state.mass = mass;
        state.radius = radius;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        return add_particle(state);
    }
    
    void _get_state(int id, int * id_out,  double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) {
        dynamics_state state = get_state(id);
        *id_out = state.id;
        *mass = state.mass;
        *radius = state.radius;
        *x = state.x;
        *y = state.y;
        *z = state.z;
        *vx = state.vx;
        *vy = state.vy;
        *vz = state.vz;
    }
    """
    
    class dynamics_state(object):
        _attributes = ['mass','radius','x','y','z','vx','vy','vz']
        def __init__(self, id = 0, doubles = [0.0 for x in range(8)]):
            self.id = id
            for i, name in enumerate(self._attributes):
                setattr(self, name, doubles[i])
                
        def to_doubles(self):
            result = [0.0 for x in range(8)]
            for i, name in enumerate(self._attributes):
                result[i] = getattr(self, name)
            return result
            
        def to_keyword_args(self):
            result = {}
            for i, name in enumerate(self._attributes):
                result[name] = getattr(self, name)
            return result
    
        
    
    timestep = legacy_global(name='timestep',id=21,dtype='d')
    eps2_for_gravity = legacy_global(name='eps2_for_gravity',id=22,dtype='d')
    theta_for_tree = legacy_global(name='theta_for_tree',id=23,dtype='d')
    
    use_self_gravity = legacy_global(name='use_self_gravity',id=24,dtype='i')
    ncrit_for_tree = legacy_global(name='ncrit_for_tree',id=25,dtype='i')
    
    dt_dia = legacy_global(name='dt_dia',id=246,dtype='d')
    
    
            
    def __init__(self, convert_nbody = None):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker')
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
        print "rank_parent",self.intercomm.rank
        self.channel = core.MpiChannel(self.intercomm)
        self.convert_nbody = convert_nbody
        
    def __del__(self):
        self.stop_worker()
        
    @core.legacy_function
    def stop_worker():
        function = RemoteFunction()  
        function.id = 0
        return function

    @core.legacy_function   
    def setup_module():
        function = RemoteFunction()  
        function.id = 1
        function.result_type = 'i'
        return function
    
    
    @core.legacy_function      
    def cleanup_module():
        function = RemoteFunction()  
        function.id = 2
        function.result_type = 'i'
        return function
    
    @core.legacy_function    
    def initialize_particles():
        function = RemoteFunction()  
        function.id = 3
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
        
    @core.legacy_function    
    def _add_particle():
        function = RemoteFunction()  
        function.id = 5
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @core.legacy_function    
    def _get_state():
        function = RemoteFunction()  
        function.id = 8
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('id_out', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = None
        return function
        
    @core.legacy_function    
    def evolve():
        function = RemoteFunction()  
        function.id = 6
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @core.legacy_function  
    def reinitialize_particles():
        function = RemoteFunction()  
        function.id = 4
        function.result_type = 'i'
        return function
        
    @core.legacy_function   
    def get_number():
        function = RemoteFunction()  
        function.id = 7
        function.result_type = 'i'
        return function;
     
    @core.legacy_function
    def set_mass():
        function = RemoteFunction()  
        function.id = 9
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;
    
    @core.legacy_function
    def get_time():
        function = RemoteFunction()  
        function.id = 10
        function.result_type = 'd'
        return function;
     
    def add_particle(self, state):
        return self._add_particle(state.id, **state.to_keyword_args())
        
    def get_state(self,id):
        name_to_value = self._get_state(id)
        result = self.dynamics_state(name_to_value['id_out'])
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            setattr(result, x, name_to_value[x])  
        return result
           
    def add_star(self, star):
        state = self.dynamics_state()
        state.id = star.id
        mass = self.convert_nbody.to_nbody(star.mass.value())
        position = self.convert_nbody.to_nbody(star.position.value())
        velocity = self.convert_nbody.to_nbody(star.velocity.value())
        
        state.mass = mass.number
        state.x = position.number[0]
        state.y = position.number[1]
        state.z = position.number[2]
        state.vx = velocity.number[0]
        state.vy = velocity.number[1]
        state.vz = velocity.number[2]
        state.radius = self.convert_nbody.to_nbody(star.radius.value()).number
        self.add_particle(state)
        
    def update_star(self, star):
        state = self.get_state(star.id)
        time = self.convert_nbody.to_si( self.get_time() | nbody_system.time)
        #star.mass.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.mass(state.mass)))
        star.position.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.length(numpy.array((state.x, state.y, state.z)))))
        star.velocity.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.speed(numpy.array((state.vx, state.vy, state.vz)))))
        return star
         
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).number, 1)
        return result
        
    def set_eps2(self, eps2):
        self.eps2_for_gravity = self.convert_nbody.to_nbody(eps2).number
    
    
    def add_particles(self, particles):
        for x in particles:
            self.add_star(x)
            
    def update_particles(self, particles):
        for x in particles:
            self.update_star(x)
            
    def update_attributes(self, attributes):
        for id, x in attributes:
            if x.name == 'mass':
                self.set_mass(id, self.convert_nbody.to_nbody(x.value()).number)
        
    
  
