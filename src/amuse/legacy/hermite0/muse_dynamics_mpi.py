import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system

class Hermite(object):
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    
    
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
    
        
    
    t = legacy_global(name='t',id=20,dtype='d')
    
    dt_param = legacy_global(name='dt_param',id=21,dtype='d')
    
    dt_dia = legacy_global(name='dt_dia',id=22,dtype='d')
    
    eps2 = legacy_global(name='eps2',id=23,dtype='d')
    
    flag_collision = legacy_global(name='flag_collision',id=24,dtype='i')
            
    def __init__(self, convert_nbody = None):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker')
        
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
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
    def add_particle():
        function = RemoteFunction()  
        function.id = 5
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @core.legacy_function    
    def _get_state():
        function = RemoteFunction()  
        function.id = 8
        function.can_handle_array = True
        function.name = "get_state"
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
        
     
    def get_state(self,id):
        name_to_value = self._get_state(id)
        result = self.dynamics_state(name_to_value['id_out'])
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            setattr(result, x, name_to_value[x])  
        return result
           
    def add_star(self, star):
        id = star.id
        mass = self.convert_nbody.to_nbody(star.mass.value())
        position = self.convert_nbody.to_nbody(star.position.value())
        velocity = self.convert_nbody.to_nbody(star.velocity.value())
        
        mass = mass.number
        x = position.number[0]
        y = position.number[1]
        z = position.number[2]
        vx = velocity.number[0]
        vy = velocity.number[1]
        vz = velocity.number[2]
        radius = self.convert_nbody.to_nbody(star.radius.value()).number
        self.add_particle(id, mass, radius, x, y, z, vx, vy, vz)
        
    def update_star(self, star):
        state = self.get_state(star.id)
        time = self.convert_nbody.to_si( self.t | nbody_system.time)
        #star.mass.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.mass(state.mass)))
        star.position.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.length(numpy.array((state.x, state.y, state.z)))))
        star.velocity.set_value_at_time(time, self.convert_nbody.to_si(nbody_system.speed(numpy.array((state.vx, state.vy, state.vz)))))
        return star
         
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).number, 0.0)
        return result
        
    def set_eps2(self, eps2):
        self.eps2 = self.convert_nbody.to_nbody(eps2).number
    
    
    def add_particles(self, particles):
        for x in particles:
            self.add_star(x)
            
    def update_particles(self, particles):
        for x in particles:
            self.update_star(x)
            
    def xadd_particles(self, particles):
        mass = []
        x_ = []
        y = []
        z = []
        vx = []
        vy = []
        vz = []
        radius = []
        id = []
        for x in particles:
            #self.update_star(x)
            id.append(x.id)
            mass_ = self.convert_nbody.to_nbody(x.mass.value())
            position = self.convert_nbody.to_nbody(x.position.value())
            velocity = self.convert_nbody.to_nbody(x.velocity.value())
            
            mass.append(mass_.number)
            x_.append(position.number[0])
            y.append(position.number[1])
            z.append(position.number[2])
            vx.append(velocity.number[0])
            vy.append(velocity.number[1])
            vz.append(velocity.number[2])
            radius.append(self.convert_nbody.to_nbody(x.radius.value()).number)
        self._add_particle(id, mass, radius, x_, y, z, vx, vy, vz)
            
    def update_attributes(self, attributes):
        for id, x in attributes:
            if x.name == 'mass':
                self.set_mass(id, self.convert_nbody.to_nbody(x.value()).number)
        
    
  
