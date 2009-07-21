import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system

class Hermite(object):
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
    flag_collision = legacy_global(name='t',id=24,dtype='i')
            
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
        return function;
        
    def __del__(self):
        self.stop_worker()

         
        
    def do_call(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_arg=[], ints_arg=[]):
        #self.send_request(tag,id, int_arg1, int_arg2, doubles_arg, ints_arg)
        #return self.recieve_result()
        return (0,[],[])
   
    @core.legacy_function   
    def setup_module():
        function = RemoteFunction()  
        function.id = 1
        function.result_type = 'i'
        return function;
    
    
    @core.legacy_function      
    def cleanup_module():
        function = RemoteFunction()  
        function.id = 2
        function.result_type = 'i'
        return function;
    
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
        return function;
    @core.legacy_function    
    def _get_state():
        function = RemoteFunction()  
        function.id = 8
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('id_out', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;
        
    @core.legacy_function    
    def evolve():
        function = RemoteFunction()  
        function.id = 6
        function.addParameter('time-end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @core.legacy_function  
    def reinitialize_particles():
        function = RemoteFunction()  
        function.id = 4
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    def add_particle(self, state):
        return self._add_particle(state.id, **state.to_keyword_args())
    def get_state(self,id):
        (doubles,ints) = self._get_state(id)
        return self.dynamics_state(ints[0], doubles)
    @core.legacy_function   
    def get_number():
        function = RemoteFunction()  
        function.id = 7
        function.result_type = 'i'
        return function;
        
    def add_star(self, star):
        state = self.dynamics_state()
        state.id = star.id
        mass = self.convert_nbody.to_nbody(star.mass)
        position = self.convert_nbody.to_nbody(star.position)
        velocity = self.convert_nbody.to_nbody(star.velocity)
        
        state.mass = mass.number
        state.x = position.number[0]
        state.y = position.number[1]
        state.z = position.number[2]
        state.vx = velocity.number[0]
        state.vy = velocity.number[1]
        state.vz = velocity.number[2]
        state.radius = self.convert_nbody.to_nbody(star.radius).number
        self.add_particle(state)
        
    def update_star(self, star):
        state = self.get_state(star.id)
        
        star.mass = self.convert_nbody.to_si(nbody_system.mass(state.mass))
        star.position = self.convert_nbody.to_si(nbody_system.length(numpy.array((state.x, state.y, state.z))))
        star.velocity = self.convert_nbody.to_si(nbody_system.speed(numpy.array((state.vx, state.vy, state.vz))))
        return star
        
    
  
