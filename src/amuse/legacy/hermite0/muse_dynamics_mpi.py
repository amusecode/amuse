import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction

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
    
    class double_property(object):
        def __init__(self, tag):
            self.tag = tag
        def __get__(self, instance, owner):
            if instance is None:
                return self
            id, int_result, double_result = instance.do_call(self.tag)
            return double_result[0]
        def __set__(self, instance, value):
            if instance is None:
                return self
            instance.do_call(self.tag, doubles_arg=[value])
        
    class int_property(object):
        def __init__(self, tag):
            self.tag = tag
        def __get__(self, instance, owner):
            if instance is None:
                return self
            id, int_result, double_result = instance.do_call(self.tag, int_arg2 = 0)
            return int_result
        def __set__(self, instance, value):
            if instance is None:
                return self
            instance.do_call(self.tag, int_arg1 = value, int_arg2 = 1)
        
    t = double_property(20)
    dt_param = double_property(21)
    dt_dia = double_property(22)
    eps2 = double_property(23)
    flag_collision = int_property(24)
            
    def __init__(self):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker')
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
        self.channel = core.MpiChannel(self.intercomm)
        
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
        function.id = 3
        function.addParameter('time-end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    def reinitialize_particles(self):
        id, int_result, double_result = self.do_call(4, doubles_arg=[time])
        return int_result
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
        
    
  
