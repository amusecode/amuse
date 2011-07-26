from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

class HuaynoInterface(CodeInterface,GravitationalDynamicsInterface):
    include_headers = ['worker_code.h']
    
    MODE_OPENCL='opencl'
    MODE_OPENMP='openmp'
        
    def name_of_worker(self,mode):
        if mode==self.MODE_OPENCL:
            return 'huayno_worker_cl'
        if mode==self.MODE_OPENMP:
            return 'huayno_worker_mp'
        return 'huayno_worker'
      
    def __init__(self, mode=None, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_worker(mode), **options) 
        
    @legacy_function    
    def new_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_potential_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def evolve_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    def set_eps2(self, e):
        return self.set_eps2_parameter(e)

    def get_eps2(self):
        return self.get_eps2_parameter()
    
class Huayno(GravitationalDynamics):

    class inttypes:
        # http://stackoverflow.com/questions/36932/whats-the-best-way-to-implement-an-enum-in-python
        SHARED2=1
        EXTRAPOLATE=5
        PASS_KDK=2
        PASS_DKD=7
        HOLD_KDK=3
        HOLD_DKD=8
        PPASS_DKD=9
        BRIDGE_KDK=4
        CC=11
        CC_KEPLER=12
        OK=13
        KEPLER=14
        SHARED4=15

    def __init__(self, convert_nbody = None, **options):
        legacy_interface = HuaynoInterface(**options)
#        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter", 
            "timestep_parameter", 
            "timestep parameter for gravity calculations", 
            default_value = 0.03 | units.none
        )

        object.add_method_parameter(
            "get_inttype_parameter",
            "set_inttype_parameter", 
            "inttype_parameter", 
            "integrator method to use", 
            default_value = 8 | units.none
        )

        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "current simulation time", 
            default_value = 0.0 | nbody_system.time
        )


    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_timestep_parameter",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_timestep_parameter",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_inttype_parameter",
            (),
            (units.none, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_inttype_parameter",
            (units.none, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )
        
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('PARAMETER_CHANGE_B', 'set_eps2_parameter')
        object.add_method('RUN', 'get_kinetic_energy')
        object.add_method('RUN', 'get_potential_energy')
