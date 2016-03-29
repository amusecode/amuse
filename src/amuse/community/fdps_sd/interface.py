from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class FDPSInterface(
        CodeInterface,
        LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        #StoppingConditionInterface,
        #SinglePointGravityFieldInterface,
        ):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(
                self,
                name_of_the_worker="fdps_sd_worker", 
                **options)
    
#    @legacy_function
#    def get_number_of_particles():
#        function = LegacyFunctionSpecification()  
#        function.can_handle_array = True 
#        function.addParameter('value', dtype='int32', direction=function.OUT)
#        function.result_type = 'int32'
#        return function    
    
    
class FDPS(GravitationalDynamics, GravityFieldCode):

    __interface__ = FDPSInterface

    def __init__(self, convert_nbody = None, **options):

        legacy_interface = self.__interface__(**options)

        GravitationalDynamics.__init__(
                self,
                legacy_interface,
                convert_nbody,
                **options
                )

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)
        #self.stopping_conditions.define_state(object)

    def define_parameters(self, object):
        object.add_method_parameter(
                "get_begin_time",
                "set_begin_time",
                "begin_time",
                "model time to start the simulation at",
                default_value = 0.0 | nbody_system.time
                )

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

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
        
        #self.stopping_conditions.define_particle_set(object)
