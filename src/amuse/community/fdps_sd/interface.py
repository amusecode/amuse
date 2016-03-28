from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class FDPSInterface(
        CodeInterface,
        #LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        #StoppingConditionInterface,
        #SinglePointGravityFieldInterface,
        ):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(
                self,
                name_of_the_worker="fdpsgravity_worker", 
                **options)
    
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    
    
class FDPS(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  FDPSGravityInterface())
    
