from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function
from amuse.community import LiteratureReferencesMixIn

from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

class ArepoInterface(
    CodeInterface,
    GravitationalDynamicsInterface,
    LiteratureReferencesMixIn):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="arepo_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        
    
class Arepo(GravitationalDynamics):

    def __init__(self, **options):
        GravitationalDynamics.__init__(self,  ArepoInterface(**options), **options)
    
