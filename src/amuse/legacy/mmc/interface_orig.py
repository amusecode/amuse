from amuse.legacy import *

class mmc2Interface(LegacyInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, name_of_the_worker="mmc2_worker", **keyword_arguments)
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        
    
class mmc2(CodeInterface):

    def __init__(self):
        CodeInterface.__init__(self,  mmc2Interface())
    
