from amuse.community import *

class supportInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="support_worker", **keyword_arguments)
    
    @legacy_function
    def add():
        function = LegacyFunctionSpecification()  
        function.addParameter('term1', dtype='float64', direction=function.IN)
        function.addParameter('term2', dtype='float64', direction=function.IN)
        function.addParameter('sum', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def many_points_on_sphere():
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='float64', direction=function.INOUT)
        function.addParameter('y', dtype='float64', direction=function.INOUT)
        function.addParameter('z', dtype='float64', direction=function.INOUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def rnd_points_on_sphere():
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='float64', direction=function.INOUT)
        function.addParameter('y', dtype='float64', direction=function.INOUT)
        function.addParameter('z', dtype='float64', direction=function.INOUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function


class support(InCodeComponentImplementation):

    def __init__(self):
        CodeInterface.__init__(self,  supportInterface())
    
