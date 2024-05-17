from amuse.rfi.core import CodeInterface, legacy_function, LegacyFunctionSpecification
from amuse.support.interface import InCodeComponentImplementation


class ForTestingInterface(CodeInterface):
    include_headers = ['c_interface.h']

    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_long_long_int():
        function = LegacyFunctionSpecification()
        function.addParameter('in', dtype='int64', direction=function.IN)
        function.addParameter('out', dtype='int64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_double():
        function = LegacyFunctionSpecification()
        function.addParameter('double_in', dtype='float64', direction=function.IN)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_float():
        function = LegacyFunctionSpecification()
        function.addParameter('float_in', dtype='float32', direction=function.IN)
        function.addParameter('float_out', dtype='float32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_strings():
        function = LegacyFunctionSpecification()
        function.addParameter('string_inout1', dtype='string', direction=function.INOUT)
        function.addParameter('string_inout2', dtype='string', direction=function.INOUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_string_int():
        function = LegacyFunctionSpecification()
        function.addParameter('inint', dtype='int32', direction=function.IN)
        function.addParameter('in', dtype='string', direction=function.IN, default="echo")
        function.addParameter('out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_string_two():
        function = LegacyFunctionSpecification()
        function.addParameter('in1', dtype='string', direction=function.IN)
        function.addParameter('in2', dtype='string', direction=function.IN, default="echo")
        function.addParameter('out1', dtype='string', direction=function.OUT)
        function.addParameter('out2', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_array():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = None
        function.must_handle_array = True
        return function

    @legacy_function
    def echo_array_with_result():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    # @legacy_function
    def return_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'string'
        function.can_handle_array = True
        return function

    @legacy_function
    def echo_2_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN)
        function.addParameter('int_in2', dtype='int32', direction=function.IN, default=1)
        function.addParameter('int_out1', dtype='int32', direction=function.OUT)
        function.addParameter('int_out2', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def echo_3_int():
        function = LegacyFunctionSpecification()
        function.addParameter('i', dtype='int32', direction=function.IN)
        function.addParameter('j', dtype='int32', direction=function.IN)
        function.addParameter('k', dtype='int32', direction=function.IN)
        function.addParameter('l', dtype='int32', direction=function.IN, default=0)
        function.addParameter('m', dtype='int32', direction=function.IN, default=1)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def echo_inout_array_with_result():
        function = LegacyFunctionSpecification()
        function.addParameter('in_out', dtype='int32', direction=function.INOUT)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def echo_logical():
        function = LegacyFunctionSpecification()
        function.addParameter('input', dtype='bool', direction=function.IN)
        function.addParameter('output', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def print_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def dummy_3_int():
        function = LegacyFunctionSpecification()
        function.addParameter('i', dtype='int32', direction=function.IN)
        function.addParameter('j', dtype='int32', direction=function.IN)
        function.addParameter('k', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def print_error_string():
        function = LegacyFunctionSpecification()
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def sum_doubles():
        function = LegacyFunctionSpecification()
        function.addParameter('double_in1', dtype='float64', direction=function.IN)
        function.addParameter('double_in2', dtype='float64', direction=function.IN, default=1.0)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
