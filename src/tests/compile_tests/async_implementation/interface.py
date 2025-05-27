from amuse.rfi.core import *
from compile_tests.c_implementation import interface as c_interface
from amuse.units import units


class ForTestingInterface(c_interface.ForTestingInterface):
    @legacy_function
    def do_sleep():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def return_error():
        function = LegacyFunctionSpecification()
        function.addParameter('out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def echo_2_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN, unit=units.m)
        function.addParameter('int_in2', dtype='int32', direction=function.IN, default=1, unit=units.kg)
        function.addParameter('int_out1', dtype='int32', direction=function.OUT, unit=units.m)
        function.addParameter('int_out2', dtype='int32', direction=function.OUT, unit=units.kg)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_x():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float32', direction=function.OUT, unit=units.m)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def set_x():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float32', direction=function.IN, unit=units.m)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def dummy():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
