from amuse.community.interface.stopping_conditions import StoppingConditionInterface

from amuse.rfi.core import CodeInterface, legacy_function, LegacyFunctionSpecification
from amuse.community import NO_UNIT


class ForTestingInterface(CodeInterface, StoppingConditionInterface):
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    include_headers = ['c_interface.h']

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def reset_stopping_conditions():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def next_index_for_stopping_condition():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_unit = NO_UNIT
        function.can_handle_array = False
        return function

    @legacy_function
    def set_stopping_condition_info():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def set_stopping_condition_particle_index():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def mpi_setup_stopping_conditions():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def mpi_collect_stopping_conditions():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function

    @legacy_function
    def mpi_distribute_stopping_conditions():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function


class ForTestingInterfaceFortranModule(ForTestingInterface):
    use_modules = ['StoppingConditions', 'AmuseInterface']

    include_headers = ['c_ext_interface.h']

    @legacy_function
    def fire_condition():
        function = LegacyFunctionSpecification()
        function.addParameter('condition_to_set', dtype='int32', direction=function.IN)
        function.addParameter('particle_index_1', dtype='int32', direction=function.IN)
        function.addParameter('particle_index_2', dtype='int32', direction=function.IN)
        function.addParameter('rank', dtype='int32', direction=function.IN, default=-1)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
