from amuse.community import legacy_function
from amuse.rfi.core import CodeInterface, LegacyFunctionSpecification
from amuse.support.interface import InCodeComponentImplementation


class {interface_class}(CodeInterface):
    """Description of the interface with the community code.

    This class describes a set of functions in C++ or Fortran that are part of the AMUSE
    wrapper of the community code. These functions in turn will call functions in the
    community code, or update variables in it directly.

    These functions will be called by the worker program when it receives the
    corresponding request from the Python part of AMUSE.
    """

    {includes}

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="{code}_worker", **keyword_arguments)

    # here you must specify the prototypes of the interface functions:

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN, unit=None)
        function.addParameter('int_out', dtype='int32', direction=function.OUT, unit=None)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    # optionally, this can be shortened to:

    # @remote_function(can_handle_array=True)
    # def echo_int(int_in='i'):
    #     returns (int_out='i')


class {user_class}(InCodeComponentImplementation):
    """One line description of this code

    Some more details about what it does, any special features it has beyond the
    standard interfaces, and anything else the user needs to know.
    """
    def __init__(self, **options):
        """Create a {user_class} instance to run simulations with."""
        super().__init__(self,  {interface_class}(**options), **options)

    # the following alternative __init__ is appropiate for codes that use an unspecified
    # unit system (i.e. the quantities have dimension but no definite scale)
    #
    # def __init__(self, unit_converter=None, **options):
    #     self.unit_converter = unit_converter
    #     super().__init__(self,  {interface_class}(**options), **options)
    #
    # in this case you also need to use the define_converter below

    # typically the high level specification also contains the following:

    def define_state(self, handler):
        """Define the state model of the code."""
        # for example:
        # handler.set_initial_state('UNINITIALIZED')
        # handler.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        # handler.add_transition('END', 'STOPPED', 'stop', False)
        # handler.add_transition(
        #     'UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        # handler.add_method('STOPPED', 'stop')
        pass

    def define_properties(self, handler):
        # handler.add_property('name_of_the_getter', public_name="name_of_the_property")
        pass

    def define_parameters(self, handler):
        """Define model parameters.

        These have a native function for getting their value, another one for setting,
        and a name, description and default value. Functions with the appropriate names
        must be defined in the native wrapper code.
        """
        # handler.add_method_parameter(
        #     "name_of_the_getter",
        #     "name_of_the_setter",
        #     "parameter_name",
        #     "description",
        #     default_value = <default value>
        # )
        pass

    def define_particle_sets(self, handler):
        """Define any particle sets inside the model."""
        # handler.define_set('particles', 'index_of_the_particle')
        # handler.set_new('particles', 'new_particle')
        # handler.set_delete('particles', 'delete_particle')
        # handler.add_setter('particles', 'set_state')
        # handler.add_getter('particles', 'get_state')
        # handler.add_setter('particles', 'set_mass')
        # handler.add_getter('particles', 'get_mass', names=('mass',))
        pass

    def define_grids(self, handler):
        """Define any grids inside the model."""
        # handler.define_grid('grid',axes_names = ["x", "y"], grid_class=StructuredGrid)
        # handler.set_grid_range('grid', '_grid_range')
        # handler.add_getter('grid', 'get_grid_position', names=["x", "y"])
        # handler.add_getter('grid', 'get_rho', names=["density"])
        # handler.add_setter('grid', 'set_rho', names=["density"])
        pass

    # def define_converter(self, handler):
        # """Handle unit conversion if an (optional) unit converter is specified."""
        #     if self.unit_converter is not None:
        #         handler.set_converter(
        #             self.unit_converter.as_converter_from_si_to_generic()
        #         )
