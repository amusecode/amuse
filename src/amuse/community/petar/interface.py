from amuse.rfi.core import legacy_function, LegacyFunctionSpecification
from amuse.community import (
    CodeInterface,
    LiteratureReferencesMixIn,
    StoppingConditionInterface,
    StoppingConditions,
)
from amuse.community.interface.gd import (
    GravitationalDynamics,
    GravitationalDynamicsInterface,
    GravityFieldInterface,
    GravityFieldCode,
)
from amuse.units import nbody_system


class petarInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    GravitationalDynamicsInterface,
    StoppingConditionInterface,
    GravityFieldInterface
):

    """
    Parallel, Particle-Particle & Particle-Tree & Few-body integration module

    .. [#] Namekata D., et al., 2018, PASJ, 70, 70
    .. [#] Iwasawa M., Tanikawa A., Hosono N., Nitadori K., Muranushi T., Makino J., 2016, PASJ, 68, 54
    .. [#] Iwasawa M., Portegies Zwart S., Makino J., 2015, ComAC, 2, 6
    .. [#] Wang, L., Nitadori, K., Makino, J., 2020, MNRAS, 493, 3398
    .. [#] Wang, L., Iwasawa, M., Nitadori, K., Makino, J., 2020, MNRAS, accepted, ArXiv: 2006.16560 [astro-ph]
    """

    include_headers = ['interface.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="petar_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def get_theta():
        """
        Get theta, the opening angle for building the tree
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'theta_for_tree', dtype='float64', direction=function.OUT,
            description=(
                "theta, the opening angle for building the tree:"
                " between 0 and 1"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_theta():
        """
        Set theta, the opening angle for building the tree
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'theta', dtype='float64', direction=function.IN,
            description=(
                "theta, the opening angle for building the tree:"
                " between 0 and 1"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    #@legacy_function
    #def get_gravitational_constant():
    #    """
    #    Get the value of gravitational constant
    #    """
    #    function = LegacyFunctionSpecification()
    #    function.addParameter(
    #        'gravitational_constant', dtype='float64', direction=function.OUT,
    #        description=(
    #            "gravitational constant:"
    #            " positive"
    #        )
    #    )
    #    function.result_type = 'int32'
    #    function.result_doc = """
    #    0 - OK
    #        the parameter was retrieved
    #    -1 - ERROR
    #        could not retrieve parameter
    #    """
    #    return function
    # 
    #@legacy_function
    #def set_gravitational_constant():
    #    """
    #    Set gravitational constant
    #    """
    #    function = LegacyFunctionSpecification()
    #    function.addParameter(
    #        'gravitational_constant', dtype='float64', direction=function.IN,
    #        description=(
    #            "gravitational constant:"
    #            " between 0 and 1"
    #        )
    #    )
    #    function.result_type = 'int32'
    #    function.result_doc = """
    #    0 - OK
    #        the parameter was set
    #    -1 - ERROR
    #        could not set parameter
    #    """
    #    return function

    @legacy_function
    def get_changeover_rout():
        """
        Get r_out, the changeover radius outer boundary reference for
        switching long/short-range interactions (if zero, auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_out', dtype='float64', direction=function.OUT,
            description=(
                "r_out, the changeover radius outer boundary reference"
                " for switching long/short-range interactions (if zero,"
                " auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_changeover_rout():
        """
        Set r_out, the changeover radius outer boundary reference for
        switching long/short-range interactions (if zero, auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_out', dtype='float64', direction=function.IN,
            description=(
                "r_out, the changeover radius outer boundary reference"
                " for switching long/short-range interactions (if zero,"
                " auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_changeover_ratio():
        """
        Get ratio_r_cut, the changeover radius ratio between the inner and
        outer boundaries
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ratio_r_cut', dtype='float64', direction=function.OUT,
            description=(
                "ratio_r_cut, the changeover radius ratio between the inner"
                " and outer boundaries"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_changeover_ratio():
        """
        Set ratio_r_cut, the changeover radius ratio between the inner and
        outer boundaries
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'ratio_r_cut', dtype='float64', direction=function.IN,
            description=(
                "ratio_r_cut, the changeover radius ratio between the inner"
                " and outer boundaries"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_group_radius():
        """
        Get r_bin, the group detection maximum radius to switch on AR
        (if zero, auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_bin', dtype='float64', direction=function.OUT,
            description=(
                "r_bin, the group detection maximum radius to switch on AR"
                " (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_group_radius():
        """
        Set r_bin, the group detection maximum radius to switch on AR
        (if zero, auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_bin', dtype='float64', direction=function.IN,
            description=(
                "r_bin, the group detection maximum radius to switch on AR"
                " (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_rsearch_min():
        """
        Get r_search_min, the minimum neighbor searching radius (if zero,
        auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_search_min', dtype='float64', direction=function.OUT,
            description=(
                "r_search_min, the minimum neighbor searching radius"
                " (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_rsearch_min():
        """
        Set r_search_min, the minimum neighbor searching radius (if zero,
        auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_search_min', dtype='float64', direction=function.IN,
            description=(
                "r_search_min, the minimum neighbor searching radius"
                " (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_tree_step():
        """
        Get dt_soft, the tree time step (if zero, auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dt_soft', dtype='float64', direction=function.OUT,
            description=(
                "dt_soft, the tree time step (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def set_tree_step():
        """
        Set dt_soft, the minimum neighbor searching radius (if zero,
        auto-determine)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'dt_soft', dtype='float64', direction=function.IN,
            description=(
                "dt_soft, the tree time step (if zero, auto-determine)"
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function


class petar(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody=None, **keyword_arguments):
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            petarInterface(**keyword_arguments),
            convert_nbody,
            **keyword_arguments)

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
        self.stopping_conditions.define_state(handler)

    def define_parameters(self, handler):
        GravitationalDynamics.define_parameters(self, handler)
        self.stopping_conditions.define_parameters(handler)

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing  parameter for gravity calculations",
            default_value=0.0 | nbody_system.length * nbody_system.length
        )

        handler.add_method_parameter(
            "get_changeover_rout",
            "set_changeover_rout",
            "r_out",
            ("changeover radius outer boundary reference for switching"
             " long/short-range interactions (if zero, auto-determine)"),
            default_value=0.0 | nbody_system.length
        )

        handler.add_method_parameter(
            "get_changeover_ratio",
            "set_changeover_ratio",
            "ratio_r_cut",
            "Changeover radius ratio between the inner and outer boundaries",
            default_value=0.1
        )

        handler.add_method_parameter(
            "get_group_radius",
            "set_group_radius",
            "r_bin",
            ("Group detection maximum radius to switch on AR (if zero,"
             " auto-determine)"),
            default_value=0.0 | nbody_system.length
        )

        handler.add_method_parameter(
            "get_rsearch_min",
            "set_rsearch_min",
            "r_search_min",
            "Minimum neighbor searching radius (if zero, auto-determine)",
            default_value=0.0 | nbody_system.length
        )

        handler.add_method_parameter(
            "get_theta",
            "set_theta",
            "theta",
            "Tree opening angle",
            default_value=0.3
        )

        #handler.add_method_parameter(
        #    "get_gravitational_constant",
        #    "set_gravitational_constant",
        #    "gravitational_constant",
        #    "Gravitational constant",
        #    default_value=0.3
        #)

        handler.add_method_parameter(
            "get_tree_step",
            "set_tree_step",
            "dt_soft",
            "Tree time step (if zero, auto-determine)",
            default_value=0.0 | nbody_system.time
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        self.stopping_conditions.define_methods(handler)

        handler.add_method(
            "set_eps2",
            (
                nbody_system.length * nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_eps2",
            (),
            (
                nbody_system.length * nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_changeover_rout",
            (
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_changeover_rout",
            (),
            (
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_group_radius",
            (
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_group_radius",
            (),
            (
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_rsearch_min",
            (
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_rsearch_min",
            (),
            (
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_tree_step",
            (
                nbody_system.time,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_tree_step",
            (),
            (
                nbody_system.time,
                handler.ERROR_CODE,
            )
        )

    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)


PetarInterface = petarInterface
Petar = petar
