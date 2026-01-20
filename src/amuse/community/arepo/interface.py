from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function
from amuse.community import LiteratureReferencesMixIn

from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics


class ArepoInterface(
    CodeInterface,
    GravitationalDynamicsInterface,
    LiteratureReferencesMixIn
):
    """
    Arepo is a cosmological magnetohydrodynamical moving-mesh simulation code,
    descended from GADGET.

    References:
        .. [#] Springel, V., 2010, MNRAS, 401, 791 (Arepo) [2010MNRAS.401..791S]
        .. [#] Pakmor, R., Bauer, A., Springel, V., 2011, MNRAS, 418, 1392 (Magnetohydrodynamics Module) [2011MNRAS.418.1392P]
        .. [#] Pakmor, R. et al., 2016, MNRAS, 455, 1134 (Gradient Estimation) [2016MNRAS.455.1134P]
        .. [#] Weinberger, R., Springel, V., Pakmor, R., 2020, ApJS, 248, 32 (Public Code Release) [2020ApJS..248...32W]
    """

    include_headers = ["worker_code.h", "interface.h"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="arepo_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        # TODO: Determine whether need to inherit from CodeWithDataDirectories.

    # This function has been kept as a basic template for future functions.
    # @legacy_function
    # def set_parameters():
    #     function = LegacyFunctionSpecification()
    #     function.addParameter("param_file", dtype="string", direction=function.IN)
    #     function.result_type = "int32"
    #     return function


class Arepo(GravitationalDynamics):

    def __init__(self, **options):
        GravitationalDynamics.__init__(self, ArepoInterface(**options), **options)

    def initialize_code(self):
        result = self.overridden().initialize_code()

        return result

    def define_methods(self, builder):
        # TODO: Determine how to link this to Arepo's run() - the main simulation loop.
        builder.add_method(
            "run_sim",
            (),
            (builder.ERROR_CODE)
        )
        # When simulation is finished, shutdown HDF5 & MPI, and exit(0)
        builder.add_method(
            "cleanup_code",
            (),
            (builder.ERROR_CODE)
        )
