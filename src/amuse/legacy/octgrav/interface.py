from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class OctgravInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface):
    """
        .. [#] Gaburov, Nitadori, Harfst, Portegies Zwart & Makino,"A gravitational tree code on graphics processing units:
               Implementation in CUDA", in preparetion; and main MUSE paper, arXiv/0807.1996
    """

    include_headers = ['octgrav_code.h', 'parameters.h', 'worker_code.h', 'local.h']

    def __init__(self, convert_nbody = None, **kwargs):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code", **kwargs)
        """
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()

        self.convert_nbody = convert_nbody
        """
        LiteratureRefs.__init__(self)

    def setup_module(self):
        self.initialize_code()
        self.commit_parameters()
        self.commit_particles()

    def cleanup_module(self):
        self.cleanup_code()

class Octgrav(GravitationalDynamics):

    def __init__(self, convert_nbody = None):

        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()


        legacy_interface = BHTreeInterface()

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
        )
