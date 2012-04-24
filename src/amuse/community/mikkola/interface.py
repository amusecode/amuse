from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from numpy import pi

class MikkolaInterface(CodeInterface,
                       GravitationalDynamicsInterface):
    
#    include_headers = ['worker_code.h']
    use_modules = ['Mikkola',]
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="mikkola_worker", **keyword_arguments)
    
class Mikkola(GravitationalDynamics):

    def __init__(self, convert_nbody=None, **options):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        GravitationalDynamics.__init__(self,  MikkolaInterface(),
            convert_nbody,
            **options
        )

