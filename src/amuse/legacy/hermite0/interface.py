from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class Hermite(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """ 
        .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """		
    include_headers = ['hermite_code.h', 'parameters.h', 'local.h']

    def __init__(self):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")
        LiteratureRefs.__init__(self)

