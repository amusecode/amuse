from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class Hermite(LegacyInterface, LiteratureRefs, GravitationalDynamics):
		
    def __init__(self):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")
        