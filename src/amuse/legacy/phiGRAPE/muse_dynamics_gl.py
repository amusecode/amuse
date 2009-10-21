import numpy

from amuse.legacy import *
from muse_dynamics_mpi import PhiGRAPE as PhiGRAPE_basic


class PhiGRAPE(PhiGRAPE_basic):            
    def __init__(self, convert_nbody = None,nworker = 1):
        LegacyInterface.__init__(self,name_of_the_worker = 'muse_worker_gl', \
          number_of_workers = nworker)
        self.convert_nbody = convert_nbody
        
    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()  
        return function

