import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core
from amuse.legacy.support.core import RemoteFunction, legacy_global
from muse_dynamics_mpi import PhiGRAPE as PhiGRAPE_basic


class PhiGRAPE(PhiGRAPE_basic):            
    def __init__(self, convert_nbody = None,nworker = None):
        if nworker is None: nworker=1
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker_gl')
        
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, nworker)
        self.channel = core.MpiChannel(self.intercomm)
        self.convert_nbody = convert_nbody
        
    @core.legacy_function
    def start_viewer():
        function = RemoteFunction()  
        function.id = 10
        return function

