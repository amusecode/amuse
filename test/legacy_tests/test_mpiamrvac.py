from amuse.legacy import *
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.mpiamrvac.interface import MpiAmrVacInterface
from amuse.legacy.mpiamrvac.interface import MpiAmrVac

class MpiAmrVacInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = MpiAmrVacInterface()
        instance.initialize_code()
        instance.stop()
    
