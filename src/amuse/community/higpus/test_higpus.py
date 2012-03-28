from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import HiGPUsInterface
from .interface import HiGPUs

class HiGPUsInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = HiGPUsInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
