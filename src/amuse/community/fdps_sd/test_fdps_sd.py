from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import fdps_sdInterface
from .interface import fdps_sd

class fdps_sdInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = fdps_sdInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
