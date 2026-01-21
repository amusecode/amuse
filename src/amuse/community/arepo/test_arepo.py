from amuse.test.amusetest import TestWithMPI

from .interface import ArepoInterface
from .interface import Arepo

class ArepoInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = ArepoInterface()
        result,error = instance.echo_int(12)  # TODO: Update test and add more...
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
