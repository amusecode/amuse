from amuse.test.amusetest import TestWithMPI

from .interface import arepoInterface
from .interface import arepo

class arepoInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = arepoInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
