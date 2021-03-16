from amuse.test.amusetest import TestWithMPI

from .interface import gpuhermite8Interface
from .interface import gpuhermite8

class gpuhermite8InterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = gpuhermite8Interface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
