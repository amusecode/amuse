from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.higpus.interface import HiGPUsInterface
from amuse.community.higpus.interface import HiGPUs

class HiGPUsInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
