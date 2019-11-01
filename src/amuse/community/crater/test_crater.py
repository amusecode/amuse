from amuse.test.amusetest import TestWithMPI
from amuse.lab import units

from interface import CraterInterface
from interface import Crater

class CraterInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = CraterInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()

    def test1(self):
        instance = Crater()
        instance.set_target_density(3.0|units.g/units.cm**3)
        rho = instance.get_target_density()
        self.assertAlmostRelativeEqual(rho, 3.0 | units.g/units.cm**3, 8)


        
        
