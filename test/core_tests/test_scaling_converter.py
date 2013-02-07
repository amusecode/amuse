from amuse.units import scaling_converter
from amuse.units import nbody_system

from amuse.test import amusetest

class TestScalingConverter(amusetest.TestCase):

    def test1(self):
        converter = scaling_converter.ScalingConverter(
            length = 0.2,
            time   = 0.1, 
        )
        input = 1 | nbody_system.time
        output = converter.convert(input)
        self.assertAlmostRelativeEquals(output, 0.1 | nbody_system.time)
    
        
    def test2(self):
        converter = scaling_converter.ScalingConverter(
            length = 0.2,
            time   = 0.1, 
        )
        input = 1 | nbody_system.length ** 2
        output = converter.convert(input)
        self.assertAlmostRelativeEquals(output, 0.2 * 0.2 | nbody_system.length ** 2)