import unittest
import numpy

from amuse.support.data import core
from amuse.support.units import units

class TestMeasurement(unittest.TestCase):
    def test1(self):
        
        instance = core.Measurement(10 | units.s,
             ['position','mass'] , [units.m, units.kg], [1,2,3], 
             numpy.zeros((2, 3)))
        self.assertEquals(str(instance) , 'id\tposition\tmass\n-\tm\tkg\n========\t========\t========\n1\t0.0\t0.0\n2\t0.0\t0.0\n3\t0.0\t0.0' )
    
    def test2(self):
        
        instance = core.Measurement(10 | units.s,
             ['position','mass'] , [units.m, units.kg], [1,2,3], 
             numpy.array([[4,5,6],[7,8,9]], 'i'))
        self.assertEquals(str(instance) , 'id\tposition\tmass\n-\tm\tkg\n========\t========\t========\n1\t4\t7\n2\t5\t8\n3\t6\t9' )
        self.assertEquals(str(instance.get_value('mass',2)) , '8 kg')
        self.assertEquals(str(instance.get_value('position',3)) , '6 m')
