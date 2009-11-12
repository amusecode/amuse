import unittest

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import attributes

class TestScalarAttributeDefinition(unittest.TestCase):
    
    def test1(self):
        definition = attributes.ScalarAttributeDefinition_Next(
            "set_test",
            "get_test",
            "test_argument",
            "test_attribute",
            "for testing purposes only", 
            units.m, 
            10 | units.m)
        
        class TestModule(object):
            def __init__(self):
                self.values = [1.0, 2.0, 3.0]
            def get_test(self, ids):
                return ([self.values[x] for x in ids],0)
            def set_test(self, ids, values):
                for index, value in zip(ids, values):
                    self.values[index] = value
                return 0
                
        o = TestModule()
        values = definition.get_values(o, [0,2])
        self.assertEquals(1.0 | units.m, values[0])
        self.assertEquals(3.0 | units.m, values[1])
        
        definition.set_values(o, [2,0], [1.0 | units.km, 5.0 | units.m])
        
        self.assertEquals(5, o.values[0])
        self.assertEquals(1000, o.values[2])
        
