from amuse.support import interface
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.binding import *

import unittest

class CodeInterfaceWithConvertedUnitsTests(unittest.TestCase):
    class TestClass(interface.CodeInterface):
        parameter_definitions = [
        
        ]
        
        total_mass = CodeProperty("get_mass", nbody_system.mass)
        
        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + (10.0 | nbody_system.length)
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithConvertedUnits(
                original,  
                convert_nbody.as_converter_from_si_to_nbody()
        )
        
        self.assertEquals(original.total_mass, 10.0 | nbody_system.mass)
        
        self.assertAlmostEquals(instance.total_mass.value_in(units.kg), 100.0, 10)
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithConvertedUnits(
                original,  
                convert_nbody.as_converter_from_si_to_nbody()
        )
        
        self.assertEquals(original.total_mass, 10.0 | nbody_system.mass)
        
        self.assertAlmostEquals(instance.total_mass.value_in(units.kg), 100.0, 10)
        
        self.assertEquals(original.add_to_length(5|nbody_system.length), 15.0 | nbody_system.length)
        self.assertAlmostEquals(instance.add_to_length(5|units.m).value_in(units.m), 55.0, 10)
        
        

