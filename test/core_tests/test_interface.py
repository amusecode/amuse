from amuse.support import interface
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.binding import *
from amuse.support.data.parameters    import *

import unittest

class CodeInterfaceWithConvertedUnitsTests(unittest.TestCase):
    class TestClass(interface.CodeInterface):
        parameter_definitions = [
            ModuleAttributeParameterDefinition(
                "eps",
                "epsilon_squared", 
                "smoothing parameter for gravity calculations", 
                nbody_system.length * nbody_system.length, 
                0.3 | nbody_system.length * nbody_system.length
            )
        ]
        
        total_mass = CodeProperty("get_mass", nbody_system.mass)
        masses = (2.0 | nbody_system.mass, 3.0 | nbody_system.mass)
        
        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + (10.0 | nbody_system.length)
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithNBodyUnitsConverted(
                original,  
                convert_nbody
        )
        
        self.assertEquals(original.total_mass, 10.0 | nbody_system.mass)
        
        self.assertAlmostEquals(instance.total_mass.value_in(units.kg), 100.0, 10)
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithNBodyUnitsConverted(
                original,  
                convert_nbody
        )
        
        self.assertEquals(original.total_mass, 10.0 | nbody_system.mass)
        
        self.assertAlmostEquals(instance.total_mass.value_in(units.kg), 100.0, 10)
        
        self.assertEquals(original.add_to_length(5|nbody_system.length), 15.0 | nbody_system.length)
        self.assertAlmostEquals(instance.add_to_length(5|units.m).value_in(units.m), 55.0, 10)
        
        
        
    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithNBodyUnitsConverted(
                original,  
                convert_nbody
        )
        
        original.parameters.epsilon_squared =  2.0 | nbody_system.length ** 2
        
        self.assertAlmostEquals(original.eps,  2.0)
        
        instance.parameters.epsilon_squared = 100.0 | units.m ** 2
        
        self.assertAlmostEquals(instance.parameters.epsilon_squared.value_in(units.m**2),  100.0, 6)
        self.assertAlmostEquals(original.parameters.epsilon_squared.value_in(nbody_system.length ** 2),  4.0, 6)
        self.assertAlmostEquals(original.eps,  4.0, 6)
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterfaceWithNBodyUnitsConverted(
                original,  
                convert_nbody
        )
        
        self.assertEquals(original.masses[0], 2.0 | nbody_system.mass)
        
        masses = list(instance.masses)
        
        self.assertAlmostEquals(masses[0].value_in(units.kg), 20.0, 10)
        self.assertAlmostEquals(masses[1].value_in(units.kg), 30.0, 10)
        
        
        

