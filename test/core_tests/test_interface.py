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
        
        

class CodeInterfaceWithCodeMethodsTests(unittest.TestCase):
    class TestClass(interface.CodeInterface):
       
        def _add_10_to_length(self, length):
            return length + 10
            
        add_10_to_length = CodeMethod('_add_10_to_length', (units.m,), units.m)
        
    def test1(self):
        instance = self.TestClass()
        
        self.assertEquals(20.0, instance._add_10_to_length(10.0))
        
        self.assertEquals(20.0 | units.m, instance.add_10_to_length(10.0 | units.m))
        self.assertEquals(1010.0 | units.m, instance.add_10_to_length(1.0| units.km))



class CodeInterfaceWithStateEngineTests(unittest.TestCase):
    class TestClass(interface.CodeInterface):
        def __init__(self):
            self.state = 0
            
        def always_works(self):
            return True
        
        def move_to_state_1(self):
            self.state = 1
            
        def move_to_state_2(self):
            self.state = 2
            
        def move_to_state_3(self):
            self.state = 3
            
        def move_to_state_4(self):
            self.state = 4
            
        def returns_1(self):
            return self.state
            
        def returns_2(self):
            return self.state
            
        def returns_3(self):
            return self.state
            
        def returns_4(self):
            return self.state
            
    def test1(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        self.assertTrue(instance.always_works())
        instance.move_to_state_1()
        self.assertEquals(1, instance.returns_1())
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_method('ONE', 'returns_1')
        instance.set_initial_state('ZERO')
        
        
        self.assertTrue(instance.always_works())
        self.assertEquals(instance._current_state.name, 'ZERO')
        instance.move_to_state_1()
        
        self.assertEquals(instance._current_state.name, 'ONE')
        self.assertEquals(instance.returns_1(), 1)
        
    def test3(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_method('ONE', 'returns_1')
        instance.set_initial_state('ZERO')
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_1(), 1)        
        self.assertEquals(instance._current_state.name, 'ONE')
        
    def test4(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_transition('ONE', 'TWO', 'move_to_state_2')
        instance.add_method('ONE', 'returns_1')
        instance.add_method('TWO', 'returns_2')
        instance.set_initial_state('ZERO')
        
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_2(), 2)        
        self.assertEquals(instance._current_state.name, 'TWO')
        
    def test5(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_transition('ONE', 'TWO', 'move_to_state_2')
        instance.add_transition('TWO', 'THREE', 'move_to_state_3')
        instance.add_transition('TWO', 'FOUR', 'move_to_state_4')
        instance.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        instance.add_method('ONE', 'returns_1')
        instance.add_method('TWO', 'returns_2')
        instance.add_method('THREE', 'returns_3')
        instance.add_method('FOUR', 'returns_4')
        instance.set_initial_state('ZERO')
        
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_4(), 4)        
        self.assertEquals(instance._current_state.name, 'FOUR')
        
        
        instance.set_initial_state('ZERO')
        self.assertEquals(instance._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_3(), 3)        
        self.assertEquals(instance._current_state.name, 'THREE')
        self.assertEquals(instance.returns_4(), 4)        
        self.assertEquals(instance._current_state.name, 'FOUR')
        
    
    def test6(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_transition('ONE', 'TWO', 'move_to_state_2')
        instance.add_transition('TWO', 'THREE', 'move_to_state_3')
        instance.add_transition('TWO', 'FOUR', 'move_to_state_4')
        instance.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        instance.add_method('ONE', 'returns_1')
        instance.add_method('TWO', 'returns_2')
        instance.add_method('THREE', 'returns_3')
        instance.add_method('FOUR', 'returns_4')
        instance.set_initial_state('ZERO')
        instance.do_automatic_state_transitions(False)
        
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        try:
            self.assertEquals(instance.returns_4(), 4)
            self.fail("Automatic state transitions is OFF, this method should fail")  
        except Exception, ex:
            print ex
            
            self.assertEquals(len(ex.transitions), 3)
            
            for x in ex.transitions:
                x.do()
            
            self.assertEquals(instance.returns_4(), 4)
            
    
    def test7(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_transition('ONE', 'TWO', 'move_to_state_2')
        
        instance.add_method('ONE', 'returns_1')
        instance.add_method('TWO', 'returns_2')
        instance.add_method('THREE', 'returns_3')
        instance.set_initial_state('ZERO')
        
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_2(), 2)
        try:
            self.assertEquals(instance.returns_3(), 3)
            self.fail("No transition to state 3 possible, function should error")  
        except Exception, ex:
            print ex
        
    def test8(self):
        original = self.TestClass()
        
        instance = interface.CodeInterfaceWithStateEngine(original)
        
        instance.add_transition('ZERO', 'ONE', 'move_to_state_1')
        instance.add_transition('ONE', 'TWO', 'move_to_state_2')
        instance.add_transition('THREE', 'ONE', 'move_to_state_1')
        instance.add_transition('TWO', 'THREE', 'move_to_state_1')
        
        
        instance.add_method('ONE', 'returns_1')
        instance.add_method('TWO', 'returns_2')
        instance.add_method('THREE', 'returns_3')
        instance.set_initial_state('ZERO')
        
        
        self.assertEquals(instance._current_state.name, 'ZERO')
        instance.move_to_state_1()
        self.assertEquals(instance._current_state.name, 'ONE')
        instance.move_to_state_2()
        self.assertEquals(instance._current_state.name, 'TWO')
        instance.move_to_state_1()
        self.assertEquals(instance._current_state.name, 'THREE')
        instance.move_to_state_1()
        self.assertEquals(instance._current_state.name, 'ONE')
        
        self.assertEquals(instance.returns_3(), 1)
        
