from amuse.support import interface
from amuse.support.exceptions import AmuseException


from amuse.datamodel.binding import *
from amuse.datamodel.parameters import *
from amuse.support.core import OrderedDictionary
from amuse.support import exceptions

from amuse.test import amusetest
import numpy
import pickle
from amuse.units import units
from amuse.units import nbody_system
from amuse import datamodel
from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification

class CodeInterfaceWithConvertedUnitsTests(amusetest.TestCase):
    class TestClass(object):
        
        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + 10.0
            
        def return_an_errorcode(self, error):
            return error
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (nbody_system.length,), nbody_system.length)
        
        self.assertAlmostEqual(instance.mass.value_in(units.kg), 100.0, 10)
        
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (nbody_system.length,), nbody_system.length)
                
        self.assertAlmostEqual(instance.add_to_length(5|units.m).value_in(units.m), 55.0, 10)
        

    def test3(self):
        original = self.TestClass()
        instance = interface.InCodeComponentImplementation(original)

        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add_10')
           
        self.assertFalse(instance.add_10.is_async_supported)
        
        
    def test4(self):
        original = self.TestClass()
        instance = interface.InCodeComponentImplementation(original)

        handler = instance.get_handler('METHOD')
        handler.add_method('return_an_errorcode', (handler.NO_UNIT,), handler.ERROR_CODE)
           
        self.assertEqual(instance.return_an_errorcode(0), None)
        self.assertRaises(Exception, lambda : instance.return_an_errorcode(-1),
            expected_message="Error when calling 'return_an_errorcode' of a 'InCodeComponentImplementation', errorcode is -1"
        )
        
        
class CodeInterfaceWithUnitsOnLegacyFunctionTests(amusetest.TestCase):
        
    def test1(self): 
        class TestImplementation(object):
        
            def get_mass(self):
                return 10.0, 0
                
        class TestInterface(object):
           
            @legacy_function
            def get_mass():
                function = LegacyFunctionSpecification()
                function.addParameter('input1', dtype='d', direction=function.OUT, unit=units.kg)
                function.result_type = 'i'
                return function
            
        original = TestInterface()
        instance = interface.InCodeComponentImplementation(original)
        
        instance.get_handler("LEGACY").legacy_interface = TestImplementation()
        self.assertAlmostRelativeEquals(instance.get_mass(), 10 | units.kg)
    
    def test2(self): 
        class TestImplementation(object):
        
            def echo_one(self, input):
                return input, 0
                
        class TestInterface(object):
           
            @legacy_function
            def echo_one():
                function = LegacyFunctionSpecification()
                function.addParameter('input', dtype='d', direction=function.IN, unit=units.kg)
                function.addParameter('output', dtype='d', direction=function.OUT, unit=units.g)
                function.result_type = 'i'
                return function
            
        original = TestInterface()
        instance = interface.InCodeComponentImplementation(original)
        
        instance.get_handler("LEGACY").legacy_interface = TestImplementation()
        self.assertAlmostRelativeEquals(instance.echo_one(1|units.kg), 1 | units.g)
    
    def test3(self): 
        class TestImplementation(object):
        
            def return_error(self):
                return  -1
                
        class TestInterface(object):
           
            @legacy_function
            def return_error():
                function = LegacyFunctionSpecification()
                function.result_type = 'i'
                return function
            
        original = TestInterface()
        instance = interface.InCodeComponentImplementation(original)
        
        instance.get_handler("LEGACY").legacy_interface = TestImplementation()
        self.assertRaises(Exception,instance.return_error, " Error when calling 'echo_one' of a 'InCodeComponentImplementation', errorcode is -1")

    def test4(self): 
        class TestImplementation(object):
        
            def echo_one(self, input):
                return input, 0
                
        class TestInterface(object):
           
            @legacy_function
            def echo_one():
                function = LegacyFunctionSpecification()
                function.addParameter('input', dtype='d', direction=function.IN, unit=units.deg)
                function.addParameter('output', dtype='d', direction=function.OUT, unit=units.deg)
                function.result_type = 'i'
                return function
            
        original = TestInterface()
        instance = interface.InCodeComponentImplementation(original)
        
        instance.get_handler("LEGACY").legacy_interface = TestImplementation()
        self.assertAlmostRelativeEquals(instance.echo_one(1. | units.rad), 1. | units.rad)
        # this was: self.assertAlmostRelativeEquals(instance.echo_one(1.), 1. | units.deg)
        # but after 85bd5d99 or 945bc46, it needs to be:
        self.assertAlmostRelativeEquals(instance.echo_one(1.), 1. | units.rad)
        # this is indeed proper behaviour!! 
        
class CodeInterfaceWithMethodsAndPropertiesTests(amusetest.TestCase):
    class TestClass(object):
       
        def add_10_to_length(self, length):
            return length + 10

        def get_one(self):
            return 1.0, 0.0
            
        def get_state(self, id):
            return (1.0, 2.0, 3.0, 0.0)
        
        def get_state_error(self, id):
            return (1.0, 2.0, 3.0, -1.0)
        
    def test1(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)

        handler = instance.get_handler('METHOD')
        handler.add_method('add_10_to_length', (units.m,), units.m)
        
        self.assertEqual(20.0, original.add_10_to_length(10.0))
        
        self.assertEqual(20.0 | units.m, instance.add_10_to_length(10.0 | units.m))
        self.assertEqual(1010.0 | units.m, instance.add_10_to_length(1.0| units.km))
        
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_10_to_length', (units.m,), units.m, public_name = 'add_10')
           
        self.assertEqual(20.0 | units.m, instance.add_10(10.0 | units.m))
        self.assertEqual(1010.0 | units.m, instance.add_10(1.0| units.km))
        
    def test3(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_one', units.m)
        
        self.assertEqual(1.0 | units.m, instance.one)
        
    def test4(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_one', units.m, public_name = 'get_one')
        
        self.assertEqual(1.0 | units.m, instance.get_one)
        
    
    def test5(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_state', (handler.NO_UNIT,), (units.m, units.m, units.kg, handler.ERROR_CODE))
        
        result = instance.get_state(1)
        self.assertEqual(3, len(result))
        self.assertEqual(1.0 | units.m, result[0])
        self.assertEqual(3.0 | units.kg, result[2])
        
    
    def test6(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_state_error', (handler.NO_UNIT,), (units.m, units.m, units.kg, handler.ERROR_CODE))
        
        self.assertRaises(AmuseException, instance.get_state_error, 1, 
            expected_message = "Error when calling 'get_state_error' of a 'InCodeComponentImplementation', errorcode is -1.0")
            

class ClassWithState(object):
    
    def __init__(self):
        self.state = 0
        self.number_of_times_move_to_state_1_called = 0
        self.number_of_times_move_to_state_2_called = 0
        self.number_of_times_move_to_state_3_called = 0
        self.number_of_times_move_to_state_4_called = 0
        
    def always_works(self):
        return True
    
    def move_to_state_1(self):
        self.number_of_times_move_to_state_1_called += 1
        self.state = 1
        
    def move_to_state_2(self):
        self.number_of_times_move_to_state_2_called += 1
        self.state = 2
        
    def move_to_state_3(self):
        self.number_of_times_move_to_state_3_called += 1
        self.state = 3
        
    def move_to_state_4(self):
        self.number_of_times_move_to_state_4_called += 1
        self.state = 4
        
    def returns_1(self):
        return self.state
        
    def returns_2(self):
        return self.state
        
    def returns_3(self):
        return self.state
        
    def returns_4(self):
        return self.state
class CodeInterfaceTests(amusetest.TestCase):
            
    def test1(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        self.assertTrue(instance.always_works())
        instance.move_to_state_1()
        self.assertEqual(1, instance.returns_1())
        
    def test2(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        handler = instance.get_handler('STATE')
        
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        self.assertTrue(instance.always_works())
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.move_to_state_1()
        
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        self.assertEqual(instance.returns_1(), 1)
        
    
    def test3(self):
        class TestCodeInterface(interface.InCodeComponentImplementation):
            
            def move_to_state_1(self):
                self.overridden().move_to_state_1()
                self.traced = True
                
        original = ClassWithState()
        
        instance = TestCodeInterface(original)
        instance.traced = False
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_1(), 1)
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        self.assertTrue(instance.traced)
        
        
    def test4(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_1(), 1)    
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        
    def test5(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_2(), 2)    
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
        
    def test6(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'THREE', 'move_to_state_3')
        handler.add_transition('TWO', 'FOUR', 'move_to_state_4')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.add_method('FOUR', 'returns_4')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_4(), 4)    
        self.assertEqual(instance.get_name_of_current_state(), 'FOUR')
        
        
        handler.set_initial_state('ZERO')
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_3(), 3)    
        self.assertEqual(instance.get_name_of_current_state(), 'THREE')
        self.assertEqual(instance.returns_4(), 4)    
        self.assertEqual(instance.get_name_of_current_state(), 'FOUR')
        
    
    def test7(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'THREE', 'move_to_state_3')
        handler.add_transition('TWO', 'FOUR', 'move_to_state_4')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.add_method('FOUR', 'returns_4')
        handler.set_initial_state('ZERO')
        handler.do_automatic_state_transitions(False)
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        try:
            self.assertEqual(instance.returns_4(), 4)
            self.fail("Automatic state transitions is OFF, this method should fail")
        except Exception as ex:
            print(ex)
            
            self.assertEqual(len(ex.transitions), 3)
            
            for x in ex.transitions:
                x.do()
            
            self.assertEqual(instance.returns_4(), 4)
            
    
    def test8(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_2(), 2)
        self.assertRaises(Exception, instance.returns_3, 
            expected_message = "While calling returns_3 of InCodeComponentImplementation: No transition from current state state 'TWO' to state 'THREE' possible")
        
    def test9(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        handler.add_transition('TWO', 'THREE', 'move_to_state_1')
        
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.move_to_state_1()
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        instance.move_to_state_2()
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
        instance.move_to_state_1()
        self.assertEqual(instance.get_name_of_current_state(), 'THREE')
        instance.move_to_state_1()
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        
        self.assertEqual(instance.returns_3(), 1)
        



    def test10(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        handler.add_transition('TWO', 'THREE', 'move_to_state_1')
        
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.set_initial_state('ZERO')
        
        self.assertTrue(instance.state_machine.is_enabled)
        self.assertEqual(instance.state_machine.get_name_of_current_state(), 'ZERO')
        instance.move_to_state_1()
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        instance.state_machine.disable()
        self.assertFalse(instance.state_machine.is_enabled)
        instance.move_to_state_2()
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        
    def test11(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        handler = instance.get_handler('STATE')
        
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('!ZERO', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'ONE', 'move_to_state_1')
        handler.add_method('TWO', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.returns_1()
        
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
        self.assertEqual(instance.number_of_times_move_to_state_1_called, 1)
        self.assertEqual(instance.number_of_times_move_to_state_2_called, 1)
    
    
        
    def test12(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original, log_transitions = True)
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'THREE', 'move_to_state_3')
        handler.add_transition('THREE', 'FOUR', 'move_to_state_4')
        handler.add_transition('!ZERO!ONE!THREE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'ONE', 'move_to_state_1')
        handler.add_method('TWO', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.returns_1()
        
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
        self.assertEqual(instance.number_of_times_move_to_state_1_called, 1)
        self.assertEqual(instance.number_of_times_move_to_state_2_called, 1)
        self.assertEqual(instance.number_of_times_move_to_state_3_called, 1)
        self.assertEqual(instance.number_of_times_move_to_state_4_called, 1)
        
    
        
    def test13(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original, log_transitions = True)
        handler = instance.get_handler('STATE')
        handler.add_transition('!ZERO!ONE!TWO', 'FOUR', 'move_to_state_4')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'THREE', 'move_to_state_3')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('!ZERO', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.returns_1()
        
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        self.assertEqual(instance.number_of_times_move_to_state_1_called, 1)
        self.assertEqual(instance.number_of_times_move_to_state_2_called, 0)
        self.assertEqual(instance.number_of_times_move_to_state_3_called, 0)
        self.assertEqual(instance.number_of_times_move_to_state_4_called, 0)
        instance.move_to_state_2()
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
        instance.returns_1()
        
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
    
    def test14(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.set_initial_state('ZERO')
        handler.remove_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.remove_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('ZERO', 'TWO', 'move_to_state_1')
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        instance.returns_2()
        self.assertEqual(instance.get_name_of_current_state(), 'TWO')
    
    def test15(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'THREE', 'move_to_state_3')
        handler.add_transition('TWO', 'FOUR', 'move_to_state_4')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        handler.add_method('THREE', 'returns_1')
        handler.add_method('TWO', 'returns_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_1(), 1)    
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        
    
    def test9(self):
        original = ClassWithState()
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('TWO', 'THREE', 'move_to_state_3')
        handler.add_transition('TWO', 'FOUR', 'move_to_state_4')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        
        handler.add_method('THREE', 'returns_3')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertEqual(instance.returns_1(), 1)    
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')
        
        pickled_instance = pickle.dumps(instance)
        unpickled_instance = pickle.loads(pickled_instance)
        self.assertEqual(unpickled_instance.get_name_of_current_state(), 'ONE')
        self.assertEqual(unpickled_instance.returns_2(), 2)    
        self.assertEqual(unpickled_instance.get_name_of_current_state(), 'TWO')
        self.assertEqual(instance.get_name_of_current_state(), 'ONE')

    def test10(self):
        original = ClassWithState()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
# remove to generate exception
#        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.set_initial_state('ZERO')
        
        
        self.assertEqual(instance.get_name_of_current_state(), 'ZERO')
        self.assertRaises( Exception, instance.returns_2, expected_message=
         "While calling returns_2 of InCodeComponentImplementation: No transition from current state state 'ZERO' to state 'TWO' possible")

        
class CodeInterfaceWithUnitsAndStateTests(amusetest.TestCase):
    class TestClass(object):
       
        def __init__(self):
            self.value = 10.0
            
        def add_to_length(self, length):
            return length + self.value
        
        def move_to_20(self):
            self.value = 20
        
    def test1(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m)
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add_to_length')
        
        self.assertEqual(40.0 | units.m, instance.add_to_length(20.0 | units.m))
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add')
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add_to_length')
        
        self.assertEqual(40.0 | units.m, instance.add(20.0 | units.m))
        
    def test3(self):
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add')
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add')
        
        self.assertEqual(40.0 | units.m, instance.add(20.0 | units.m))
        

class CodeInterfaceWithErrorHandlingTests(amusetest.TestCase):
    class TestClass(object):
        errorcode = 0
        
        def get_mass(self):
            return 10.0, self.errorcode
            
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass', (), (units.m, handler.ERROR_CODE,))
        handler = instance.get_handler('ERRORCODE')
        handler.add_errorcode(-2, "no such method")
        handler.add_errorcode(-3, "not available")
        
        self.assertEqual(instance.get_mass(), 10.0 | units.m)
        original.errorcode = -2
        self.assertRaises(AmuseException, instance.get_mass, expected_message = 
            "Error when calling 'get_mass' of a 'InCodeComponentImplementation', errorcode is -2, error is 'no such method'")
            
        original.errorcode = -1
        self.assertRaises(AmuseException, instance.get_mass, expected_message = 
            "Error when calling 'get_mass' of a 'InCodeComponentImplementation', errorcode is -1")
            

class CodeInterfaceWithParticlesTests(amusetest.TestCase):
    class TestClass(object):

        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + 10.0
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        
class ParticlesWithBindingInterface(object):

    def __init__(self):
        self.masses = {}
        
    def get_mass(self, id):
        masses = []
        errors = []
        for x in id:
            masses.append(self.masses[x])
            errors.append(0)
        result = OrderedDictionary()
        result["mass"] = masses
        result["__result"] = errors
        return result
    
    def set_mass(self, id, mass):
        for i,m in zip(id,mass):
            self.masses[i] = m
            
        return ( [0] * len(id),)
        
    def new_particle(self, mass):
        ids = []
        errors = []
        
        for x in mass:
            id = len(self.masses)
            self.masses[len(self.masses)]  = x
            ids.append(id)
            errors.append(0)
            
        return (ids, errors)
    
    def delete_particle(self, id):
        errors = []
        for x in id:
            self.masses[x].stop()
            errors.append(0)
        return errors
        
    def get_number_of_particles(self):
        return (len(self.masses), 0)
        
    def get_colliding_indices(self):
        return (1,2,0)
        
    
    def get_next(self, id):
        return ([x+1 for x in id], [0 for x in id])
        
    def add_1_to_mass(self, id):
        if isinstance(id, (int, numpy.int64, numpy.int32)):
            self.masses[id] += 1.0
            return [0]
        for i in id:
            self.masses[i] += 1.0
        return [0 for x in id]
    

    def get_heaviest_particle(self):
        max = -1
        id = -1
        for index, mass in self.masses.items():
            if mass > max:
                id = index
                max = mass
        return id
        
    def get_number_of_lightest_particles(self):
        return 2
        
    def get_lightest_particles(self, offset):
        sorted_masses = sorted(list(self.masses.items()), key = lambda x : x[1])
        return [sorted_masses[x][0] for x in offset]
        
    def get_number_of_particles_with_same_mass(self, id):
        result = []
        for index in id:
            mass = self.masses[index]
            count = 0
        
            for otherindex, m in self.masses.items():
                if otherindex == index:
                    continue
                if abs(m - mass) < 1.0:
                    count += 1
            result.append(count)
        return result
        
    def get_particle_with_same_mass(self, id, offset):
        result = []
        for index, index_in_list in zip(id, offset):
            mass = self.masses[index]
            count = 0
            for otherindex, m in self.masses.items():
                if otherindex == index:
                    continue
                if abs(m - mass) < 1.0:
                    count += 1
                    if count == (index_in_list + 1):
                        result.append(otherindex)
                        break
        return result
class TestParticlesWithBinding(amusetest.TestCase):
    
            
            
    def test1(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        
        
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(original.masses), 4)
        self.assertEqual(original.masses[0], 3.0)
        self.assertEqual(original.masses[3], 6.0)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles[0].mass, 3.0 | units.kg)
        
    def test2(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method('get_colliding_indices',(), (handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE,))
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        
        
        handler.add_query('particles', 'get_colliding_indices', public_name = 'select_colliding')
        handler.add_select_from_particle('particles', 'get_next', public_name = 'next')
        
        
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        colliding_particles = instance.particles.select_colliding()
        
        self.assertEqual(len(colliding_particles), 2)
        self.assertEqual(colliding_particles.mass[0], 4.0 | units.kg)
        self.assertEqual(colliding_particles.mass[1], 5.0 | units.kg)
        
        attribute_names = dir(instance)
        
        self.assertTrue('particles' in attribute_names)
        self.assertTrue('get_colliding_indices' in attribute_names)
        
    def test3(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('get_next',(handler.INDEX,), (handler.LINK("particles"), handler.ERROR_CODE))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_getter('particles', 'get_next', names = ('next_particle',))
        
        
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles[0].next_particle, instance.particles[1])
        self.assertEqual(instance.particles[1].next_particle, instance.particles[2])
        self.assertEqual(instance.particles[2].next_particle, instance.particles[3])
        self.assertEqual(instance.particles[3].next_particle, None)
        
        
    def test4(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        
        handler.add_method('add_1_to_mass',(handler.INDEX,), (handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_method('particles', 'add_1_to_mass', 'add_one')
        
        
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles.mass, units.kg.new_quantity([3.0, 4.0, 5.0, 6.0]))
        instance.particles[0].add_one()
        self.assertEqual(instance.particles.mass, units.kg.new_quantity([4.0, 4.0, 5.0, 6.0]))
        instance.particles.add_one()
        self.assertEqual(instance.particles.mass, units.kg.new_quantity([5.0, 5.0, 6.0, 7.0]))
        
        
    def test5(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (nbody_system.mass, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, nbody_system.mass,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(nbody_system.mass,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        
        
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)

        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(original.masses), 4)
        self.assertAlmostEqual(original.masses[0], 0.3, 5)
        self.assertAlmostEqual(original.masses[3], 0.6, 5)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles[0].mass, 3.0 | units.kg)
        
        
        
        
    def test6(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (nbody_system.mass, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, nbody_system.mass,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(nbody_system.mass,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_subselect_in_set('particles', 'get_heaviest_particle', public_name = 'heaviest')
        handler.add_subselect_in_set('particles', 'get_lightest_particles', get_number_of_particles_name = 'get_number_of_lightest_particles', public_name = 'lightest')
        
        
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg, 5.0 | units.m )
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
    
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 6.0, 5.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles[0].mass, 3.0 | units.kg)
        subset = instance.particles.heaviest()
        self.assertEqual(len(subset), 1)
        self.assertEqual(subset[0], local_particles[2])
        
        subset = instance.particles.lightest()
        self.assertEqual(len(subset), 2)
        self.assertEqual(subset[0], local_particles[0])
        self.assertEqual(subset[1], local_particles[1])
    

    def test7(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (nbody_system.mass, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, nbody_system.mass,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(nbody_system.mass,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_subselect_from_particle('particles', 'get_particle_with_same_mass', get_number_of_particles_name = 'get_number_of_particles_with_same_mass', public_name = 'samemass')
        
        
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kg, 5.0 | units.m )
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
    
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.5, 4.0, 6.0, 4.5])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(instance.particles[1].mass, 4.0 | units.kg)
        subset = instance.particles[1].samemass()
        self.assertEqual(len(subset), 2)
        self.assertEqual(subset[0], local_particles[0])
        self.assertEqual(subset[1], local_particles[3])
    
    def test8(self):
        original = ParticlesWithBindingInterface()
        original.masses[1] = 10
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('get_next',(handler.INDEX,), (handler.LINK("particles"), handler.ERROR_CODE))
        pickled_handler = pickle.dumps(handler)
        unpickled_handler = pickle.loads(pickled_handler)
        method1 = handler.get_attribute('get_mass', original.get_mass)
        method2 = unpickled_handler.get_attribute('get_mass', unpickled_handler.interface.legacy_interface.get_mass)
        self.assertTrue(method1 is not None)
        self.assertTrue(method2 is not None)
        self.assertAlmostRelativeEquals(method1([1])[0], 10 | units.kg)
        self.assertAlmostRelativeEquals(method2([1])[0], 10 | units.kg)
        original.masses[1] = 20
        self.assertAlmostRelativeEquals(method1([1]), 20 | units.kg)
        self.assertAlmostRelativeEquals(method2([1]), 10 | units.kg)
        
    def test9(self):
        original = ParticlesWithBindingInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        
        local_particles = datamodel.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEqual(len(original.masses), 4)
        self.assertEqual(original.masses[0], 3.0)
        self.assertEqual(original.masses[3], 6.0)
        
        pickled_instance = pickle.dumps(instance)
        unpickled_instance = pickle.loads(pickled_instance)
       
        self.assertEqual(len(unpickled_instance.particles), 4)
        self.assertEqual(unpickled_instance.particles[1].mass, 4.0 | units.kg)
        
        
class TestGridWithBinding(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (11,5,5)
        
        def __init__(self):
            self.storage = numpy.arange(
                self.shape[0]*self.shape[1]*self.shape[2]).reshape(self.shape)
            
        def get_range(self):
            return (0,self.shape[0]-1,0,self.shape[1]-1,0,self.shape[2]-1)
        
        def get_position(self, i_s, j_s, k_s):
            return [numpy.asarray(i_s), numpy.asarray(j_s)]
        
        def get_a(self,i_s,j_s,k_s):
            return [numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)]),]
            
        def set_a(self, i_s, j_s, k_s, values):
            index = 0
            for i,j,k in zip(i_s, j_s, k_s):
                self.storage[i][j][k] = values[index].value_in(units.m)
                index += 1
        
            
        
            
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX,handler.INDEX, units.kg,), ())
              
        self.assertEqual(instance.get_a([1],[2],[3]), [38] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass',))
        handler.add_getter('grid', 'get_a', names = ('mass',))
        
        grid = instance.grid
        
        
        self.assertEqual(grid[1][2][3].mass, 38 | units.kg)
        self.assertEqual(grid[0:2][1][2][3].mass, 38 | units.kg)
        self.assertEqual(len(grid[1][2].mass), 5)
        self.assertTrue(numpy.all(grid[1][2].mass == [35,36,37,38,39] | units.kg))
        
    
            
    def test2(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_position',(handler.INDEX, handler.INDEX,handler.INDEX,), (units.m, units.m,))
      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid', axes_names = ['x','y'])
        handler.add_getter('grid', 'get_position', names = ('x','y',))
        
        grid = instance.grid
        
        
        self.assertEqual(grid[1][2][3].position, [1,2] |units.m)
        self.assertEqual(grid[1][2][1].position, [1,2] |units.m)
        self.assertEqual(grid[1][2][2].position, [1,2] |units.m)
        self.assertEqual(grid[1][2][0].position, [1,2] |units.m)
        self.assertEqual(grid[1][2].position, [[1,2],[1,2],[1,2],[1,2],[1,2]] |units.m)
        self.assertEqual(grid[0][1][1].position, [0,1] |units.m)

    def test3(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_position',(handler.INDEX, handler.INDEX,handler.INDEX,), (units.m, units.m,))
      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid', axes_names = ['x','y'])
        handler.add_getter('grid', 'get_position', names = ('x','y',))
        
        grid = instance.grid
        
        self.assertEqual(grid.__class__.__name__, "Grid")
        self.assertTrue(isinstance(grid, datamodel.RegularGrid))
        
    def test4(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_position',(handler.INDEX, handler.INDEX,handler.INDEX,), (units.m, units.m,))
      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid', axes_names = ['x','y'], grid_class=datamodel.CartesianGrid)
        handler.add_getter('grid', 'get_position', names = ('x','y',))
        
        grid = instance.grid
        
        self.assertEqual(grid.__class__.__name__, "CartesianGrid")
        

class TestGridWithBinding2(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (11,5)
        
        def __init__(self):
            self.storage = numpy.arange(
                self.shape[0]*self.shape[1]).reshape(self.shape)
            
        def get_range(self):
            return (0,self.shape[0]-1,0,self.shape[1]-1)
        
        def get_position(self, i_s):
            return [numpy.asarray(i_s)]
        
        def get_a(self,i_s,j_s):
            #print "indices:", i_s, j_s, k_s
            #print "values:", numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)])
            return [numpy.asarray([(self.storage[i][j]) for i,j in zip(i_s, j_s)]),]
            
        def set_a(self, i_s, j_s, values):
            index = 0
            for i,j in zip(i_s, j_s):
                self.storage[i][j] = values[index]
                index += 1
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX, units.kg,), ())
              
        self.assertEqual(instance.get_a([1],[2]), [7] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass',))
        handler.add_getter('grid', 'get_a', names = ('mass',))
        
        grid = instance.grid
  
        self.assertEqual(grid[1,2].mass, original.storage[1,2] | units.kg)
        self.assertEqual(len(grid[1].mass), len(original.storage[1]))
  
        grid[1,1:].mass=[5. | units.kg]*4

        self.assertEqual(original.storage[1,1:],5)
        self.assertEqual(grid[1,1:].mass, original.storage[1,1:] | units.kg)

class TestGridWithBinding3(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (11,5)
        
        def __init__(self):
            self.storage = numpy.arange(
                self.shape[0]*self.shape[1]).reshape(self.shape)
            
        def get_range(self):
            return (0,self.shape[0]-1)
        def get_range2(self):
            return (0,self.shape[1]-1)
                
        def get_a(self,i_s,j_s):
            return [numpy.asarray([(self.storage[i][j]) for i,j in zip(i_s, j_s)]),]
            
        def set_a(self, i_s, j_s, values):
            for i,j,val in zip(i_s, j_s,values):
                self.storage[i][j] = val
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX, units.kg,), ())
              
        self.assertEqual(instance.get_a([1],[2]), [7] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_gridded_setter('grid', 'set_a','get_range2', names = ('mass',))
        handler.add_gridded_getter('grid', 'get_a','get_range2', names = ('mass',))
        
        grid = instance.grid
  
        self.assertEqual(grid[1].mass, original.storage[1] | units.kg)
        self.assertEqual(len(grid[1].mass), len(original.storage[1]))

        self.assertEqual(grid[2:4].mass, original.storage[2:4] | units.kg)

        self.assertEqual(grid[-3:].mass, original.storage[-3:] | units.kg)
  
        grid[1].mass=[5.]*5 | units.kg

        self.assertEqual(original.storage[1],5)
        self.assertEqual(grid[1].mass, original.storage[1] | units.kg)

class TestGridWithBinding4(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (3,4,5)
        
        def __init__(self):
            self.storage = numpy.arange(
                self.shape[0]*self.shape[1]*self.shape[2]).reshape(self.shape)
            
        def get_range(self):
            return (0, self.shape[0]-1, 0, self.shape[1]-1)
        def get_range2(self):
            return (0, self.shape[2]-1)
                
        def get_a(self, i_s, j_s, k_s):
            return [numpy.asarray([self.storage[i,j,k] for i,j,k in zip(i_s, j_s, k_s)]),]
            
        def set_a(self, i_s, j_s, k_s, values):
            for i,j,k,val in zip(i_s, j_s,k_s,values):
                self.storage[i,j,k] = val
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX, handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX,handler.INDEX, units.kg,), ())
              
        self.assertEqual(instance.get_a([1],[2],[3]), [33] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_gridded_setter('grid', 'set_a','get_range2', names = ('mass',))
        handler.add_gridded_getter('grid', 'get_a','get_range2', names = ('mass',))
        
        grid = instance.grid
  
        self.assertEqual(grid[1,1].mass, original.storage[1,1] | units.kg)
        self.assertEqual(grid[0:2,1:3].mass, original.storage[0:2,1:3] | units.kg)
        self.assertEqual(len(grid[1].mass), len(original.storage[1]))
        self.assertEqual(grid[2:4].mass, original.storage[2:4] | units.kg)
        self.assertEqual(grid[-3:].mass, original.storage[-3:] | units.kg)
  
        grid[1].mass=numpy.ones((4,5))*5 | units.kg

        self.assertEqual(original.storage[1],5)
        self.assertEqual(grid[1].mass, original.storage[1] | units.kg)

class TestGridWithBinding5(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (3,4,5,6)
        
        def __init__(self):
            self.storage = numpy.arange(
                numpy.product(self.shape)).reshape(self.shape)
            
        def get_range(self):
            return (0, self.shape[0]-1, 0, self.shape[1]-1)
        def get_range2(self):
            return (0, self.shape[2]-1, 0, self.shape[3]-1)
                
        def get_a(self, i_s, j_s, k_s,l_s):
            return [numpy.asarray([self.storage[i,j,k,l] for i,j,k,l in zip(i_s, j_s, k_s,l_s)]),]
            
        def set_a(self, i_s, j_s, k_s,l_s, values):
            for i,j,k,l,val in zip(i_s, j_s,k_s,l_s,values):
                self.storage[i,j,k,l] = val
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,handler.INDEX,handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX,handler.INDEX,handler.INDEX, units.kg,), ())
              
        self.assertEqual(instance.get_a([1],[2],[3],[4]), original.storage[1,2,3,4] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_gridded_setter('grid', 'set_a','get_range2', names = ('mass',))
        handler.add_gridded_getter('grid', 'get_a','get_range2', names = ('mass',))
        
        grid = instance.grid
  
        self.assertEqual(grid[1,1].mass, original.storage[1,1] | units.kg)
        self.assertEqual(grid[0:2,1:3].mass, original.storage[0:2,1:3] | units.kg)
        self.assertEqual(len(grid[1].mass), len(original.storage[1]))
        self.assertEqual(grid[2:4].mass, original.storage[2:4] | units.kg)
        self.assertEqual(grid[-3:].mass, original.storage[-3:] | units.kg)
  
        grid[1].mass=numpy.ones((4,5,6))*5 | units.kg

        self.assertEqual(original.storage[1],5)
        self.assertEqual(grid[1].mass, original.storage[1] | units.kg)

# attempt to test string returns, needs improvement
class TestGridWithBinding6(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (3,4,5,6)
        
        def __init__(self):
            self.storage = numpy.arange(
                numpy.product(self.shape)).reshape(self.shape).astype(numpy.dtype(str))
            
        def get_range(self):
            return (0, self.shape[0]-1, 0, self.shape[1]-1)
        def get_range2(self):
            return (0, self.shape[2]-1, 0, self.shape[3]-1)
                
        def get_a(self, i_s, j_s, k_s,l_s):
            return [numpy.asarray([self.storage[i,j,k,l] for i,j,k,l in zip(i_s, j_s, k_s,l_s)]),]
            
        def set_a(self, i_s, j_s, k_s,l_s, values):
            for i,j,k,l,val in zip(i_s, j_s,k_s,l_s,values):
                self.storage[i,j,k,l] = val
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,handler.INDEX,handler.INDEX,), (handler.NO_UNIT,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX,handler.INDEX,handler.INDEX, handler.NO_UNIT,), ())
              
        self.assertEqual(instance.get_a([1],[2],[3],[4]), original.storage[1,2,3,4] )
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_gridded_setter('grid', 'set_a','get_range2', names = ('mass',))
        handler.add_gridded_getter('grid', 'get_a','get_range2', names = ('mass',))
        
        grid = instance.grid
  
        self.assertEqual(grid[1,1].mass, original.storage[1,1] )
        self.assertEqual(grid[0:2,1:3].mass, original.storage[0:2,1:3])
        self.assertEqual(len(grid[1].mass), len(original.storage[1]))
        self.assertEqual(grid[2:4].mass, original.storage[2:4] )
        self.assertEqual(grid[-3:].mass, original.storage[-3:] )
  
        grid[1].mass=numpy.ones((4,5,6))*5 

#        self.assertEquals(original.storage[1],'5')
        self.assertEqual(grid[1].mass, original.storage[1])

class TestGridWithBinding7(amusetest.TestCase):
    class TestInterface(object):
        
        shape = ()
        
        def __init__(self):
            self.storage = 123.
            
        def get_range(self):
            return ()
                
        def get_a(self):
            return self.storage
            
        def set_a(self, value):
            self.storage = value
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(), (units.kg,))
        handler.add_method('set_a',(units.kg,), ())
                      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass',))
        handler.add_getter('grid', 'get_a', names = ('mass',))
        
        grid = instance.grid
        self.assertEqual(grid.mass, 123 | units.kg)

    def test2(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(), (units.kg,))
        handler.add_method('set_a',(units.kg,), ())
                      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass',))
        handler.add_getter('grid', 'get_a', names = ('mass',))
        
        grid = instance.grid

        grid.mass=321| units.kg
        self.assertEqual(original.storage,321)
        self.assertEqual(grid.mass, 321 | units.kg)

class TestGridWithBinding8(amusetest.TestCase):
    class TestInterface(object):
        
        shape = ()
        
        def __init__(self):
            self.storage1 = 12.
            self.storage2 = 123.
            
        def get_range(self):
            return ()
                
        def get_a(self):
            return self.storage1, self.storage2
            
        def set_a(self, value1, value2):
            self.storage1 = value1
            self.storage2 = value2
        
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(), (units.kg,units.m))
        handler.add_method('set_a',(units.kg,units.m), ())
                      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass','l'))
        handler.add_getter('grid', 'get_a', names = ('mass','l'))
        
        grid = instance.grid
        self.assertEqual(grid.mass, 12 | units.kg)
        self.assertEqual(grid.l, 123 | units.m)

    def test2(self):
        original = self.TestInterface()
        
        instance = interface.InCodeComponentImplementation(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(), (units.kg,units.m))
        handler.add_method('set_a',(units.kg,units.m), ())
                      
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass','l'))
        handler.add_getter('grid', 'get_a', names = ('mass','l'))
        
        grid = instance.grid

        grid1=grid.copy()

        grid1.mass=321| units.kg
        grid1.l=32| units.m
        
        grid1.new_channel_to(grid).copy_all_attributes()
        
        self.assertEqual(original.storage1,321)
        self.assertEqual(original.storage2,32)
        self.assertEqual(grid.mass, 321 | units.kg)
        self.assertEqual(grid.l, 32 | units.m)
        
        

class CodeInterfaceAndLegacyFunctionsTest(amusetest.TestCase):
    
        
    def test1(self):
        class TestClass(object):
           
            @legacy_function
            def echo_inputs():
                function = LegacyFunctionSpecification()
                function.addParameter('input1', dtype='d', direction=function.OUT)
                function.addParameter('input2', dtype='d', direction=function.OUT)
                function.addParameter('output1', dtype='d', direction=function.IN)
                function.addParameter('output2', dtype='d', direction=function.IN)
                function.result_type = 'i'
                return function
            
        original = TestClass()
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('echo_inputs', (units.m, units.s), (units.s, handler.ERROR_CODE))
        self.assertRaises(exceptions.AmuseException, lambda: instance.echo_inputs, 
            expected_message = "Incorrect definition of method 'echo_inputs' of class 'InCodeComponentImplementation', "
            "the number of outputs do not match, expected 3, actual 2.")
            
        instance = interface.InCodeComponentImplementation(original)
        handler = instance.get_handler('METHOD')
        handler.add_method('echo_inputs', (units.m, units.s), (units.s, units.m, units.s, handler.ERROR_CODE))
        self.assertRaises(exceptions.AmuseException, lambda: instance.echo_inputs, 
            expected_message = "Incorrect definition of method 'echo_inputs' of class 'InCodeComponentImplementation', "
            "the number of outputs do not match, expected 3, actual 4.")


    def test2(self):
        class TestClass(object):
           
            @legacy_function
            def echo_inputs():
                function = LegacyFunctionSpecification()
                function.addParameter('input1', dtype='d', direction=function.OUT)
                function.addParameter('input2', dtype='d', direction=function.OUT)
                function.addParameter('output1', dtype='d', direction=function.IN)
                function.addParameter('output2', dtype='d', direction=function.IN)
                function.result_type = 'i'
                return function
            
        original = TestClass()
        
        instance = interface.InCodeComponentImplementation(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('echo_inputs', (units.m), (units.s, units.m, handler.ERROR_CODE), )
        self.assertRaises(exceptions.AmuseException, lambda: instance.echo_inputs, 
            expected_message = "Incorrect definition of method 'echo_inputs' of class 'InCodeComponentImplementation', "
            "the number of inputs do not match, expected 2, actual 1.")
    
    
    def test3(self):
        class TestClass(object):
            async_request=None
           
            @legacy_function
            def echo_inputs():
                function = LegacyFunctionSpecification()
                function.addParameter('input1', dtype='d', direction=function.OUT)
                function.addParameter('input2', dtype='d', direction=function.OUT)
                function.addParameter('output1', dtype='d', direction=function.IN)
                function.addParameter('output2', dtype='d', direction=function.IN)
                function.result_type = 'i'
                return function
            
        original = TestClass()
        self.assertRaises(exceptions.CodeException, original.echo_inputs, 1, 2,
            expected_message = "Exception when calling function 'echo_inputs', of code 'TestClass', exception was ''TestClass' object has no attribute 'channel''")
