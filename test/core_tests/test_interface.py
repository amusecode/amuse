from amuse.support import interface
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.binding import *
from amuse.support.data.parameters import *
from amuse.support.data import core
from amuse.support.core import OrderedDictionary

from amuse.test import amusetest
import numpy

class CodeInterfaceWithConvertedUnitsTests(amusetest.TestCase):
    class TestClass(object):
        
        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + 10.0
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (nbody_system.length,), nbody_system.length)
        
        print instance.mass
        self.assertAlmostEquals(instance.mass.value_in(units.kg), 100.0, 10)
        
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (nbody_system.length,), nbody_system.length)
                
        self.assertAlmostEquals(instance.add_to_length(5|units.m).value_in(units.m), 55.0, 10)
        
        
        
    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        instance = interface.CodeInterface(original)

        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        handler = instance.get_handler('PARAMETER')
        handler.add_attribute_parameter(
            "eps",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.3 | nbody_system.length * nbody_system.length
        )
          
        instance.parameters.epsilon_squared = 100.0 | units.m ** 2
        self.assertAlmostEquals(instance.parameters.epsilon_squared.value_in(units.m**2),  100.0, 6)
        self.assertAlmostEquals(original.eps,  4.0, 6)

    def test4(self):
        original = self.TestClass()
        instance = interface.CodeInterface(original)

        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add_10')
           
        self.assertFalse(instance.add_10.is_async_supported)

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
        
        instance = interface.CodeInterface(original)

        handler = instance.get_handler('METHOD')
        handler.add_method('add_10_to_length', (units.m,), units.m)
        
        self.assertEquals(20.0, original.add_10_to_length(10.0))
        
        self.assertEquals(20.0 | units.m, instance.add_10_to_length(10.0 | units.m))
        self.assertEquals(1010.0 | units.m, instance.add_10_to_length(1.0| units.km))
        
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_10_to_length', (units.m,), units.m, public_name = 'add_10')
           
        self.assertEquals(20.0 | units.m, instance.add_10(10.0 | units.m))
        self.assertEquals(1010.0 | units.m, instance.add_10(1.0| units.km))
        
    def test3(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_one', units.m)
        
        self.assertEquals(1.0 | units.m, instance.one)
        
    def test4(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_one', units.m, public_name = 'get_one')
        
        self.assertEquals(1.0 | units.m, instance.get_one)
        
    
    def test5(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_state', (handler.NO_UNIT,), (units.m, units.m, units.kg, handler.ERROR_CODE))
        
        result = instance.get_state(1)
        self.assertEquals(3, len(result))
        self.assertEquals(1.0 | units.m, result[0])
        self.assertEquals(3.0 | units.kg, result[2])
        
    
    def test6(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_state_error', (handler.NO_UNIT,), (units.m, units.m, units.kg, handler.ERROR_CODE))
        
        try:
            result = instance.get_state_error(1)
            self.fail("should raise an exception")
        except:
            
            pass
            
class CodeInterfaceTests(amusetest.TestCase):
    class TestClass(object):
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
        
        instance = interface.CodeInterface(original)
        
        self.assertTrue(instance.always_works())
        instance.move_to_state_1()
        self.assertEquals(1, instance.returns_1())
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        handler = instance.get_handler('STATE')
        
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        self.assertTrue(instance.always_works())
        self.assertEquals(handler._current_state.name, 'ZERO')
        instance.move_to_state_1()
        
        self.assertEquals(handler._current_state.name, 'ONE')
        self.assertEquals(instance.returns_1(), 1)
        
    
    def test3(self):
        class TestCodeInterface(interface.CodeInterface):
            
            def move_to_state_1(self):
                self.overridden().move_to_state_1()
                self.traced = True
                
        original = self.TestClass()
        
        instance = TestCodeInterface(original)
        instance.traced = False
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
            
        
        
        self.assertEquals(instance.get_name_of_current_state(), 'ZERO')
        self.assertEquals(instance.returns_1(), 1)
        self.assertEquals(instance.get_name_of_current_state(), 'ONE')
        self.assertTrue(instance.traced)
        
        
    def test4(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('STATE')
        
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_method('ONE', 'returns_1')
        handler.set_initial_state('ZERO')
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_1(), 1)        
        self.assertEquals(handler._current_state.name, 'ONE')
        
    def test5(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.set_initial_state('ZERO')
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_2(), 2)        
        self.assertEquals(handler._current_state.name, 'TWO')
        
    def test6(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_4(), 4)        
        self.assertEquals(handler._current_state.name, 'FOUR')
        
        
        handler.set_initial_state('ZERO')
        self.assertEquals(handler._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_3(), 3)        
        self.assertEquals(handler._current_state.name, 'THREE')
        self.assertEquals(instance.returns_4(), 4)        
        self.assertEquals(handler._current_state.name, 'FOUR')
        
    
    def test7(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        try:
            self.assertEquals(instance.returns_4(), 4)
            self.fail("Automatic state transitions is OFF, this method should fail")  
        except Exception, ex:
            print ex
            
            self.assertEquals(len(ex.transitions), 3)
            
            for x in ex.transitions:
                x.do()
            
            self.assertEquals(instance.returns_4(), 4)
            
    
    def test8(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.set_initial_state('ZERO')
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        self.assertEquals(instance.returns_2(), 2)
        try:
            self.assertEquals(instance.returns_3(), 3)
            self.fail("No transition to state 3 possible, function should error")  
        except Exception, ex:
            print ex
        
    def test8(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('STATE')
        handler.add_transition('ZERO', 'ONE', 'move_to_state_1')
        handler.add_transition('ONE', 'TWO', 'move_to_state_2')
        handler.add_transition('THREE', 'ONE', 'move_to_state_1')
        handler.add_transition('TWO', 'THREE', 'move_to_state_1')
        
        
        handler.add_method('ONE', 'returns_1')
        handler.add_method('TWO', 'returns_2')
        handler.add_method('THREE', 'returns_3')
        handler.set_initial_state('ZERO')
        
        
        self.assertEquals(handler._current_state.name, 'ZERO')
        instance.move_to_state_1()
        self.assertEquals(handler._current_state.name, 'ONE')
        instance.move_to_state_2()
        self.assertEquals(handler._current_state.name, 'TWO')
        instance.move_to_state_1()
        self.assertEquals(handler._current_state.name, 'THREE')
        instance.move_to_state_1()
        self.assertEquals(handler._current_state.name, 'ONE')
        
        self.assertEquals(instance.returns_3(), 1)
        


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
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m)
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add_to_length')
        
        self.assertEquals(40.0 | units.m, instance.add_to_length(20.0 | units.m))
        
    def test2(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add')
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add_to_length')
        
        self.assertEquals(40.0 | units.m, instance.add(20.0 | units.m))
        
    def test3(self):
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        
        handler = instance.get_handler('METHOD')
        handler.add_method('add_to_length', (units.m,), units.m, public_name = 'add')
        
        
        handler = instance.get_handler('STATE')
        handler.set_initial_state('ZERO')
        handler.add_transition('ZERO', 'ONE', 'move_to_20')
        handler.add_method('ONE', 'add')
        
        self.assertEquals(40.0 | units.m, instance.add(20.0 | units.m))
        

class CodeInterfaceWithErrorHandlingTests(amusetest.TestCase):
    class TestClass(object):
        errorcode = 0
        
        def get_mass(self):
            print "EC:", self.errorcode
            return 10.0, self.errorcode
            
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass', (), (units.m, handler.ERROR_CODE,))
        handler = instance.get_handler('ERRORCODE')
        handler.add_errorcode(-2, "no such method")
        handler.add_errorcode(-3, "not available")
        
        self.assertEquals(instance.get_mass(), 10.0 | units.m)
        try:
            original.errorcode = -2
            instance.get_mass()
        except Exception, ex:
            self.assertEquals("Error when calling 'get_mass' of a 'CodeInterface', errorcode is -2, error is 'no such method'",str(ex))
        else:
            self.fail("should raise an exception")
            
        try:
            original.errorcode = -1
            instance.get_mass()
        except Exception, ex:
            self.assertEquals("Error when calling 'get_mass' of a 'CodeInterface', errorcode is -1", str(ex))
        else:
            self.fail("should raise an exception")
            

class CodeInterfaceWithParticlesTests(amusetest.TestCase):
    class TestClass(object):

        def get_mass(self):
            return 10.0, 0
            
        def add_to_length(self, length):
            return length + 10.0
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        original = self.TestClass()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('PROPERTY')
        handler.add_property('get_mass', nbody_system.mass)
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)
        
        

class TestParticlesWithBinding(amusetest.TestCase):
    class TestInterface(object):

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
            return (map(lambda x : x+1, id), map(lambda x : 0, id))
            
        def add_1_to_mass(self, id):
            if isinstance(id, int):
                self.masses[id] += 1.0
                return [0]
            for i in id:
                self.masses[i] += 1.0
            return map(lambda x : 0, id)
            
        
            
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        local_particles = core.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(original.masses), 4)
        self.assertEquals(original.masses[0], 3.0)
        self.assertEquals(original.masses[3], 6.0)
        
        self.assertEquals(len(instance.particles), 4)
        self.assertEquals(instance.particles[0].mass, 3.0 | units.kg)
        
    def test2(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        local_particles = core.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        colliding_particles = instance.particles.select_colliding()
        
        self.assertEquals(len(colliding_particles), 2)
        self.assertEquals(colliding_particles.mass[0], 4.0 | units.kg)
        self.assertEquals(colliding_particles.mass[1], 5.0 | units.kg)
        
        attribute_names = dir(instance)
        
        self.assertTrue('particles' in attribute_names)
        self.assertTrue('get_colliding_indices' in attribute_names)
        
    def test3(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_mass',(handler.NO_UNIT,), (units.kg, handler.ERROR_CODE))
        handler.add_method('set_mass',(handler.NO_UNIT, units.kg,), (handler.ERROR_CODE,))
        handler.add_method('get_next',(handler.INDEX,), (handler.INDEX, handler.ERROR_CODE))
        handler.add_method('new_particle',(units.kg,), (handler.NO_UNIT, handler.ERROR_CODE))
        handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
        handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler = instance.get_handler('PARTICLES')
        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_getter('particles', 'get_next', names = ('next',))
        
        
        local_particles = core.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        
        self.assertEquals(len(instance.particles), 4)
        self.assertEquals(instance.particles[0].next, instance.particles[1])
        self.assertEquals(instance.particles[1].next, instance.particles[2])
        self.assertEquals(instance.particles[2].next, instance.particles[3])
        self.assertEquals(instance.particles[3].next, None)
        
        
    def test4(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        local_particles = core.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(instance.particles), 4)
        self.assertEquals(instance.particles.mass, units.kg.new_quantity([3.0, 4.0, 5.0, 6.0]))
        instance.particles[0].add_one()
        self.assertEquals(instance.particles.mass, units.kg.new_quantity([4.0, 4.0, 5.0, 6.0]))
        instance.particles.add_one()
        self.assertEquals(instance.particles.mass, units.kg.new_quantity([5.0, 5.0, 6.0, 7.0]))
        
        
    def test5(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
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
        
        
        convert_nbody = nbody_system.nbody_to_si(10.0 | units.kg , 5.0 | units.m )
        
        handler = instance.get_handler('UNIT')
        handler.set_nbody_converter(convert_nbody)

        local_particles = core.Particles(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = instance.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(original.masses), 4)
        self.assertAlmostEquals(original.masses[0], 0.3, 5)
        self.assertAlmostEquals(original.masses[3], 0.6 , 5)
        
        self.assertEquals(len(instance.particles), 4)
        self.assertEquals(instance.particles[0].mass, 3.0 | units.kg)
        
        
        
class TestGridWithBinding(amusetest.TestCase):
    class TestInterface(object):
        
        shape = (11,5,5)
        
        def __init__(self):
            self.storage = numpy.arange(
                self.shape[0]*self.shape[1]*self.shape[2]).reshape(self.shape)
            
        def get_range(self):
            return (0,self.shape[0]-1,0,self.shape[1]-1,0,self.shape[2]-1)
            
        def get_a(self,i_s,j_s,k_s):
            #print "indices:", i_s, j_s, k_s
            #print "values:", numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)])
            return [numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)]),]
            
        def set_a(self, i_s, j_s, k_s, values):
            index = 0
            for i,j,k in zip(i_s, j_s, k_s):
                self.storage[i][j][k] = values[index].value_in(units.m)
                index += 1
                print index
        
            
        
            
    def test1(self):
        original = self.TestInterface()
        
        instance = interface.CodeInterface(original)
        
        handler = instance.get_handler('METHOD')
        handler.add_method('get_a',(handler.INDEX, handler.INDEX,handler.INDEX,), (units.kg,))
        handler.add_method('set_a',(handler.INDEX, handler.INDEX,handler.INDEX, units.kg,), ())
      
        print instance.get_a([1],[2],[3])
        print [38] | units.kg
        
        self.assertEquals(instance.get_a([1],[2],[3]), [38] | units.kg)
        
        handler = instance.get_handler('PARTICLES')
        handler.define_grid('grid',)
        handler.add_setter('grid', 'set_a', names = ('mass',))
        handler.add_getter('grid', 'get_a', names = ('mass',))
        
        grid = instance.grid
        
        
        self.assertEquals(grid[1][2][3].mass, 38 | units.kg)
        self.assertEquals(grid[0:2][1][2][3].mass, 38 | units.kg)
        self.assertEquals(len(grid[1][2].mass), 5)
        self.assertTrue(numpy.all(grid[1][2].mass == [35,36,37,38,39] | units.kg))
        
