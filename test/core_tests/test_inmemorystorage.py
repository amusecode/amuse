import numpy
import time

from amuse.test import amusetest
from amuse.units import units
from amuse.datamodel.memory_storage import InMemoryGridAttributeStorage
from amuse.datamodel.memory_storage import get_in_memory_attribute_storage_factory
from amuse.datamodel.memory_storage import InMemoryVectorQuantityAttribute

class TestInMemoryAttributeStorage(amusetest.TestCase):
    
    def test1(self):
        particles = [0,1,2]
        attributes = "a", "b"    
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        
        self.assertEquals(2.0 | units.m, instance.get_value_of(particles[1], "a"))
        

    def test2(self):
        particles = [0,1,2]
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        
        indices = instance.get_indices_of([
            particles[2], 
            particles[1], 
            particles[0]
        ])
        
        for index, wanted in zip(indices, [2,1,0]):
            self.assertEquals(index,wanted)
            
        all_values = instance.get_values_in_store(
            [particles[2], particles[0]],
            ["b","a"]
        )
        self.assertEquals(all_values[0][0],6.0 | units.g)
        self.assertEquals(all_values[0][1],4.0 | units.g)
        self.assertEquals(all_values[1][0],0.003 | units.km)
        self.assertEquals(all_values[1][1],0.001 | units.km)
        
    
    def test3(self):
        particles = [0,1,2]
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        self.assertEquals(values[0][0], 1.0 | units.m)
        
        
        instance.set_values_in_store(
            [particles[0], particles[2]],
            ["b","a"],
            [
                units.kg.new_quantity(numpy.array([9.0,11.0])),  
                units.km.new_quantity(numpy.array([1.0,2.0]))
            ]
        )
        values = instance.get_values_in_store(None, ["a", "b"])
        self.assertEquals(values[0][0], 1000.0 | units.m)
        self.assertEquals(values[0][2], 2000.0 | units.m)
        self.assertEquals(values[1][0], 9.0 | units.kg)
        self.assertEquals(values[1][2], 11.0 | units.kg)
        
    
    def test4(self):
        particles = [0,1,2, 3]
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0,4.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0,7.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        
        instance.remove_particles_from_store([particles[0], particles[2]])
        
      
        all_values = instance.get_values_in_store(
            [particles[1], particles[3]],
            ["a","b"]
        )
        self.assertEquals(all_values[0][0], 0.002 | units.km)
        self.assertEquals(all_values[0][1], 0.004 | units.km)
        self.assertEquals(all_values[1][0], 5.0 | units.g)
        self.assertEquals(all_values[1][1], 7.0 | units.g)
        
    
    def test5(self):
        particles = [0,1,2,3]
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0,4.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0,7.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        
        self.assertEquals(len(instance), 4)
        
        particles = [4,5,6,7]
        instance.add_particles_to_store(particles, attributes, values)
        
        self.assertEquals(len(instance), 8)
        
        all_values = instance.get_values_in_store([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0 | units.m)
        self.assertEquals(all_values[0][1], 2.0 | units.m)
        
        instance.remove_particles_from_store((0, 4))
        
        self.assertEquals(len(instance), 6)
        
        all_values = instance.get_values_in_store([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0 | units.m)
        self.assertEquals(all_values[0][1], 2.0 | units.m)
        
    def test6(self):
        particles = [0,1,2,3]
        attributes = "a", "b"
        values = [
            numpy.array([1.0,2.0,3.0,4.0]), 
            numpy.array([4.0,5.0,6.0,7.0])
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(particles, attributes, values)
        
        self.assertEquals(len(instance), 4)
        
        particles = [4,5,6,7]
        instance.add_particles_to_store(particles, attributes, values)
        
        self.assertEquals(len(instance), 8)
        
        all_values = instance.get_values_in_store([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0)
        self.assertEquals(all_values[0][1], 2.0)
        
        instance.remove_particles_from_store((0, 4))
        
        self.assertEquals(len(instance), 6)
        
        all_values = instance.get_values_in_store([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0)
        self.assertEquals(all_values[0][1], 2.0)
        
        instance.set_values_in_store([1,5], ["a"], [[4.0, 5.0]])
        all_values = instance.get_values_in_store([1,5], ["a"])[0]
        self.assertEquals(all_values[0], 4.0)
        self.assertEquals(all_values[1], 5.0)
        
        instance2 = instance.copy()
        instance2.set_values_in_store([1,5], ["a"], [[3.0, 1.0]])
        all_values = instance2.get_values_in_store([1,5], ["a"])[0]
        self.assertEquals(all_values[0], 3.0)
        self.assertEquals(all_values[1], 1.0)
        all_values = instance.get_values_in_store([1,5], ["a"])[0]
        self.assertEquals(all_values[0], 4.0)
        self.assertEquals(all_values[1], 5.0)
        

    def test7(self):
        n = 1000000
        keys = numpy.arange(0, n)
        attributes = "a", "b"
        values = [
            numpy.ones(n),
            numpy.ones(n) * 2,
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance.add_particles_to_store(keys, attributes, values)
        
        self.assertEquals(len(instance), n)
        
        t0 = time.time()
        all_values = instance.get_values_in_store(keys, ["a"], by_key = True)
        t1 = time.time()
        dt_by_key = t1 - t0
        
        t0 = time.time()
        all_values = instance.get_values_in_store(keys, ["a"], by_key = False)
        t1 = time.time()
        dt_by_index = t1 - t0
        self.assertTrue(dt_by_index < dt_by_key)
        
        
        t0 = time.time()
        all_values = instance.set_values_in_store(keys, ["a"], [values[1]],by_key = True)
        t1 = time.time()
        dt_by_key = t1 - t0
        
        t0 = time.time()
        all_values = instance.set_values_in_store(keys, ["a"], [values[1]], by_key = False)
        t1 = time.time()
        dt_by_index = t1 - t0
        print dt_by_index, dt_by_key
        self.assertTrue(dt_by_index < dt_by_key)
        
class TestInMemoryGridAttributeStorage(amusetest.TestCase):
    
    def test1(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x.set_values_in_store(
            (i,j,k), 
            ['a','b'], 
            [2.0 | units.kg, 1.0 | units.m]
        )
        
        b, a = x.get_values_in_store(None, ['b','a'])
        print a
        self.assertEquals(b[0][1][0], 1.0 | units.m)
        self.assertEquals(b[0][0][0], 0.0 | units.m)
        self.assertEquals(a[0][1][0], 2.0 | units.kg)
        self.assertEquals(a[1][3][2], 2.0 | units.kg)
        self.assertEquals(a[1][2][2], 0.0 | units.kg)
        
        (b,) = x.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0 | units.kg)
        self.assertEquals(b[0][0][2], 0.0 | units.kg)
        self.assertEquals(b[1][2][2], 2.0 | units.kg)
        self.assertEquals(b[2][0][0], 2.0 | units.kg)
        self.assertEquals(b[3][2][2], 2.0 | units.kg)
        
        self.assertEquals(b.sum(), 8.0 | units.kg)
        
        self.assertEquals(sorted(x.get_defined_attribute_names()), ["a", "b"])
    
    def test2(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x.set_values_in_store(
            (i,j,k), 
            ['a','b'], 
            [2.0 , 1.0]
        )
        
        b, a = x.get_values_in_store(None, ['b','a'])
        print a
        self.assertEquals(b[0][1][0], 1.0 )
        self.assertEquals(b[0][0][0], 0.0 )
        self.assertEquals(a[0][1][0], 2.0 )
        self.assertEquals(a[1][3][2], 2.0  )
        self.assertEquals(a[1][2][2], 0.0  )
        
        (b,) = x.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0  )
        self.assertEquals(b[0][0][2], 0.0  )
        self.assertEquals(b[1][2][2], 2.0  )
        self.assertEquals(b[2][0][0], 2.0  )
        self.assertEquals(b[3][2][2], 2.0  )
        
        self.assertEquals(b.sum(), 8.0  )
        
        self.assertEquals(sorted(x.get_defined_attribute_names()), ["a", "b"])
        
        y = x.copy()
        
        (b,) = y.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0  )
        self.assertEquals(b[0][0][2], 0.0  )
        self.assertEquals(b[1][2][2], 2.0  )
        
        
        
class TestInMemoryVectorQuantityAttribute(amusetest.TestCase):
    
    def test1(self):
        quantity  = units.m.new_quantity(numpy.array([1.0,2.0,3.0]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        
        self.assertEquals(attribute.get_length(), 3)
        self.assertEquals(attribute.get_shape(), (3,) )
        self.assertEquals(attribute.get_values([1,2]), [2.0,3.0] | units.m)
    
        attribute.increase_to_length(5)
        self.assertEquals(attribute.get_values(None), [1.0,2.0,3.0,0.0,0.0] | units.m)
        
    def test2(self):
        quantity  = units.m.new_quantity(numpy.array([1.0,2.0,3.0]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        attribute.set_values([1,2], [4.0,5.0] | units.m)
    
        attribute.increase_to_length(5)
        self.assertEquals(attribute.get_values(None), [1.0,4.0,5.0,0.0,0.0] | units.m)

    def test3(self):
        quantity  = units.m.new_quantity(numpy.array([[1.0,2.0,3.0], [4.0,5.0,6.0]]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        self.assertEquals( attribute.get_values([1]),  [4.0,5.0,6.0] | units.m)
        self.assertEquals( attribute.get_shape(), (2,3))
        attribute.increase_to_length(4)
        self.assertEquals( attribute.get_shape(), (4,3))
        self.assertEquals(attribute.get_values(None), [[1.0,2.0,3.0], [4.0,5.0,6.0], [0.0,0.0,0.0], [0.0,0.0,0.0]] | units.m)

