from amuse.test import amusetest
import numpy

from amuse.support.data.memory_storage import InMemoryGridAttributeStorage
from amuse.support.data.memory_storage import get_in_memory_attribute_storage_factory

from amuse.support.units  import units

class TestInMemoryAttributeStorage(amusetest.TestCase):
    
    def test1(self):
        particles = [0,1,2]
        attributes = "a", "b"    
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance._add_particles(particles, attributes, values)
        
        self.assertEquals(2.0 | units.m, instance.get_value_of(particles[1], "a"))
        

    def test2(self):
        particles = [0,1,2]
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0,2.0,3.0])), 
            units.g.new_quantity(numpy.array([4.0,5.0,6.0]))
        ]
        
        instance = get_in_memory_attribute_storage_factory()()
        instance._add_particles(particles, attributes, values)
        
        indices = instance.get_indices_of([
            particles[2], 
            particles[1], 
            particles[0]
        ])
        
        for index, wanted in zip(indices, [2,1,0]):
            self.assertEquals(index,wanted)
            
        all_values = instance._get_values(
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
        instance._add_particles(particles, attributes, values)
        self.assertEquals(values[0][0], 1.0 | units.m)
        
        
        instance._set_values(
            [particles[0], particles[2]],
            ["b","a"],
            [
                units.kg.new_quantity(numpy.array([9.0,11.0])),  
                units.km.new_quantity(numpy.array([1.0,2.0]))
            ]
        )
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
        instance._add_particles(particles, attributes, values)
        
        instance._remove_particles([particles[0], particles[2]])
        
      
        all_values = instance._get_values(
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
        instance._add_particles(particles, attributes, values)
        
        self.assertEquals(len(instance), 4)
        
        particles = [4,5,6,7]
        instance._add_particles(particles, attributes, values)
        
        self.assertEquals(len(instance), 8)
        
        all_values = instance._get_values([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0 | units.m)
        self.assertEquals(all_values[0][1], 2.0 | units.m)
        
        instance._remove_particles((0, 4))
        
        self.assertEquals(len(instance), 6)
        
        all_values = instance._get_values([1,5], ["a"])
        
        self.assertEquals(all_values[0][0], 2.0 | units.m)
        self.assertEquals(all_values[0][1], 2.0 | units.m)
        

class TestInMemoryGridAttributeStorage(amusetest.TestCase):
    
    def test1(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x._set_values(
            (i,j,k), 
            ['a','b'], 
            [2.0 | units.kg, 1.0 | units.m]
        )
        
        b, a = x._get_values(None, ['b','a'])
        print a
        self.assertEquals(b[0][1][0], 1.0 | units.m)
        self.assertEquals(b[0][0][0], 0.0 | units.m)
        self.assertEquals(a[0][1][0], 2.0 | units.kg)
        self.assertEquals(a[1][3][2], 2.0 | units.kg)
        self.assertEquals(a[1][2][2], 0.0 | units.kg)
        
        (b,) = x._get_values((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0 | units.kg)
        self.assertEquals(b[0][0][2], 0.0 | units.kg)
        self.assertEquals(b[1][2][2], 2.0 | units.kg)
        self.assertEquals(b[2][0][0], 2.0 | units.kg)
        self.assertEquals(b[3][2][2], 2.0 | units.kg)
        
        self.assertEquals(b.sum(), 8.0 | units.kg)
        
        self.assertEquals(sorted(x._get_attribute_names()), ["a", "b"])
        
        
        
