import unittest
import numpy

from amuse.support.data.core import AttributeList, Particle
from amuse.support.units  import units

class TestAttributeList(unittest.TestCase):
    
    def test1(self):
        particles = [Particle(0), Particle(1), Particle(2)]
        attributes = "a", "b"
        values = [[1,2,3], [4,5,6]]
        units_of_attributes = [units.m, units.kg]
        
        instance = AttributeList(particles, attributes, values, units_of_attributes)
        
        self.assertEquals(2.0 | units.m, instance.get_value_of(particles[1], "a"))
        

    def test2(self):
        particles = [Particle(0), Particle(1), Particle(2)]
        attributes = "a", "b"
        values = [numpy.array([1.0,2.0,3.0]),  numpy.array([4.0,5.0,6.0])]
        units_of_attributes = [units.m, units.kg]
        
        instance = AttributeList(particles, attributes, values, units_of_attributes)
        
        indices = instance.get_indices_of([
            particles[2], 
            particles[1], 
            particles[0]
        ])
        
        for index, wanted in zip(indices, [2,1,0]):
            self.assertEquals(index,wanted)
            
        all_values = instance.get_values_of_particles_in_units(
            [particles[2], particles[0]],
            ["b","a"],
            [units.kg, units.km]
        )
        self.assertEquals(all_values[0][0],6.0)
        self.assertEquals(all_values[0][1],4.0)
        self.assertEquals(all_values[1][0],0.003)
        self.assertEquals(all_values[1][1],0.001)
        
    
    def test3(self):
        particles = [Particle(0), Particle(1), Particle(2)]
        attributes = "a", "b"
        values = [numpy.array([1.0,2.0,3.0]),  numpy.array([4.0,5.0,6.0])]
        units_of_attributes = [units.m, units.kg]
        
        instance = AttributeList(particles, attributes, values, units_of_attributes)
        
        
        instance.set_values_of_particles_in_units(
            [particles[0], particles[2]],
            ["b","a"],
            [numpy.array([9.0,11.0]),  numpy.array([1.0,2.0])],
            [units.kg, units.km]
        )
        self.assertEquals(values[0][0], 1000.0)
        self.assertEquals(values[0][2], 2000.0)
        self.assertEquals(values[1][0], 9.0)
        self.assertEquals(values[1][2], 11.0)
        
    
    def test4(self):
        particles = [Particle(0), Particle(1), Particle(2), Particle(3)]
        attributes = "a", "b"
        values = [numpy.array([1.0,2.0,3.0, 4.0]),  numpy.array([4.0,5.0,6.0, 7.0])]
        units_of_attributes = [units.m, units.kg]
        
        instance = AttributeList(particles, attributes, values, units_of_attributes)
        
        instance.remove_particles([particles[0], particles[2]])
        
      
        all_values = instance.get_values_of_particles_in_units(
            [particles[1], particles[3]],
            ["a","b"],
            [units.km, units.kg]
        )
        self.assertEquals(all_values[0][0], 0.002)
        self.assertEquals(all_values[0][1], 0.004)
        self.assertEquals(all_values[1][0], 5.0)
        self.assertEquals(all_values[1][1], 7.0)
        
