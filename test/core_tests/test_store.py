from amuse.support.io import store
from amuse.support.units import units
from amuse.support.data.core import Stars

import unittest
import os

class TestStoreHDF(unittest.TestCase):
    
    def test1(self):
        if os.path.exists('test.hdf5'):
            os.remove('test.hdf5')
            
        instance = store.StoreHDF("test.hdf5")
        
        number_of_particles = 10
        p = Stars(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s
        
        instance.store(p)
        
        loaded_particles = instance.load()
        
        loaded_mass_in_kg = loaded_particles.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEquals(expected, actual)
        
        
        
        
