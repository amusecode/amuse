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
        
        print p.attributelist
        
        instance.store(p)
        #self.assertTrue(False)
        
        
        
