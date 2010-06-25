from amuse.support.io import store
from amuse.support import io
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.core import Particles
from amuse.test import amusetest

import os

class TestStoreHDF(amusetest.TestCase):
    
    def test1(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = store.StoreHDF(output_file)
        
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s
        
        instance.store(p)
        
        loaded_particles = instance.load()
        
        loaded_mass_in_kg = loaded_particles.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEquals(expected, actual)
        
    def test2(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
            
        instance = store.StoreHDF(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.savepoint(1 | units.Myr)
        instance.store(p.previous_state())
        
        p.mass = [x * 4.0 for x in range(number_of_particles)] | units.kg
        p.savepoint(2 | units.Myr)
        instance.store(p.previous_state())
        instance.close()
        
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load().previous_state()
        masses = loaded_particles[1].get_timeline_of_attribute("mass")
        self.assertEquals(len(masses), 2)
        
    
    def test3(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
            
        instance = store.StoreHDF(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | nbody_system.mass
        
        instance.store(p.savepoint(1 | nbody_system.time))
        instance.close()
        
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(p.mass[1], 2.0 | nbody_system.mass)
        
    
    
    def test4(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        
        particles = Particles(10)
        particles.mass = 1.0 | units.kg
        particles.x = 2.0 | units.m
        x = particles.savepoint(2.0 | units.s)
        
        io.write_set_to_file(x,"test.hdf5", format='amuse')
        
        
        particles_from_file = io.read_set_from_file("test.hdf5", format='amuse')
        
        particles_in_memory = particles_from_file.copy()
        
        self.assertAlmostRelativeEquals(particles_in_memory.mass[1], 1.0 | units.kg)
        
        particles_from_file.savepoint(4.0 | units.s)
        
        self.assertAlmostRelativeEquals(particles_from_file.mass[2], 1.0 | units.kg)
        
        

        
        
