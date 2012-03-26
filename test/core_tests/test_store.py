from amuse.test import amusetest

import os
import numpy

from amuse import io
from amuse.io import store
from amuse.units import units
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import Grid

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
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(loaded_particles.previous_state().get_timestamp(), 2 | units.Myr)
        self.assertAlmostRelativeEquals(loaded_particles.previous_state().previous_state().get_timestamp(), 1 | units.Myr)
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

        io.write_set_to_file(x,output_file, format='amuse')

        particles_from_file = io.read_set_from_file(output_file, format='amuse')
        particles_in_memory = particles_from_file.copy()

        self.assertAlmostRelativeEquals(particles_in_memory.mass[1], 1.0 | units.kg)

        particles_from_file.savepoint(4.0 | units.s)

        self.assertAlmostRelativeEquals(particles_from_file.mass[2], 1.0 | units.kg)
    
    def test5(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "testgrid.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = store.StoreHDF(output_file)
        
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        instance.store_grid(p)

        loaded_grid = instance.load_grid()

        self.assertEquals(loaded_grid.shape, shape)
        
        loaded_mass_in_kg = loaded_grid.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEquals(expected, actual)
    
    def test6(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "testgrid.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = store.StoreHDF(output_file)
        
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s
        
        instance.store_grid(p.savepoint(1| units.Myr))

        loaded_grid = instance.load_grid().previous_state()
        
        self.assertAlmostRelativeEquals(loaded_grid.get_timestamp(), 1| units.Myr)
        
        self.assertEquals(loaded_grid.shape, shape)
        self.assertEquals(loaded_grid[0][0][0].mass, 0 | units.kg)
        self.assertAlmostRelativeEquals(loaded_grid[...,1,1].mass, p[...,1,1].mass)

    def test7(self):
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
        instance.close()
    
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load().previous_state()

        loaded_particles.mass = [x * 3.0 for x in range(number_of_particles)] | units.kg
        previous_mass_in_kg = [x * 3.0 for x in range(number_of_particles)]
        instance.close()
        
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load()
        loaded_mass_in_kg = loaded_particles.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEquals(expected, actual)
        
        instance.close()
        
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load().previous_state()
        loaded_particles[2].mass = 44 | units.kg
        instance.close()
        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load().previous_state()
        self.assertEquals( 44 | units.kg, loaded_particles[2].mass)
        instance.close()
        
    def test8(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        instance = store.StoreHDF(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = numpy.asarray([x * 2.0 for x in range(number_of_particles)])

        instance.store(p.savepoint(1 | nbody_system.time))
        instance.close()

        instance = store.StoreHDF(output_file)
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(p.mass[1], 2.0)
        
    

    def test9(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test9.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p, output_file, "hdf5",  timestamp = 2 | units.Myr, scale = 1 | units.kg)
        
        loaded_particles = io.read_set_from_file(output_file, "hdf5").previous_state()
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)


    def test10(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test10.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p, output_file, "hdf5",  timestamp = 2 | units.Myr, scale = 1 | units.kg)
        
        loaded = io.read_set_from_file(output_file, "hdf5").previous_state()
        a = loaded.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)


    def test11(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test11.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p.savepoint(timestamp = 2 | units.Myr, scale = 1 | units.kg), output_file, "hdf5")
        
        loaded_particles = io.read_set_from_file(output_file, "hdf5").previous_state()
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)
        
    
    def test12(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test10.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p.savepoint(timestamp = 2 | units.Myr, scale = 1 | units.kg), output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5").previous_state()
        a = loaded.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)

