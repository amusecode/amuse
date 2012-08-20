from amuse.test import amusetest

import os
import numpy
import time

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
        output_file = os.path.join(test_results_path, "test12.hdf5")
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
    
    def test13(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test13.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[1,2,3],[4,5,6]] | units.km
       

        io.write_set_to_file(stars, output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5")
        self.assertEquals(loaded[0].md, [1,2,3] | units.km)
        self.assertEquals(loaded[1].md, [4,5,6] | units.km)
        
        self.assertEquals(loaded.md[0], [1,2,3] | units.km)
        self.assertEquals(loaded.md[1], [4,5,6] | units.km)
         
        self.assertEquals(loaded.previous_state()[0].md, [1,2,3] | units.km)
        #self.assertEquals(loaded.previous_state()[0].md,  [7,8,9] | units.km)
    
      
    def test14(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test14.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5")
        self.assertEquals(loaded[0].md, [[1,3],[2,4],[3,5]])
        self.assertEquals(loaded[1].md, [[4,6],[5,7],[6,8]])
        
        self.assertEquals(loaded.md[0], [[1,3],[2,4],[3,5]])
        self.assertEquals(loaded.md[1], [[4,6],[5,7],[6,8]])
    
    def test15(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test15.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5")
        processor = store.StoreHDF(output_file, True, open_for_writing = True)
        loaded = processor.load()
        loaded.previous_state()[0].md = [[3,1],[3,4],[5,2]] 
        self.assertEquals(loaded.previous_state()[0].md, [[3,1],[3,4],[5,2]] )
        processor.close()
        
        loaded = io.read_set_from_file(output_file, "hdf5")
        self.assertEquals(loaded[0].md, [[3,1],[3,4],[5,2]])
        self.assertEquals(loaded[1].md, [[4,6],[5,7],[6,8]])
        
    def test16(self):
        import h5py
        print h5py.version.version
        print h5py
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test16.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        stars[0].nn = stars[1]
        stars[1].nn = stars[0]
       

        self.assertEqual(stars[0].nn,stars[1])
        self.assertEqual(stars[1].nn,stars[0])
        
        io.write_set_to_file(stars, output_file, "hdf5")
        processor = store.StoreHDF(output_file, True, open_for_writing = True)
        loaded = processor.load()
        
        self.assertEqual(loaded[0].nn,loaded[1])
        self.assertEqual(loaded[1].nn,loaded[0])
        self.assertEqual(loaded[0].nn.key,stars[1].key)
        self.assertEqual(loaded[0].nn.x,stars[1].x)
        
    
    def test17(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test17.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(4)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        stars[2].x = 3.0  | units.km
        stars[3].x = 4.0  | units.km
        
        binaries = Particles(2)
        binaries[0].y = 1.0 | units.km
        binaries[1].y = 2.0 | units.km
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]
        binaries[1].child1 = stars[2]
        binaries[1].child2 = stars[3]

        self.assertEqual(binaries[0].child1,stars[0])
        self.assertEqual(binaries[1].child1,stars[2])
        
        io.write_set_to_file(binaries, output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5")
        
        
        self.assertEqual(loaded[0].child1.key,stars[0].key)
        self.assertEqual(loaded[1].child1.key,stars[2].key)
        
    def test18(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test18.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        binaries = Particles(1)
        binaries[0].y = 1.0 | units.km
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]

        stars[0].parent = binaries[0]
        stars[1].parent = binaries[0]
        
        self.assertEqual(binaries[0].child1,stars[0])
        self.assertEqual(binaries[0].child2,stars[1])
        self.assertEqual(binaries[0].child1.parent,binaries[0])
        self.assertEqual(binaries[0].child2.parent,binaries[0])
        
        io.write_set_to_file(binaries, output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5")
        
        self.assertEqual(loaded[0].child1.key,stars[0].key)
        self.assertEqual(loaded[0].child2.key,stars[1].key)
        self.assertEqual(loaded[0].child1,stars[0])
        self.assertEqual(loaded[0].child2,stars[1])
        self.assertEqual(loaded[0].child1.parent,loaded[0])
        self.assertEqual(loaded[0].child2.parent,loaded[0])
    
    def test19(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test18.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        binaries = Particles(1)
        binaries[0].y = 1.0 | units.km
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]

        stars[0].parent = binaries[0]
        stars[1].parent = binaries[0]
        
        self.assertEqual(binaries[0].child1,stars[0])
        self.assertEqual(binaries[0].child2,stars[1])
        self.assertEqual(binaries[0].child1.parent,binaries[0])
        self.assertEqual(binaries[0].child2.parent,binaries[0])
        
        io.write_set_to_file([binaries,stars], output_file, "hdf5", names = ['binaries', 'children'])
        
        loader_binaries, loaded_stars = io.read_set_from_file(output_file, "hdf5", names = ['binaries', 'children'])
        
        self.assertEqual(loader_binaries[0].child1.key,stars[0].key)
        self.assertEqual(loader_binaries[0].child2.key,stars[1].key)
        self.assertEqual(loader_binaries[0].child1,loaded_stars[0])
        self.assertEqual(loader_binaries[0].child2,loaded_stars[1])
        self.assertEqual(loader_binaries[0].child1.parent,loader_binaries[0])
        self.assertEqual(loader_binaries[0].child2.parent,loader_binaries[0])

    def test20(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test19.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = (numpy.asarray([x * 2.0 for x in range(p.size)])).reshape(shape)

        io.write_set_to_file(p, output_file, "hdf5")
        
        loaded = io.read_set_from_file(output_file, "hdf5").previous_state()
        self.assertAlmostRelativeEquals(p.mass[0][1][2], 24)
        self.assertAlmostRelativeEquals(p[0][1][2].mass, 24)
        
    def test21(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test15.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5")
       
        
        loaded = io.read_set_from_file(output_file, "hdf5", close_file = True)
        self.assertEquals(loaded[0].md, [[1,3],[2,4],[3,5]])
        self.assertEquals(loaded[1].md, [[4,6],[5,7],[6,8]])
        
        previous = loaded.previous_state()
        self.assertEquals(previous, None)
