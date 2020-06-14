from amuse.test import amusetest

import os
import numpy
import time

from amuse import io
from amuse.io import store
from amuse.io import store_v1
from amuse.io import store_v2
from amuse.units import units
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesOverlay
from amuse.datamodel import Grid
from amuse.datamodel import ParticlesSuperset




class _AbstractTestStoreHDF(object):
    
    def store_factory(self):
        raise NotImplementedError
    
    def store_version(self):
        raise NotImplementedError
        
    def get_version_in_store(self, container):
        raise NotImplementedError

    def test1(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test1."+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = self.store_factory()(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        instance.store(p)

        loaded_particles = instance.load()

        loaded_mass_in_kg = loaded_particles.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEqual(expected, actual)
        
        instance.close()

    def test2(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test2."+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        instance = self.store_factory()(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        instance.store(p.savepoint(1 | units.Myr))

        p.mass = [x * 4.0 for x in range(number_of_particles)] | units.kg
        instance.store(p.savepoint(2 | units.Myr))
        instance.close()

        instance = self.store_factory()(output_file)
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(self.get_version_in_store(loaded_particles).get_timestamp(), 2 | units.Myr)
        self.assertAlmostRelativeEquals(self.get_version_in_store(loaded_particles).previous_state().get_timestamp(), 1 | units.Myr)
        masses = loaded_particles[1].get_timeline_of_attribute("mass")
        self.assertEqual(len(masses), 2)

        instance.close()

    def test3(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test3"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        instance = self.store_factory()(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | nbody_system.mass

        instance.store(p.savepoint(1 | nbody_system.time))
        instance.close()

        instance = self.store_factory()(output_file)
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(p.mass[1], 2.0 | nbody_system.mass)

        instance.close()

    def test4(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test4"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Particles(10)
        particles.mass = 1.0 | units.kg
        particles.x = 2.0 | units.m
        x = particles.savepoint(2.0 | units.s)

        io.write_set_to_file(x,output_file, format='amuse', version = self.store_version())

        particles_from_file = io.read_set_from_file(output_file, format='amuse', version = self.store_version(), close_file = True)
        particles_in_memory = particles_from_file.copy()

        self.assertAlmostRelativeEquals(particles_in_memory.mass[1], 1.0 | units.kg)

        particles_in_memory.savepoint(4.0 | units.s)

        self.assertAlmostRelativeEquals(particles_from_file.mass[2], 1.0 | units.kg)
    
    def test5(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test5_grid"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = self.store_factory()(output_file)
        
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        instance.store_grid(p)

        loaded_grid = instance.load_grid()

        self.assertEqual(loaded_grid.shape, shape)
        
        loaded_mass_in_kg = loaded_grid.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEqual(expected, actual)
        
        instance.close()
    
    def test6(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test6_grid"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = self.store_factory()(output_file)
        
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s
        
        instance.store_grid(p.savepoint(1| units.Myr))

        loaded_grid = self.get_version_in_store(instance.load_grid())
        
        self.assertAlmostRelativeEquals(loaded_grid.get_timestamp(), 1| units.Myr)
        
        self.assertEqual(loaded_grid.shape, shape)
        self.assertEqual(loaded_grid[0][0][0].mass, 0 | units.kg)
        self.assertAlmostRelativeEquals(loaded_grid[...,1,1].mass, p[...,1,1].mass)

        instance.close()
        
    def test7(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test7"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = self.store_factory()(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        instance.store(p)
        instance.close()
    
        instance = self.store_factory()(output_file)
        loaded_particles =  self.get_version_in_store(instance.load())

        loaded_particles.mass = [x * 3.0 for x in range(number_of_particles)] | units.kg
        previous_mass_in_kg = [x * 3.0 for x in range(number_of_particles)]
        instance.close()
        
        instance = self.store_factory()(output_file)
        loaded_particles = instance.load()
        loaded_mass_in_kg = loaded_particles.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEqual(expected, actual)
        
        instance.close()
        
        instance = self.store_factory()(output_file)
        loaded_particles = self.get_version_in_store(instance.load())
        loaded_particles[2].mass = 44 | units.kg
        instance.close()
        instance = self.store_factory()(output_file)
        loaded_particles = self.get_version_in_store(instance.load())
        self.assertEqual( 44 | units.kg, loaded_particles[2].mass)
        instance.close()
        
    def test8(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test8"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        instance = self.store_factory()(output_file)
        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = numpy.asarray([x * 2.0 for x in range(number_of_particles)])

        instance.store(p.savepoint(1 | nbody_system.time))
        instance.close()

        instance = self.store_factory()(output_file)
        loaded_particles = instance.load()
        self.assertAlmostRelativeEquals(p.mass[1], 2.0)
        instance.close()
    

    def test9(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test9"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        io.write_set_to_file(
            p,
            output_file,
            "hdf5", 
            timestamp = 2 | units.Myr,
            scale = 1 | units.kg,
            version = self.store_version(),
            close_file = True
        )
        
        loaded_particles = io.read_set_from_file(
            output_file, 
            "hdf5",
            version = self.store_version()
        )
        loaded_particles = self.get_version_in_store(loaded_particles)
        
        
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)

    def test10(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test10"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p, output_file, "hdf5",  timestamp = 2 | units.Myr, scale = 1 | units.kg, version = self.store_version())
        
        loaded_particles = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        loaded_particles = self.get_version_in_store(loaded_particles)
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)

    def test11(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test11"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p.savepoint(timestamp = 2 | units.Myr, scale = 1 | units.kg), output_file, "hdf5", version = self.store_version())
        
        loaded_particles = io.read_set_from_file(output_file, "hdf5", version = self.store_version())        
        loaded_particles = self.get_version_in_store(loaded_particles)
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)
        
    
    def test12(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test12"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        io.write_set_to_file(p.savepoint(timestamp = 2 | units.Myr, scale = 1 | units.kg), output_file, "hdf5", version = self.store_version())
        
        loaded_particles = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        loaded_particles = self.get_version_in_store(loaded_particles)
        a = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(a.timestamp, 2 | units.Myr)
        self.assertAlmostRelativeEquals(a.scale, 1 | units.kg)
    
    def test13(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test13"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[1,2,3],[4,5,6]] | units.km
       

        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        self.assertEqual(loaded[0].md, [1,2,3] | units.km)
        self.assertEqual(loaded[1].md, [4,5,6] | units.km)
        
        self.assertEqual(loaded.md[0], [1,2,3] | units.km)
        self.assertEqual(loaded.md[1], [4,5,6] | units.km)
         
        self.assertEqual(self.get_version_in_store(loaded)[0].md, [1,2,3] | units.km)
        #self.assertEquals(self.get_version_in_store(loaded)[0].md,  [7,8,9] | units.km)
    
      
    def test14(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test14"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        self.assertEqual(loaded[0].md, [[1,3],[2,4],[3,5]])
        self.assertEqual(loaded[1].md, [[4,6],[5,7],[6,8]])
        
        self.assertEqual(loaded.md[0], [[1,3],[2,4],[3,5]])
        self.assertEqual(loaded.md[1], [[4,6],[5,7],[6,8]])
    
    def test15(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test15"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        processor = self.store_factory()(output_file, True, open_for_writing = True)
        loaded = processor.load()
        self.get_version_in_store(loaded)[0].md = [[3,1],[3,4],[5,2]] 
        self.assertEqual(self.get_version_in_store(loaded)[0].md, [[3,1],[3,4],[5,2]] )
        processor.close()
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        self.assertEqual(loaded[0].md, [[3,1],[3,4],[5,2]])
        self.assertEqual(loaded[1].md, [[4,6],[5,7],[6,8]])
        
    def test16(self):
        import h5py
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test16"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        stars[0].nn = stars[1]
        stars[1].nn = stars[0]
       

        self.assertEqual(stars[0].nn,stars[1])
        self.assertEqual(stars[1].nn,stars[0])
        
        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        processor = self.store_factory()(output_file, True, open_for_writing = True)
        loaded = processor.load()
        
        self.assertEqual(loaded[0].nn,loaded[1])
        self.assertEqual(loaded[1].nn,loaded[0])
        self.assertEqual(loaded[0].nn.key,stars[1].key)
        self.assertEqual(loaded[0].nn.x,stars[1].x)
        
    
    def test17(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test17"+self.store_version()+".hdf5")
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
        
        io.write_set_to_file(binaries, output_file, "hdf5", version = self.store_version())
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertEqual(loaded[0].child1.key,stars[0].key)
        self.assertEqual(loaded[1].child1.key,stars[2].key)
        
    def test18(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test18"+self.store_version()+".hdf5")
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
        
        io.write_set_to_file(binaries, output_file, "hdf5", version = self.store_version())
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertEqual(loaded[0].child1.key,stars[0].key)
        self.assertEqual(loaded[0].child2.key,stars[1].key)
        self.assertEqual(loaded[0].child1,stars[0])
        self.assertEqual(loaded[0].child2,stars[1])
        self.assertEqual(loaded[0].child1.parent,loaded[0])
        self.assertEqual(loaded[0].child2.parent,loaded[0])
    
    def test19(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test19"+self.store_version()+".hdf5")
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
        
        io.write_set_to_file([binaries,stars], output_file, "hdf5", names = ['binaries', 'children'], version = self.store_version())
        
        loader_binaries, loaded_stars = io.read_set_from_file(output_file, "hdf5", names = ['binaries', 'children'], version = self.store_version())
        
        self.assertEqual(loader_binaries[0].child1.key,stars[0].key)
        self.assertEqual(loader_binaries[0].child2.key,stars[1].key)
        self.assertEqual(loader_binaries[0].child1,loaded_stars[0])
        self.assertEqual(loader_binaries[0].child2,loaded_stars[1])
        self.assertEqual(loader_binaries[0].child1.parent,loader_binaries[0])
        self.assertEqual(loader_binaries[0].child2.parent,loader_binaries[0])

    def test20(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test20"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        shape = 10, 10, 10
        p = Grid(*shape)
        p.mass = (numpy.asarray([x * 2.0 for x in range(p.size)])).reshape(shape)

        io.write_set_to_file(p, output_file, "hdf5", version = self.store_version())
        
        loaded = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        loaded = self.get_version_in_store(loaded)
         
        self.assertAlmostRelativeEquals(p.mass[0][1][2], 24)
        self.assertAlmostRelativeEquals(p[0][1][2].mass, 24)
        
    def test21(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test21"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[[1,3],[2,4],[3,5]],[[4,6],[5,7],[6,8]]]
       

        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
       
        loaded = io.read_set_from_file(output_file, "hdf5", close_file = True, version = self.store_version(), copy_history=False)
        
        self.assertEqual(loaded[0].md, [[1,3],[2,4],[3,5]])
        self.assertEqual(loaded[1].md, [[4,6],[5,7],[6,8]])
        
        previous = loaded.previous_state()
        self.assertEqual(previous, None)


    def test22(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test22"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars.x = 1.0  | units.km
       
        overlay = ParticlesOverlay(stars)
        overlay.y = 2.0  | units.km


        io.write_set_to_file(overlay, output_file, "hdf5", version = self.store_version())
       
        loaded = io.read_set_from_file(output_file, "hdf5", close_file = True, version = self.store_version())
        
        self.assertEqual(loaded[0].x,  1.0 | units.km)
        self.assertEqual(loaded[1].y,  2.0 | units.km)
        
    
    def test23(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test23"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Particles()

        io.write_set_to_file(particles, output_file, format='amuse', version = self.store_version())

        particles_from_file = io.read_set_from_file(output_file, format='amuse', version = self.store_version())
        
        self.assertEqual(len(particles_from_file), 0)
    
    def test24(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test24"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Particles(20)
#~        particles = Particles(1000000) # For testing memory usage
        particles.mass = 1.0 | units.kg
        particles.position = [0.0, 0.0, 0.0] | units.m
        
        for i in range(10):
            particles.position += [1.0, 2.0, 3.0] | units.m
            io.write_set_to_file(particles, output_file, format='amuse', version = self.store_version(), append_to_file=True)
        
        particles_from_file = io.read_set_from_file(output_file, format='amuse', version = self.store_version())
        os.remove(output_file)
        self.assertEqual(len(list(particles_from_file.history)), 10)
        for i, snap in enumerate(particles_from_file.history):
            self.assertEqual(len(snap), 20)
            self.assertEqual((i+1) * ([1.0, 2.0, 3.0] | units.m), snap.center_of_mass())
            
    def test25(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test25"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Particles(20)
#~        particles = Particles(1000000) # For testing memory usage
        particles.mass = 1.0 | units.kg
        particles.position = [0.0, 0.0, 0.0] | units.m
        
        for i in range(10):
            particles.position += [1.0, 2.0, 3.0] | units.m
            io.write_set_to_file(particles, output_file, format='amuse', version = self.store_version(), append_to_file=True)
        
        particles_from_file = io.read_set_from_file(output_file, format='amuse', version = self.store_version(), copy_history = True, close_file = True)
        
        history = list(particles_from_file.history)
        self.assertEqual(len(history), 10)
        self.assertFalse("HDF" in str(type(history[1]._private.attribute_storage)))
        for i, snap in enumerate(particles_from_file.history):
            self.assertEqual(len(snap), 20)
            self.assertEqual((i+1) * ([1.0, 2.0, 3.0] | units.m), snap.center_of_mass())
        os.remove(output_file)
        
    def test26(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test26"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        p = Grid(20)
        
        io.write_set_to_file(p, output_file, format='amuse', version = self.store_version())
        
        os.remove(output_file)

    def test27(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test27"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
        instance = self.store_factory()(output_file)
        
        shape = 10,
        p = Grid(*shape)
        p.mass = ([x * 2.0 for x in range(p.size)] | units.kg).reshape(shape)
        p.model_time = 2.0 | units.s

        instance.store_grid(p)

        loaded_grid = instance.load_grid()

        self.assertEqual(loaded_grid.shape, shape)
        
        loaded_mass_in_kg = loaded_grid.mass.value_in(units.kg)
        previous_mass_in_kg = p.mass.value_in(units.kg)
        for expected, actual in zip(previous_mass_in_kg, loaded_mass_in_kg):
            self.assertEqual(expected, actual)
        
        instance.close()

    def test28(self):
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test28"+self.store_version()+".hdf5")
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
        
        io.write_set_to_file(binaries, output_file, "hdf5", version = self.store_version())
        
        with io.read_set_from_file(output_file, "hdf5", version = self.store_version(), return_context = True) as loaded:
        
            self.assertEqual(loaded[0].child1.key,stars[0].key)
            self.assertEqual(loaded[0].child2.key,stars[1].key)
            self.assertEqual(loaded[0].child1,stars[0])
            self.assertEqual(loaded[0].child2,stars[1])
            self.assertEqual(loaded[0].child1.parent,loaded[0])
            self.assertEqual(loaded[0].child2.parent,loaded[0])

    def test29(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test29"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Grid(10)
        particles.attribute1 = "normal"
        particles.attribute2 = "unicode"

        io.write_set_to_file(particles,output_file, format='amuse', version = self.store_version())
        
        output = io.read_set_from_file(output_file, format='amuse')
        self.assertEqual(output[0].attribute2, "unicode")
        
    def test30(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test30"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles = Particles(10)
        particles.mass = 10 | units.kg

        io.write_set_to_file(particles,output_file, format='amuse', version = self.store_version())
        
        output = io.read_set_from_file(output_file, format='amuse', allow_writing = True)
        output.new_attribute= 2 * output.mass

        self.assertEqual(output.new_attribute, 2 * output.mass)
        output = next(output.iter_history())
        output.new_attribute= 2 * output.mass

        self.assertEqual(output.new_attribute, 2 * output.mass)
        
    def test31(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test31"+self.store_version()+".h5")
        if os.path.exists(output_file):
            os.remove(output_file)

        p=Particles(10)
        p.mass=numpy.arange(10)
        q=Particles(3)
        q[0].link=p[1]
        q[1].link=p[3]
        q[2].link=p[4]
        io.write_set_to_file(q,output_file,"amuse",version=self.store_version())
        r=io.read_set_from_file(output_file,"amuse")
        self.assertEqual(len(r.link.shape),1)
        self.assertEqual(r[0].link.mass,1)
        self.assertEqual(r[1].link.mass,3)
        self.assertEqual(r[2].link.mass,4)



class TestStoreHDFV1(amusetest.TestCase, _AbstractTestStoreHDF):
    
    def store_factory(self):
        return store_v1.StoreHDF
    
    def store_version(self):
        return '1.0'
        
    def get_version_in_store(self, container):
        return container.previous_state()
    
class TestStoreHDFV2(amusetest.TestCase, _AbstractTestStoreHDF):
    
    def store_factory(self):
        return store_v2.StoreHDF
    
    def store_version(self):
        return '2.0'
        
    def get_version_in_store(self, container):
        return container
        
    def test52(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test52"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        gas = Particles(1)
        gas[0].y = 1.0 | units.km
        
        stars[0].gas = gas
        
        
        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        
        loaded_stars = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertEqual(loaded_stars[0].gas[0].key,gas[0].key)
        self.assertAlmostRelativeEquals(loaded_stars[0].gas[0].y,1.0 | units.km)
        self.assertEqual(loaded_stars[1].gas, None)

        
    def test53(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test53"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        gas = Grid(2,3)
        gas.y = [[1.0, 2.0, 3.0], [4.0,5.0,6.0]] | units.km
        
        stars[0].gas = gas
        
        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        loaded_stars = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertAlmostRelativeEquals(loaded_stars[0].gas[0][0].y,1.0 | units.km)
        self.assertAlmostRelativeEquals(loaded_stars[0].gas[0][2].y,3.0 | units.km)
        self.assertEqual(loaded_stars[1].gas, None)

    def test54(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test23"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        gas = Grid(2,3)
        gas.y = [[1.0, 2.0, 3.0], [4.0,5.0,6.0]] | units.km
        
        stars[0].gas = gas[0][0]
        stars[1].gas = gas[1][0]
        
        self.assertAlmostRelativeEquals(stars[0].gas.y,1.0 | units.km)
        self.assertAlmostRelativeEquals(stars[1].gas.y,4.0 | units.km)
        
        io.write_set_to_file(stars, output_file, "hdf5", version = self.store_version())
        loaded_stars = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertAlmostRelativeEquals(loaded_stars[0].gas.y,1.0 | units.km)
        self.assertAlmostRelativeEquals(loaded_stars[1].gas.y,4.0 | units.km)

    def test54b(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test24"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        stars = Particles(2)
        stars[0].x = 1.0  | units.km
        stars[1].x = 2.0  | units.km
        
        gas = Grid(2,3)
        gas.y = [[1.0, 2.0, 3.0], [4.0,5.0,6.0]] | units.km
        gas[0][0].particle = stars[0]
        gas[1][2].particle = stars[0]
        gas[1][0].particle = stars[1]
                
        #stars[0].gas = gas[0][0]
        #stars[1].gas = gas[1][0]
        
        self.assertAlmostRelativeEquals(gas[0][0].particle.x,1.0 | units.km)
        self.assertAlmostRelativeEquals(gas[1][2].particle.x,1.0 | units.km)
        self.assertAlmostRelativeEquals(gas[1][0].particle.x,2.0 | units.km)
        
        io.write_set_to_file(gas, output_file, "hdf5", version = self.store_version())
        loaded_gas = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertAlmostRelativeEquals(loaded_gas[0][0].particle.x, 1.0 | units.km)
        self.assertAlmostRelativeEquals(loaded_gas[1][2].particle.x, 1.0 | units.km)
        self.assertAlmostRelativeEquals(loaded_gas[1][0].particle.x, 2.0 | units.km)
        self.assertEqual(id(loaded_gas[0][0].particle.get_containing_set()), id(loaded_gas[1][2].particle.get_containing_set()))
        
        self.assertEqual(loaded_gas.particle.shape, loaded_gas.shape)
                
        gas_copy = loaded_gas.copy()

        gas_copy[0][0].particle.x = 3 | units.km
        self.assertAlmostRelativeEquals(loaded_gas[0][0].particle.x, 1.0 | units.km)
        self.assertAlmostRelativeEquals(gas_copy[0][0].particle.x, 3.0 | units.km)
        self.assertAlmostRelativeEquals(gas_copy[1][2].particle.x, 3.0 | units.km)
        
    
    def test55(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test25"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        particles1 = Particles(2)
        particles1[0].x = 1.0  | units.km
        particles1[1].x = 2.0  | units.km
        particles2 = Particles(2)
        particles2[0].x = 3.0  | units.km
        particles2[1].x = 4.0  | units.km
        
        particles_superset = ParticlesSuperset([particles1, particles2])
    
        
        self.assertAlmostRelativeEquals(particles_superset[0].x,1.0 | units.km)
        self.assertAlmostRelativeEquals(particles_superset[2].x,3.0 | units.km)
        
        io.write_set_to_file(particles_superset, output_file, "hdf5", version = self.store_version())
        loaded_particles_superset = io.read_set_from_file(output_file, "hdf5", version = self.store_version())
        
        self.assertAlmostRelativeEquals(loaded_particles_superset[0].x,1.0 | units.km)
        self.assertAlmostRelativeEquals(loaded_particles_superset[2].x,3.0 | units.km)
        
    

    def test56(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test26"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        number_of_particles = 10
        p = Particles(number_of_particles)
        p.mass = [x * 2.0 for x in range(number_of_particles)] | units.kg
        p.model_time = 2.0 | units.s

        gas = Grid(2,3)
        gas.y = [[1.0, 2.0, 3.0], [4.0,5.0,6.0]] | units.km
        
        io.write_set_to_file(
            p,
            output_file,
            "hdf5", 
            particle = p[1],
            particles = p,
            gridpoint = gas[0][0],
            grid = gas,
            version = self.store_version()
        )
        
        loaded_particles = io.read_set_from_file(
            output_file, 
            "hdf5",
            version = self.store_version(), copy_history=False, close_file=False
        )
        
        attributes = loaded_particles.collection_attributes
        self.assertAlmostRelativeEquals(attributes.particle.mass, loaded_particles[1].mass)
        self.assertEquals(attributes.particle.key, loaded_particles[1].key)
        self.assertEquals(id(attributes.particles), id(loaded_particles))
        self.assertAlmostRelativeEquals(attributes.gridpoint.y, 1.0 | units.km)
        self.assertAlmostRelativeEquals(attributes.grid[0][0].y, 1.0 | units.km)

    def test57(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test26"+self.store_version()+".hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)

        x = Particles(2)
        x.sub = None
        x[0].sub = Particles(2)
        io.write_set_to_file(x, output_file,"amuse", version=self.store_version())
        y = io.read_set_from_file(output_file,"amuse")
        
        self.assertEqual(len(x), len(y))
        self.assertEqual(len(x[0].sub), len(y[0].sub))
        self.assertTrue(y[1].sub==x[1].sub==None)

        os.remove(output_file)
    def test58(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test58"+self.store_version()+".h5")
        if os.path.exists(output_file):
            os.remove(output_file)

        x = Particles(keys = (1,2))
        x.mass = [1,2] | units.kg
        y = Particles(keys = (11,12,13,14))
        y.mass = [40,50,60,70] | units.kg
        
        x.sub = None
        x[0].sub = y[0:2]
        y[0].sub = y[2:]
        
        io.write_set_to_file(x, output_file,"amuse", version=self.store_version())
        z = io.read_set_from_file(output_file,"amuse")
        
        self.assertEqual(len(x), len(z))
        self.assertEqual(len(x[0].sub), len(z[0].sub))
        self.assertEqual(x[0].sub[0].sub.key, (13,14))
        self.assertEqual(z[0].sub[0].sub.key, (13,14))
        self.assertEqual(id(x[0].sub._original_set()), id(x[0].sub[0].sub._original_set()))
        # no more subsets, could this be fixed, would be very nice to do so
        # self.assertEqual(id(z[0].sub._original_set()), id(z[0].sub[0].sub._original_set()))
        self.assertTrue(z[1].sub==x[1].sub==None)

        os.remove(output_file)

    def test59(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test59"+self.store_version()+".h5")
        if os.path.exists(output_file):
            os.remove(output_file)

        g=Grid(10,10)

        p=Particles(2)

        p[0].sub=g[0:5]

        io.write_set_to_file(p, output_file,"amuse", version=self.store_version())
        z = io.read_set_from_file(output_file,"amuse")
        
        os.remove(output_file)
    def test60(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test60"+self.store_version()+".h5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        x = Particles(2)
        x.mass = 1 | units.kg
        y = Particles(2)
        y.mass = 2 | units.kg
        io.write_set_to_file((x,y), output_file,"amuse", version=self.store_version(), names = ("x", "y"))
        z = io.read_set_from_file(output_file,"amuse")
        
        os.remove(output_file)
        self.assertTrue("x" in z)
        self.assertTrue("y" in z)
        self.assertAlmostRelativeEquals(z["x"].mass, 1 | units.kg)
        self.assertAlmostRelativeEquals(z["y"].mass, 2 | units.kg)


    def test61(self):
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "test61"+self.store_version()+".h5")
        if os.path.exists(output_file):
            os.remove(output_file)

       
        x = Particles(10)
        x.mass = 1 | units.kg
        y = Particles(20)
        y.id = list(range(20))
        x[0].y = y[5:7]

        io.write_set_to_file(x, output_file,"amuse", version=self.store_version())
        z = io.read_set_from_file(output_file,"amuse")
        self.assertEqual(x[0].y[0].id, 5)
        self.assertEqual(z[0].y[0].id, 5)

        os.remove(output_file)

        
        
