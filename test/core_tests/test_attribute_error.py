import numpy
import math
import os

from amuse.test import amusetest
from amuse.io import store
from amuse.units import units
from amuse.units import nbody_system
from amuse.datamodel import *
from amuse.datamodel import incode_storage

class TestAttributeError(amusetest.TestCase):
    
    def new_code_particles(self):
        class Code(object):
            def __init__(self):
                # mass
                self.data = []
                self.get_mass_called = False
                self.set_mass_called = False
                self.number_of_particles = 0
                
            def get_number_of_particles(self):
                return  self.number_of_particles
                
            def get_mass(self,index):
                self.get_mass_called = True
                data_to_return = [self.data[i] for i in index]
                return units.kg(data_to_return)
                
            def set_mass(self,index,mass):
                self.set_mass_called = True
                pass
                
            def new_particle(self, mass):
                mass = mass.value_in(units.kg)
                self.data = mass
                self.number_of_particles = len(self.data)
                return [i for i in range(len(mass))]
                
        code = Code()
        storage = incode_storage.InCodeAttributeStorage(
            code,
            incode_storage.NewParticleMethod(code.new_particle,("mass",)),
            None,
            code.get_number_of_particles,
            [],
            [
                incode_storage.ParticleGetAttributesMethod(code.get_mass,("mass",)),
            ],
            name_of_the_index = "index"
        )
        
        return Particles(storage = storage)
    
    def test1(self):
        print "Test1: Should get error when accessing non-existent attributes (InMemoryAttributeStorage)."
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])

        instances = [particles, particle, subset, superset]
        classes = [Particles, Particle, ParticlesSubset, ParticlesSuperset]
        lengths = [4, 1, 2, 5]
        for i, x in enumerate(instances):
            self.assertTrue(isinstance(x, classes[i]))
            self.assertEquals(len(x.as_set()), lengths[i])
            self.assertRaises(AttributeError, lambda: x.bogus, expected_message = 
                "You tried to access attribute 'bogus' but this attribute is not defined for this set.")
                
    def test2(self):
        print "Test2: Should get error when accessing non-existent attributes (in Legacy code storage)."
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])
        superset.mass = 1.0 | units.MSun
        superset.radius = 1.0 | units.RSun
        for i, x in enumerate(superset):
            x.mass = 2.0 | units.kg
        instances = [particles, particle.as_set(), subset, superset]
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        for i, x in enumerate(instances):
            code_particles = self.new_code_particles()
            code_particles.add_particles(x)
            self.assertRaises(AttributeError, lambda: code_particles.bogus, expected_message = 
                "You tried to access attribute 'bogus' but this attribute is not defined for this set.")
    
    def test3(self):
        print "Test3: Should get error when accessing non-existent attributes (HDF5 storage)."
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])
        superset.mass = 1.0 | units.MSun
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "attr_test.hdf5")
        
        instances = [particles, particle, subset]
        classes = [Particles, Particle, ParticlesSubset]
        lengths = [4, 1, 2]
        for i, x in enumerate(instances):
            self.assertTrue(isinstance(x, classes[i]))
            self.assertEquals(len(x.as_set()), lengths[i])
            if os.path.exists(output_file):
                os.remove(output_file)
            HDFstorage = store.StoreHDF(output_file)
            if isinstance(x, Particle): x = x.as_set()
            x.model_time = 2.0 | units.s
            HDFstorage.store(x)
            loaded_particles = HDFstorage.load()
            self.assertRaises(AttributeError, lambda: loaded_particles.bogus, expected_message = 
                "You tried to access attribute 'bogus' but this attribute is not defined for this set.")
            HDFstorage.close()
            del HDFstorage
    
    def bogus_func(self, x):
        x.mass = 1.0
    
    def xtest4(self):
        print "Test4: Should get error when setting attributes with non-quantities (InMemoryAttributeStorage)."
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])

        instances = [particles, particle, subset, superset]
        classes = [Particles, Particle, ParticlesSubset, ParticlesSuperset]
        lengths = [4, 1, 2, 5]
        for i, x in enumerate(instances):
            self.assertTrue(isinstance(x, classes[i]))
            self.assertEquals(len(x.as_set()), lengths[i])
            self.assertRaises(AttributeError, self.bogus_func, x, expected_message = 
                "Can only assign quantities or other particles to an attribute.")
    
    def xtest5(self):
        print "Test5: Should get error when setting attributes with non-quantities (in code storage)."        
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])
        superset.mass = 1.0 | units.MSun
        for i, x in enumerate(superset):
            x.mass = 2.0 | units.kg
        instances = [particles, particle.as_set(), subset, superset]
        classes = [Particles, Particle, ParticlesSubset, ParticlesSuperset]
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        for i, x in enumerate(instances):
            code_particles = self.new_code_particles()
            code_particles.add_particles(x)
            #self.assertRaises(AttributeError, self.bogus_func, code_particles, expected_message = 
            #    "Can only assign quantities or other particles to an attribute.")
    
    def xtest6(self):
        print "Test6: Should get error when setting attributes with non-quantities (HDF5 storage)."
        particles = Particles(4)
        particle  = Particle()
        subset    = particles[:2]
        superset  = ParticlesSuperset([particles, particle.as_set()])
        superset.mass = 1.0 | units.MSun
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "attr_test.hdf5")
        instances = [particles, particle, subset]
        classes = [Particles, Particle, ParticlesSubset]
        lengths = [4, 1, 2]
        for i, x in enumerate(instances):
            self.assertTrue(isinstance(x, classes[i]))
            self.assertEquals(len(x.as_set()), lengths[i])
            if os.path.exists(output_file):
                os.remove(output_file)
            HDFstorage = store.StoreHDF(output_file)
            if isinstance(x, Particle): x = x.as_set()
            x.model_time = 2.0 | units.s
            HDFstorage.store(x)
            loaded_particles = HDFstorage.load()
            self.assertRaises(AttributeError, self.bogus_func, loaded_particles, expected_message = 
                "Can only assign quantities or other particles to an attribute.")
            HDFstorage.close()
            del HDFstorage
    
   
