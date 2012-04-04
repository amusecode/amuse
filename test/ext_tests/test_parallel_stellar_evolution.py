import numpy

from amuse.units import units
from amuse.datamodel import Particles
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase

from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.sse.interface import SSE

from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution

class TestParallelStellarEvolution(TestCase):
    
    code_factory = EVtwin
    
    def test1(self):
        print "Testing ParallelStellarEvolution initialization"
        instance = ParallelStellarEvolution(self.code_factory, number_of_workers=4)
        instance.initialize_code()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Testing ParallelStellarEvolution particles"
        instance = ParallelStellarEvolution(self.code_factory, number_of_workers=3)
        instance.initialize_code()
        instance.commit_parameters()
        particles = Particles(7)
        particles.mass = range(1, 1+len(particles)) | units.MSun
        incode = instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertAlmostEqual(incode.mass, range(1, 1+len(particles)) | units.MSun)
        print "Note that the order of instance.particles is different from the",
        print "original particle order, since particles are distributed over 3 processes"
        self.assertAlmostEqual(instance.particles.mass, [1,4,7, 2,5, 3,6] | units.MSun)
        instance.stop()
        
    def slowtest3(self):
        print "Testing ParallelStellarEvolution evolve_model"
        particles = Particles(4)
        particles.mass = range(1, 1+len(particles)) | units.MSun
        
        serial = self.code_factory()
        inserial = serial.particles.add_particles(particles)
        self.assertAlmostEqual(inserial.mass, range(1, 1+len(particles)) | units.MSun)
        serial.evolve_model(0.2 | units.Myr)
        
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=3)
        inparallel = parallel.particles.add_particles(particles)
        self.assertAlmostEqual(inparallel.mass, range(1, 1+len(particles)) | units.MSun)
        parallel.evolve_model(0.2 | units.Myr)
        self.assertEqual(parallel.model_time, 0.2 | units.Myr)
        self.assertAlmostEqual(inparallel.age, [0.2]*4 | units.Myr)
        
        self.assertEqual(inserial.luminosity, inparallel.luminosity)
        self.assertEqual(inserial.time_step, inparallel.time_step)
        self.assertEqual(inserial.temperature, inparallel.temperature)
        serial.stop()
        parallel.stop()
    
    def test4(self):
        print "Testing ParallelStellarEvolution parameters"
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=3)
        self.assertRaises(AmuseException, getattr, parallel, "parameters", expected_message=
            "Not implemented for parallel stellar evolution")
#~        parallel.parameters.metallicity = 0.01
#~        self.assertEqual(parallel.parameters.metallicity, 0.01)
        parallel.stop()
    
