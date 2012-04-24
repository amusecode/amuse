import os
import os.path
import numpy

from amuse.units import units
from amuse.datamodel import Particles
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase, get_path_to_results

from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.sse.interface import SSE

from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution


default_options = dict(must_run_threaded = False) # Always thread-safe
#~default_options = dict() # Really really parallel


class TestParallelStellarEvolution(TestCase):
    
    code_factory = SSE
    
    def test1(self):
        print "Testing ParallelStellarEvolution initialization"
        instance = ParallelStellarEvolution(self.code_factory, number_of_workers=3, **default_options)
        instance.initialize_code()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Testing ParallelStellarEvolution particles"
        instance = ParallelStellarEvolution(self.code_factory, number_of_workers=2, **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        particles = Particles(5)
        particles.mass = range(1, 1+len(particles)) | units.MSun
        incode = instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertAlmostEqual(incode.mass, range(1, 1+len(particles)) | units.MSun)
        print "Note that the order of instance.particles is different from the",
        print "original particle order, since particles are distributed over 2 processes"
        self.assertAlmostEqual(instance.particles.mass, [1,3,5, 2,4] | units.MSun)
        
        instance.stop()
        
    def slowtest3(self):
        print "Testing ParallelStellarEvolution evolve_model"
        particles = Particles(4)
        particles.mass = range(1, 1+len(particles)) | units.MSun
        
        serial = self.code_factory()
        inserial = serial.particles.add_particles(particles)
        self.assertAlmostEqual(inserial.mass, range(1, 1+len(particles)) | units.MSun)
        serial.evolve_model(0.2 | units.Myr)
        
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=3, **default_options)
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
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=3, **default_options)
        parallel.parameters.metallicity = 0.01
        self.assertEqual(parallel.parameters.metallicity, 0.01)
        for code in parallel.code_instances:
            self.assertEqual(code.parameters.metallicity, 0.01)
        parallel.stop()
    
    def test5(self):
        print "Testing ParallelStellarEvolution individual options"
        base_name = os.path.join(get_path_to_results(), "parallel_stellar_evolution_out_")
        for filename in [base_name+str(i) for i in range(3)]:
            if os.path.exists(filename):
                os.remove(filename)
        
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=3, 
            individual_options=[dict(redirect_file=base_name+str(i)) for i in range(3)], redirection="file", **default_options)
        
        for filename in [base_name+str(i) for i in range(3)]:
            self.assertTrue(os.path.exists(filename))
        
        parallel.stop()
    
    def test6(self):
        print "Testing ParallelStellarEvolution exception handling"
        self.code_factory = MESA
        if self.code_factory == MESA:
            expected_message = ("Error when calling 'evolve_for' of a 'MESA', errorcode is -12, error is 'Evolve terminated: Maximum age reached.'")
        elif self.code_factory == EVtwin:
            expected_message = ("Error when calling 'evolve_for' of a 'EVtwin', errorcode is 5, error is 'PRINTB -- age greater than limit'")
        else:
            self.skip("Skipping test: {0} has no parameter max_age_stop_condition".format(self.code_factory))
        
        parallel = ParallelStellarEvolution(self.code_factory, number_of_workers=2, **default_options)
        parallel.parameters.max_age_stop_condition = 0.1 | units.Myr
        
#~        parallel.particles.add_particles(Particles(2, mass=[1,2]|units.MSun)) # Test speed-up:
        parallel.particles.add_particles(Particles(1, mass=1|units.MSun))
        
        self.assertRaises(AmuseException, parallel.evolve_model, 1.0|units.Myr, 
            expected_message = expected_message)
        self.assertTrue((parallel.particles.age >= 0.1 | units.Myr).all())
        self.assertTrue((parallel.particles.age-parallel.particles.time_step < 0.1 | units.Myr).all())
        parallel.stop()
    
