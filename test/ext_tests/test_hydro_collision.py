import numpy

from amuse.units import units, nbody_system
from amuse.units.quantities import zero
from amuse.datamodel import Particles, Particle
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestWithMPI
from amuse.plot import pynbody_column_density_plot, HAS_PYNBODY

from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from amuse.community.seba.interface import SeBa
from amuse.community.evtwin.interface import EVtwin

from amuse.ext.plotting_hydro import new_plotting_hydrodynamics_code
from amuse.ext.hydro_collision import StellarEncounterInHydrodynamics

class TestStellarEncounterInHydrodynamics(TestWithMPI):
    
    def new_colliders(self):
        colliders = Particles(2)
        colliders.mass = [5, 2] | units.MSun
        colliders.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.RSun
        colliders.velocity = [[0.0, 0.0, 0.0], [0.0, 2000.0, 0.0]] | units.km / units.s
        colliders.move_to_center()
        return colliders
    
    def test1(self):
        print "Test collect_required_attributes"
        in_memory = self.new_colliders()
        gravity = BHTree(nbody_system.nbody_to_si(1|units.MSun, 1.0|units.RSun))
        gravity.particles.add_particles(in_memory)
        stellar = SeBa()
        stellar.particles.add_particles(in_memory)
        
        collision = StellarEncounterInHydrodynamics(None, None, verbose=True)
        
        self.assertFalse(hasattr(in_memory, "radius"))
        collision.collect_required_attributes(in_memory, gravity, stellar)
        self.assertTrue(hasattr(in_memory, "radius"))
        self.assertAlmostRelativeEqual(in_memory.radius.sum(), 4.2458 | units.RSun, 3)
        
        from_stellar = stellar.particles.copy()
        for attribute in ["x", "y", "z", "vx", "vy", "vz"]:
            self.assertFalse(hasattr(from_stellar, attribute))
        collision.collect_required_attributes(from_stellar, gravity, stellar)
        
        gravity.stop()
        stellar.stop()
        
        for attribute in ["x", "y", "z", "vx", "vy", "vz"]:
            self.assertTrue(hasattr(from_stellar, attribute))
        self.assertAlmostEqual(from_stellar.position, in_memory.position)
        self.assertAlmostEqual(from_stellar.velocity, in_memory.velocity)
    
    def test2(self):
        print "Test backtrack_particles"
        colliders = self.new_colliders()
        colliders.radius = [1, 2] | units.RSun
        
        collision = StellarEncounterInHydrodynamics(None, None, verbose=True)
        
        relative_position = colliders[1].position - colliders[0].position
        self.assertAlmostRelativeEqual(relative_position.length(), 1.0 | units.RSun, 7)
        
        total_energy_before = colliders.kinetic_energy() + colliders.potential_energy()
        self.assertTrue(total_energy_before > zero)
        
        collision.backtrack_particles(colliders)
        collision.kepler.stop()
        
        relative_position = colliders[1].position - colliders[0].position
        self.assertAlmostRelativeEqual(relative_position.length(), 15 | units.RSun, 3)
        
        total_energy_after = colliders.kinetic_energy() + colliders.potential_energy()
        self.assertAlmostRelativeEqual(total_energy_after, total_energy_before, 7)
    
    def test3(self):
        print "Test convert_stars"
        colliders = self.new_colliders()
        colliders.position = [[-100.0, 0.0, 0.0], [100.0, 0.0, 0.0]] | units.RSun
        stellar = EVtwin()
        stellar.particles.add_particles(colliders)
        
        collision = StellarEncounterInHydrodynamics(700, None, relax_sph_models=False, verbose=True)
        gas_particles = collision.convert_stars(colliders, stellar)
        stellar.stop()
        
        self.assertEqual(gas_particles.mass, 0.01 | units.MSun)
        self.assertTrue(numpy.all(gas_particles[:500].x < zero))
        self.assertTrue(numpy.all(gas_particles[500:].x > zero))
        
        self.assertIsOfOrder((
                gas_particles[:500].position - ([-100.0, 0.0, 0.0] | units.RSun)
            ).lengths_squared().amax().sqrt(), 1 | units.RSun)
        self.assertIsOfOrder((
                gas_particles[500:].position - ([100.0, 0.0, 0.0] | units.RSun)
            ).lengths_squared().amax().sqrt(), 1 | units.RSun)
        
        self.assertAlmostEqual(gas_particles[500:].center_of_mass_velocity().y - 
            gas_particles[:500].center_of_mass_velocity().y, 2000.0 | units.km / units.s)
    
    def test4(self):
        print "Test binary_will_collide"
        collision = StellarEncounterInHydrodynamics(None, None, verbose=True)
        collision.dynamical_timescale = zero
        collision.start_kepler(7 | units.MSun, 10 | units.RSun)
        # at periastron, close enough:
        colliders = self.new_colliders()
        colliders.radius = 1 | units.RSun
        self.assertTrue(collision.binary_will_collide(colliders[0], colliders[1]))
        # at periastron, distance too large:
        colliders.radius = 0.4 | units.RSun
        self.assertFalse(collision.binary_will_collide(colliders[0], colliders[1]))
        # at apastron, will collide at periastron:
        colliders.velocity = [[0.0, 0.0, 0.0], [0.0, 1000.0, 0.0]] | units.km / units.s
        self.assertTrue(collision.binary_will_collide(colliders[0], colliders[1]))
        # hyperbolic orbits, moving away from each other:
        colliders.position = [[0.0, 0.0, 0.0], [1.0, 100.0, 0.0]] | units.RSun
        self.assertFalse(collision.binary_will_collide(colliders[0], colliders[1]))
        # hyperbolic orbits, moving towards each other:
        colliders.velocity = [[0.0, 0.0, 0.0], [0.0, -1000.0, 0.0]] | units.km / units.s
        self.assertTrue(collision.binary_will_collide(colliders[0], colliders[1]))
        collision.kepler.stop()
    
    def test5(self):
        print "Test group_bound_particles"
        colliders = self.new_colliders()
        colliders.position = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]] | units.RSun
        colliders.velocity = [[0.0, 0.0, 0.0], [10000, 0.0, 0.0]] | units.km / units.s
        stellar = EVtwin()
        stellar.particles.add_particles(colliders)
        
        collision = StellarEncounterInHydrodynamics(7000, None, relax_sph_models=False, verbose=True, 
            star_to_sph_arguments=dict(base_grid_options=dict(type="sobol")))
        collision.dynamical_timescale = zero
        gas_particles = collision.convert_stars(colliders, stellar)
        stellar.stop()
        
        collision.hop = collision.new_hop(gas_particles)
        collision.start_kepler(7 | units.MSun, 10 | units.RSun)
        self.assertTrue(collision.encounter_is_over(gas_particles))
        collision.hop.stop()
        collision.kepler.stop()
        
        groups = collision.groups_after_encounter
        self.assertEqual(len(groups), 2)
        self.assertTrue(4500 < len(groups[0]) < 5000)
        self.assertTrue(1800 < len(groups[1]) < 2000)
        self.assertEqual(len(gas_particles - groups[0] - groups[1]), 346)
        self.assertAlmostEqual(groups[0].center_of_mass()[0], 0 | units.RSun, 1)
        self.assertAlmostEqual(groups[1].center_of_mass()[0], 1.1 | units.RSun, 0)
        self.assertIsOfOrder(groups[1].center_of_mass_velocity()[0], 10000 | units.km / units.s)
    
    def test6(self):
        print "Test handle_collision"
        position_offset = [100.0, 200.0, 300.0] | units.RSun
        velocity_offset = [10000.0, 20000.0, 30000.0] | units.km / units.s
        colliders = self.new_colliders()
        colliders.position *= 1.5 # Grazing collision
        colliders.position += position_offset
        colliders.velocity += velocity_offset
        
        class GravityCodeStub(object):
            def __init__(self, particles):
                self.particles = particles
        gravity = GravityCodeStub(colliders)
        
        stellar = EVtwin()
        stellar.particles.add_particles(colliders)
        
        collision = StellarEncounterInHydrodynamics(280,  # For speed
            Gadget2,
            relax_sph_models=False,  # For speed
            verbose=True)
        collision.extra_steps_when_encounter_is_over = 0 # For speed, no extra steps...
        collision.dynamical_timescales_per_step = 1.3 # ... but then we need to evolve a bit longer in the first step
        result = collision.handle_collision(colliders[0], colliders[1], gravity_code=gravity, stellar_evolution_code=stellar)
        stellar.stop()
        print result
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(len(result), 2)
        self.assertAlmostEqual(result.mass, [5.0, 1.9] | units.MSun, 1)
        self.assertAlmostRelativeEqual(result.center_of_mass(), position_offset, 2)
        self.assertAlmostRelativeEqual(result.center_of_mass_velocity(), velocity_offset, 2)
    
    def slowtest7(self):
        print "Test handle_collision"
        position_offset = [100.0, 200.0, 300.0] | units.RSun
        velocity_offset = [10000.0, 20000.0, 30000.0] | units.km / units.s
        colliders = self.new_colliders()
        colliders.position += position_offset
        colliders.velocity += velocity_offset
        
        class GravityCodeStub(object):
            def __init__(self, particles):
                self.particles = particles
        gravity = GravityCodeStub(colliders)
        
        stellar = EVtwin()
        stellar.particles.add_particles(colliders)
        
        collision = StellarEncounterInHydrodynamics(
            #28000, # Bit too slow, even for a slowtest... 
            2800, 
            new_plotting_hydrodynamics_code(Gadget2, 0.2|units.hour, 
                plot_function = pynbody_column_density_plot if HAS_PYNBODY else None,
                plot_function_arguments = dict(width=20|units.RSun, vmin=29, vmax=35) if HAS_PYNBODY else dict(width=20|units.RSun)), 
            hydrodynamics_arguments = dict(redirection="file", redirect_file="hydro_collision_slowtest7_gadget_out.log"),
            hydrodynamics_parameters = dict(time_limit_cpu=1|units.day, gas_epsilon=0.01|units.RSun),
            verbose=True,
            debug=False
        )
        result = collision.handle_collision(colliders[0], colliders[1], gravity_code=gravity, stellar_evolution_code=stellar)
        stellar.stop()
        print result
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(len(result), 2)
        self.assertTrue((result.mass < [5.0, 2.0] | units.MSun).all())
        self.assertAlmostRelativeEqual(result.center_of_mass(), position_offset, 2)
        self.assertAlmostRelativeEqual(result.center_of_mass_velocity(), velocity_offset, 2)
    

