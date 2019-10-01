# nosetests --nocapture --nologcapture -w test/codes_tests --tests=test_multiples 

from amuse.test.amusetest import TestWithMPI

import os
import sys
import numpy
import time
import math

from amuse.community.hermite.interface import Hermite
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants

from amuse import datamodel
from amuse.ic import plummer
from amuse.couple import multiples
from amuse.couple import encounters
from amuse import io

class TestSimpleMultiples(TestWithMPI):
    previous = None
    
    def new_smalln(self):
        if not self.previous is None:
            self.previous.stop()
            
        result = SmallN()
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 2001
        self.previous = result
        return result
        
    def new_kepler_si(self):
        unit_converter = nbody_system.nbody_to_si(
            1.0 | units.MSun,
            1.0 | units.AU
        )
        kepler = Kepler(unit_converter)
        kepler.initialize_code()
        return kepler
        
    def new_kepler(self):
        kepler = Kepler()
        kepler.initialize_code()
        return kepler
        
    def new_smalln_si(self):
    
        if not self.previous is None:
            self.previous.stop()
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        result = SmallN(converter)
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 2001
        return result
        
    def new_binary(self, mass1, mass2, semi_major_axis,
                   eccentricity = 0, keyoffset = -1):
        total_mass = mass1 + mass2
        mass_fraction_particle_1 = mass1 / (total_mass)
    
        if keyoffset >= 0:
            binary = datamodel.Particles(keys=list(range(keyoffset, keyoffset+2)))
        else:
            binary = datamodel.Particles(2)
            
        binary[0].mass = mass1
        binary[1].mass = mass2
    
        mu = nbody_system.G * total_mass
    
        velocity_perihelion = numpy.sqrt( mu / semi_major_axis  * ((1.0 + eccentricity)/(1.0 - eccentricity)))
        radius_perihelion = semi_major_axis * (1.0 - eccentricity)
        
        binary[0].position = ((1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0,0.0,0.0])
        binary[1].position = -(mass_fraction_particle_1 * radius_perihelion * [1.0,0.0,0.0])
    
        binary[0].velocity = ((1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0,1.0,0.0])
        binary[1].velocity = -(mass_fraction_particle_1 * velocity_perihelion * [0.0,1.0,0.0])

        return binary
        
    
    def create_binaries(self, center_of_mass_particles, mass1, mass2, semi_major_axis,
                   eccentricity = 0):
        singles_in_binaries = datamodel.Particles()
        for binary in center_of_mass_particles:
            particles_in_binary = self.new_binary(
                mass1,
                mass2,
                semi_major_axis
            )
                
            particles_in_binary.radius = semi_major_axis
            
            binary.child1 = particles_in_binary[0]
            binary.child2 = particles_in_binary[1]
            binary.mass = mass1 + mass2
           
            particles_in_binary.position += binary.position
            particles_in_binary.velocity += binary.velocity
            singles_in_binaries.add_particles(particles_in_binary)
        return center_of_mass_particles, singles_in_binaries
        
    
    def test0(self):
        code = Hermite()
        stars = datamodel.Particles(2)
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0.0, 0,0],
            [1.2, 0, 0]
        ]|nbody_system.length
        stars.velocity = [
            [0.0,0,0],
            [0,0.1, 0]
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        encounter_code.parameters.hard_binary_factor = 1
        encounter_code.small_scale_factor = 1
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        multiples_code.evolve_model(0.6|nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.binaries), 1)
        
        
    def test1(self):
        code = Hermite()
        stars = datamodel.Particles(keys = (1,2,3, 4))
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0.0,0,0],
            [0.5, 0, 0],
            [2.0, 0, 0],
            [-10.0, 0, 0],
        ]|nbody_system.length
        stars.velocity = [
            [0.0,0,0],
            [0,0.1, 0],
            [0,-0.1, 0],
            [0,0.2, 0],
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        
        multiples_code.evolve_model(0.6|nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.binaries), 1)
        
        self.assertAlmostRelativeEquals(multiples_code.particles[:-1].radius, 0.5 | nbody_system.length)
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].radius, 0.4446| nbody_system.length, 3)
        multiples_code.evolve_model(2|nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.binaries), 1)
        multiples_code.evolve_model(3|nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.particles), 2)
        self.assertEqual(len(multiples_code.binaries), 1)
        

    def test2(self):
        code = Hermite()
        stars = datamodel.Particles(keys = (1,2,3, 4))
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0.0,0,0],
            [0.5, 0, 0],
            [3, 0, 0],
            [-10, 0, 0],
        ]|nbody_system.length
        stars.velocity = [
            [0.0,0,0],
            [0,0.1, 0],
            [0.0,-0.5, 0],
            [0,0.2, 0],
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        
        multiples_code.evolve_model(3|nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)
        print(multiples_code.multiples[0].components)
        self.assertEqual(len(multiples_code.multiples[0].components), 2)
        self.assertEqual(len(multiples_code.particles), 3)
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.singles), 2)
    
    
    def test3(self):
        code = Hermite()
        particles_in_binary = self.new_binary(
            0.1 | nbody_system.mass,
            0.1 | nbody_system.mass,
            0.01 | nbody_system.length,
            keyoffset = 1
        )
        particles_in_binary.radius = 0.001 | nbody_system.length
        binary = datamodel.Particle(key = 3)
        binary.child1 = particles_in_binary[0]
        binary.child2 = particles_in_binary[1]
        binary.radius = 0.5 | nbody_system.length
        binary.mass = 0.2 | nbody_system.mass
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(particles_in_binary)
        multiples_code.binaries.add_particle(binary)
        
        self.assertEqual(len(multiples_code.singles_in_binaries), 2)
        self.assertEqual(id(multiples_code.binaries[0].child1.particles_set), id(multiples_code.singles_in_binaries))
        
        multiples_code.commit_particles()
        
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.components_of_multiples), 2)
        
        
    
    def test4(self):
        code = Hermite()
        stars = datamodel.Particles(keys = (1,2,3, 4))
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0.0,0,0],
            [0.5, 0, 0],
            [2, 0, 0],
            [-10, 0, 0],
        ]|nbody_system.length
        stars.velocity = [
            [0,0,0],
            [0,0.2, 0],
            [0,-0.2, 0],
            [0,0.3, 0],
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        stopping_condition = multiples_code.stopping_conditions.multiples_change_detection
        stopping_condition.enable()
        
        multiples_code.evolve_model(3|nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        self.assertAlmostRelativeEquals(multiples_code.model_time , 0.0075 | nbody_system.time, 4)
        self.assertEqual(len(stopping_condition.particles(0)), 1)
        self.assertEqual(len(stopping_condition.particles(1)), 0)
        
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.multiples[0].components), 2)
        self.assertEqual(len(multiples_code.particles), 3) # 1 multiples with 2 singles, plus 2 singles free
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.singles), 2)
        
        multiples_code.evolve_model(3|nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        self.assertAlmostRelativeEquals(multiples_code.model_time , 1.2195 | nbody_system.time, 4)
        self.assertEqual(len(stopping_condition.particles(0)), 1) # 1 new multiple
        self.assertEqual(len(stopping_condition.particles(1)), 1) # 1 dissolved multiple
        
        self.assertEqual(len(multiples_code.multiples[0].components), 3)
        self.assertEqual(len(multiples_code.particles), 2) # 1 multiple, plus 1 single free
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.singles), 1)
    
    def test5(self):
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        
        code = Hermite(converter)
        stars = datamodel.Particles(keys=(1,2))
        stars.mass = converter.to_si(1 | nbody_system.mass)
        stars.position = converter.to_si([
            [0,0,0],
            [1.2, 0, 0]
        ]|nbody_system.length)
        stars.velocity = converter.to_si([
            [0,0,0],
            [0,0.1, 0]
        ]|nbody_system.speed)
        stars.radius = converter.to_si(0.5 | nbody_system.length)
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler_si(),
            resolve_collision_code = self.new_smalln_si(),
            interaction_over_code = None,
            G = constants.G
        )
        encounter_code.parameters.hard_binary_factor = 1
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code,
            G = constants.G
        )
        end_time = converter.to_si(1.0|nbody_system.time)
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        multiples_code.evolve_model(end_time)
        
        self.assertEqual(len(multiples_code.particles),1) # 1 multiples with 2 singles
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.multiples[0].components), 2)
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.singles), 0)
        
        
    
    def test6(self):
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        
        code = Hermite(converter)
        stars = datamodel.Particles(keys=(1,2,3,4))
        stars.mass = converter.to_si(1 | nbody_system.mass)
        stars.position = converter.to_si([
            [0,0,0],
            [1.2, 0, 0],
            [100, 0, 0],
            [100, 1.2, 0]
        ]|nbody_system.length)
        stars.velocity = converter.to_si([
            [0,0,0],
            [0,0.1, 0],
            [0,0,0],
            [0,0,0.1],
        ]|nbody_system.speed)
        stars.radius = converter.to_si(0.5 | nbody_system.length)
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler_si(),
            resolve_collision_code = self.new_smalln_si(),
            interaction_over_code = None,
            G = constants.G
        )
        encounter_code.small_scale_factor = 1.0
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code,
            G = constants.G
        )
        multiples_code.must_handle_one_encounter_per_stopping_condition = False
        multiples_code.particles.add_particles(stars)
        multiples_code.commit_particles()
        
        stopping_condition = multiples_code.stopping_conditions.multiples_change_detection
        stopping_condition.enable()
        
        end_time = converter.to_si(3.0|nbody_system.time)
        print(end_time.as_quantity_in(units.Myr))
        multiples_code.evolve_model(end_time)
        self.assertTrue(stopping_condition.is_set())
        print(multiples_code.model_time.as_quantity_in(units.Myr))
        self.assertAlmostRelativeEquals(multiples_code.model_time , 7.99844 | units.Myr, 4)
        self.assertEqual(len(stopping_condition.particles(0)), 2)
        self.assertEqual(len(stopping_condition.particles(1)), 0)
        
        self.assertEqual(len(multiples_code.particles), 2)             # 1 multiples with 2 singles
        self.assertEqual(len(multiples_code.multiples), 2)
        self.assertEqual(len(multiples_code.binaries), 2)
        self.assertEqual(len(multiples_code.multiples[0].components), 2)
        self.assertEqual(len(multiples_code.multiples[1].components), 2)
        self.assertEqual(len(multiples_code.singles), 0)
        self.assertEqual(len(multiples_code.all_singles), 4)

    
    def test7(self):
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        
        code = Hermite(converter)
        stars = datamodel.Particles(keys=(1,2))
        stars.mass = converter.to_si(1 | nbody_system.mass)
        stars.position = converter.to_si([
            [0,0,0],
            [1.1, 0, 0],
        ]|nbody_system.length)
        stars.velocity = converter.to_si([
            [0,0,0],
            [-0.5,1.5, 0],
        ]|nbody_system.speed)
        stars.radius = converter.to_si(0.55 | nbody_system.length)
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler_si(),
            resolve_collision_code = self.new_smalln_si(),
            interaction_over_code = None,
            G = constants.G
        )
        encounter_code.small_scale_factor = 1.0
        encounter_code.parameters.hard_binary_factor = 1
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code,
            G = constants.G
        )
        multiples_code.must_handle_one_encounter_per_stopping_condition = False
        multiples_code.singles.add_particles(stars)
        multiples_code.commit_particles()
        
        stopping_condition = multiples_code.stopping_conditions.encounter_detection
        stopping_condition.enable()
        
        end_time = converter.to_si(3.0|nbody_system.time)
        print(end_time.as_quantity_in(units.Myr))
        multiples_code.evolve_model(end_time)
        self.assertTrue(stopping_condition.is_set())
        print(multiples_code.model_time.as_quantity_in(units.Myr))
        #self.assertAlmostRelativeEquals(multiples_code.model_time , 5.96955 | units.Myr, 4)
        self.assertEqual(len(stopping_condition.particles(0)), 1)
        model = stopping_condition.particles(0)[0]
        
        self.assertEqual(len(model.particles_before_encounter), 2)
        self.assertEqual(len(model.particles_after_encounter), 2)
        
        
        before = model.particles_before_encounter
        after = model.particles_after_encounter
        
        self.assertAlmostRelativeEquals(before.center_of_mass(), after.center_of_mass(), 7)
        self.assertAlmostRelativeEquals(before.center_of_mass_velocity(), after.center_of_mass_velocity(), 7)
        
        total_energy_before = before.kinetic_energy() + before.potential_energy(G=constants.G)
        total_energy_after = after.kinetic_energy() + after.potential_energy(G=constants.G)
        
        self.assertAlmostRelativeEquals(total_energy_before, total_energy_after, 7)
        
    
    def test8(self):
        code = Hermite()
        particles_in_binary = self.new_binary(
            0.1 | nbody_system.mass,
            0.1 | nbody_system.mass,
            0.01 | nbody_system.length,
            keyoffset = 1
        )
        particles_in_binary.radius = 0.001 | nbody_system.length
        binary = datamodel.Particle(key = 3)
        binary.child1 = particles_in_binary[0]
        binary.child2 = particles_in_binary[1]
        binary.radius = 0.5 | nbody_system.length
        binary.mass = 0.2 | nbody_system.mass
        binary.position = [0.0,0.0,0.0] | nbody_system.length
        binary.velocity = [0.0,0.0,0.0] | nbody_system.speed
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
            interaction_over_code = None
        )
        encounter_code.parameters.hard_binary_factor = 1
        encounter_code.small_scale_factor = 1
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(particles_in_binary)
        multiples_code.binaries.add_particle(binary)
        multiples_code.must_handle_one_encounter_per_stopping_condition = False
        
        field_particle = datamodel.Particle(key = 4)
        field_particle.mass = 0.5  | nbody_system.mass
        field_particle.radius = 0.1 | nbody_system.length
        field_particle.position = [0.0,0.2,0.0]| nbody_system.length
        field_particle.velocity = [0.0,0.0,0.0] | nbody_system.speed
        
        multiples_code.singles.add_particle(field_particle)
        
        self.assertEqual(len(multiples_code.singles_in_binaries), 2)
        self.assertEqual(id(multiples_code.binaries[0].child1.particles_set), id(multiples_code.singles_in_binaries))
        
        multiples_code.commit_particles()
        multiples_code.multiples.radius =  0.5 | nbody_system.length
        initial_energy = multiples_code.get_total_energy()
        
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.components_of_multiples), 2)
        self.assertEqual(len(multiples_code.particles), 2)
        
        stopping_condition = multiples_code.stopping_conditions.encounter_detection
        stopping_condition.enable()
        
        singles = datamodel.Particles()
        singles.add_particles(particles_in_binary)
        singles.add_particle(field_particle)
        
        
        singles_energy = singles.kinetic_energy() + singles.potential_energy(G=nbody_system.G)
        self.assertAlmostRelativeEquals(initial_energy, singles_energy, 3)
        
        
        multiples_code.evolve_model(2 |nbody_system.time)
        
        final_energy = multiples_code.get_total_energy()
        self.assertTrue(stopping_condition.is_set())
        self.assertAlmostRelativeEquals(initial_energy, final_energy, 7)
        
    
    def test9(self):
        code = Hermite()
        
        particles_in_binary = self.new_binary(
            0.1 | nbody_system.mass,
            0.1 | nbody_system.mass,
            0.01 | nbody_system.length,
            keyoffset = 1
        )
        particles_in_binary.radius = 0.001 | nbody_system.length
        
        binary = datamodel.Particle(key = 3)
        binary.child1 = particles_in_binary[0]
        binary.child2 = particles_in_binary[1]
        binary.radius = 0.5 | nbody_system.length
        binary.mass = 0.2 | nbody_system.mass
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        
        
        others = datamodel.Particles(key = [4,5,6])
        for i in range(3):
            others[i].position = [i, 0, 0] | nbody_system.length
            others[i].velocity = [0, 0, i] | nbody_system.speed
            others[i].mass = 1 | nbody_system.mass
            others[i].radius = 0 | nbody_system.length
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(particles_in_binary)
        multiples_code.binaries.add_particle(binary)
        
        multiples_code.singles.add_particles(others)
        
        
        multiples_code.commit_particles()
        
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.components_of_multiples), 2)
        self.assertEqual(len(multiples_code.singles), 3)
        self.assertEqual(len(multiples_code.particles), 4)
        self.assertEqual(len(code.particles), 4)
        
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass,0.2 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[-1].mass,0.2 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[-1].position, [0,0,0] | nbody_system.length, 6)
        self.assertAlmostRelativeEquals(code.particles[-1].velocity, [0,0, 0] | nbody_system.speed, 6)
        
        
        multiples_code.update_model()
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass, 0.2 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[-1].mass, 0.2 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[-1].position, [0,0,0] | nbody_system.length, 6)
        self.assertAlmostRelativeEquals(code.particles[-1].velocity, [0,0, 0] | nbody_system.speed, 6)
        
        multiples_code.singles_in_binaries[0].mass = 0.2 | nbody_system.mass
        
        multiples_code.update_model()
        
        print(code.particles.mass)
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass, 0.3 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[-1].mass, 0.3 | nbody_system.mass)
        print(code.particles[-1].position)
        print(code.particles[-1].velocity)
        self.assertAlmostRelativeEquals(code.particles[-1].position, [0.00166666666667,0,0] | nbody_system.length, 6)
        self.assertAlmostRelativeEquals(code.particles[-1].velocity, [0, 0.7453559925, 0] | nbody_system.speed, 6)
        
        
        
      
    def test10(self):
        code = Hermite()
        
        particles_in_binary = self.new_binary(
            0.1 | nbody_system.mass,
            0.1 | nbody_system.mass,
            0.01 | nbody_system.length,
            keyoffset = 1
        )
        particles_in_binary.radius = 0.001 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        
        encounter_code.parameters.hard_binary_factor = 1
        encounter_code.small_scale_factor = 1
        
        others = datamodel.Particles(key = [4,5,6])
        for i in range(3):
            others[i].position = [i,  0, 0] | nbody_system.length
            others[i].velocity = [0, 0, i] | nbody_system.speed
            others[i].mass = 1 | nbody_system.mass
            others[i].radius  = 0.05 | nbody_system.length
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.must_handle_one_encounter_per_stopping_condition = False
        multiples_code.singles.add_particles(particles_in_binary)
        multiples_code.singles.add_particles(others)
        
        
        multiples_code.commit_particles()
        multiples_code.evolve_model(1 | nbody_system.time)
        
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.components_of_multiples), 2)
        self.assertEqual(len(multiples_code.singles), 3)
        self.assertEqual(len(multiples_code.particles), 4)
        self.assertEqual(len(code.particles), 4)
        
        self.assertEqual(id(multiples_code.singles_in_binaries), id(multiples_code.binaries[0].child1.particles_set))
        self.assertEqual(id(multiples_code.components_of_multiples), id(multiples_code.multiples[0].components[0].particles_set))
        #multiples_code.singles_in_binaries[0].mass = 0.2 | nbody_system.mass
        print(multiples_code.particles.mass)
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass, 1.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(multiples_code.particles.mass.sum(), 0.1 + 0.1 + 3.0 | nbody_system.mass)
        multiples_code.update_model()
        
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass, 1.1 | nbody_system.mass)
        
        index = -1
        if not code.particles[index].mass > 1.0| nbody_system.mass:
            index = -2
        self.assertAlmostRelativeEquals(code.particles[index].mass, 1.1 | nbody_system.mass)
        
        multiples_code.singles_in_binaries[0].mass += 0.2 | nbody_system.mass
        
        multiples_code.update_model()
        
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].mass, 1.3 | nbody_system.mass)
        self.assertAlmostRelativeEquals(code.particles[index].mass, 1.3 | nbody_system.mass)
    
    
      
    def test11(self):
        code = Hermite()
        
        particles_in_binary = self.new_binary(
            1.0 | nbody_system.mass,
            1.0 | nbody_system.mass,
            0.001 | nbody_system.length,
            keyoffset = 1
        )
        particles_in_binary.radius = 0.01 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        
        
        others = datamodel.Particles(keys = [4,5,6])
        for i in range(3):
            others[i].position = [i, 0, 0] | nbody_system.length
            others[i].velocity = [0, 0, 0] | nbody_system.speed
            others[i].mass = 0.2 | nbody_system.mass
            others[i].radius  = 0.05 | nbody_system.length
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles.add_particles(particles_in_binary)
        multiples_code.singles.add_particles(others)
        
        
        stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
        stopping_condition.enable()
        
        multiples_code.commit_particles()   
        multiples_code.evolve_model(1 | nbody_system.time)
        self.assertEqual(len(multiples_code.multiples), 1)        
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.components_of_multiples), 2)
        self.assertEqual(len(multiples_code.singles), 3)
        self.assertEqual(len(multiples_code.particles), 4)
        self.assertEqual(len(code.particles), 4)
        self.assertTrue(stopping_condition.is_set())
        multiples_code.particles[-1].velocity = [0, 0, 0] | nbody_system.speed
        multiples_code.update_model()
        print(multiples_code.particles.key)
        
        self.assertEqual(len(stopping_condition.particles(0)), 1)
        self.assertEqual(len(stopping_condition.particles(1)), 0)
        self.assertEqual(len(stopping_condition.particles(2)), 0)
        self.assertAlmostRelativeEquals(multiples_code.multiples[0].mass, 2.0 | nbody_system.mass)
        self.assertAlmostRelativeEquals(multiples_code.particles.mass.sum(), 2.6 | nbody_system.mass)
        print(multiples_code.particles.velocity)
        multiples_code.evolve_model(2 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        self.assertEqual(len(stopping_condition.particles(0)), 0)
        self.assertEqual(len(stopping_condition.particles(1)), 0)
        self.assertEqual(len(stopping_condition.particles(2)), 1)
        self.assertAlmostRelativeEquals(multiples_code.multiples[0].mass, 2.0 | nbody_system.mass)
        self.assertAlmostRelativeEquals(multiples_code.particles.mass.sum(), 2.6 | nbody_system.mass)
    
    
    def test12(self):
        code = Hermite()
        
        particles_in_binary = self.new_binary(
            1.0 | nbody_system.mass,
            1.0 | nbody_system.mass,
            0.001 | nbody_system.length,
            keyoffset = 10
        )
        particles_in_binary.radius = 0.01 | nbody_system.length
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        binary = datamodel.Particle(key=20)
        binary.child1 = particles_in_binary[0]
        binary.child2 = particles_in_binary[1]
        binary.position = [1,0,1] | nbody_system.length
        particles_in_binary.position += [1,0,1] | nbody_system.length
        
        others = datamodel.Particles(keys = [4,5,6])
        for i in range(3):
            others[i].position = [i*10, 0, 0] | nbody_system.length
            others[i].velocity = [0, 0, 0] | nbody_system.speed
            others[i].mass = 0.2 | nbody_system.mass
            others[i].radius  = 0.05 | nbody_system.length
            
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.particles.add_particles(others)
        multiples_code.singles_in_binaries.add_particles(particles_in_binary)
        multiples_code.binaries.add_particle(binary)
        multiples_code.commit_particles()   
        print(multiples_code.particles)
        self.assertEqual(len(multiples_code.particles), 4)
        self.assertAlmostRelativeEquals(multiples_code.particles[-1].position,  [1,0,1] | nbody_system.length)
        
        
        
    
    def test13(self):
        code = Hermite()
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        center_of_mass_particles = datamodel.Particles(5)
        center_of_mass_particles.position = (numpy.asarray(range(5))).reshape(5,1) * ([1.0, 0.0, 0.0] | nbody_system.length)
        center_of_mass_particles.velocity = [0.0, 0.0, 0.0] | nbody_system.speed
        center_of_mass_particles.radius  = 0.05 | nbody_system.length
        binaries, singles_in_binaries = self.create_binaries(
            center_of_mass_particles, 
            1 | nbody_system.mass,
            0.01 | nbody_system.mass,
            0.0001 | nbody_system.length
        )
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
        multiples_code.binaries.add_particles(binaries)
        multiples_code.commit_particles()   
        
        
        
        #stopping_condition = multiples_code.stopping_conditions.encounter_detection
        #stopping_condition.enable()
        stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
        stopping_condition.enable()
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)

        multiples_code.evolve_model(1 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(0):
            print("NEW:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(1):
            print("REMOVED:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(2):
            print("UPDATED:", x.key, x.child1.key, x.child2.key)
        for x in multiples_code.singles:
            print(x.key, x.mass)
        self.assertEqual(len(multiples_code.singles_in_binaries) + len(multiples_code.singles), 2*len(center_of_mass_particles))
        self.assertEqual(len(multiples_code.binaries) - len(stopping_condition.particles(0)) + len(stopping_condition.particles(1)),  len(center_of_mass_particles))
    
    def test14(self):
        code = Hermite()
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        center_of_mass_particles = datamodel.Particles(5)
        center_of_mass_particles.position = (numpy.asarray(range(5))).reshape(5,1) * ([1.0, 0.0, 0.0] | nbody_system.length)
        center_of_mass_particles.velocity = [0.0, 0.0, 0.0] | nbody_system.speed
        center_of_mass_particles.radius  = 0.05 | nbody_system.length
        binaries, singles_in_binaries = self.create_binaries(
            center_of_mass_particles, 
            1 | nbody_system.mass,
            0.1 | nbody_system.mass,
            0.00000001 | nbody_system.length
        )
            
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
        multiples_code.binaries.add_particles(binaries)
        multiples_code.commit_particles()   
        
        
        
        #stopping_condition = multiples_code.stopping_conditions.encounter_detection
        #stopping_condition.enable()
        stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
        stopping_condition.enable()
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)

        multiples_code.evolve_model(2 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(0):
            print("NEW:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(1):
            print("REMOVED:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(2):
            print("UPDATED:", x.key, x.child1.key, x.child2.key)
        for x in multiples_code.singles:
            print(x.key, x.mass)
        self.assertEqual(len(multiples_code.singles_in_binaries) + len(multiples_code.singles), 2*len(center_of_mass_particles))
        self.assertEqual(len(multiples_code.binaries) - len(stopping_condition.particles(0)) + len(stopping_condition.particles(1)),  len(center_of_mass_particles))
        
    def test15(self):
        code = Hermite()
        
        
        encounter_code = encounters.HandleEncounter(
            kepler_code =  self.new_kepler(),
            resolve_collision_code = self.new_smalln(),
        )
        n = 10
        center_of_mass_particles = plummer.new_plummer_model(n, random=numpy.random.mtrand.RandomState(1))
        center_of_mass_particles.radius  = 0.5 | nbody_system.length
        center_of_mass_particles.velocity *= 0
        binaries, singles_in_binaries = self.create_binaries(
            center_of_mass_particles, 
            0.999 * ((1.0 | nbody_system.mass) / n),
            0.001 * ((1.0 | nbody_system.mass) / n),
            0.00001 | nbody_system.length
        )
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounter_code
        )
        multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
        multiples_code.binaries.add_particles(binaries)
        multiples_code.commit_particles()   
        
        
        
        #stopping_condition = multiples_code.stopping_conditions.encounter_detection
        #stopping_condition.enable()
        stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
        stopping_condition.enable()
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)

        multiples_code.evolve_model(2 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        for x in multiples_code.binaries:
            print(x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(0):
            print("NEW:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(1):
            print("REMOVED:", x.key, x.child1.key, x.child2.key)
        for x in stopping_condition.particles(2):
            print("UPDATED:", x.key, x.child1.key, x.child2.key)
        for x in multiples_code.singles:
            print(x.key, x.mass)
        self.assertEqual(len(multiples_code.binaries) - len(stopping_condition.particles(0)) + len(stopping_condition.particles(1)),  len(center_of_mass_particles))
     
    
    def test16(self):
        code = Hermite()
        
        
        n = 10
        singles = datamodel.Particles(keys = list(range(1,n+1)))
        singles.mass = 1 | nbody_system.mass
        for x in range(n):
            singles[x].position = [x*x, 0, 0] | nbody_system.length
        singles.velocity = [0,0,0] | nbody_system.speed
            
        singles.radius  = 0.5 | nbody_system.length
        
        multiples_code = encounters.Multiples(
            gravity_code = code,
            handle_encounter_code = encounters.StickyHandleEncounter()
        )
        multiples_code.singles.add_particles(singles)
        multiples_code.commit_particles()   
        
        multiples_code.evolve_model(1 | nbody_system.time)
        print(len(multiples_code.multiples))
        self.assertEqual(len(multiples_code.multiples), 1)
        self.assertEqual(len(multiples_code.particles), 9)
        self.assertEqual(len(multiples_code.singles), 8)
        self.assertEqual(len(multiples_code.binaries), 1)
        self.assertEqual(len(multiples_code.singles_in_binaries), 2)
        self.assertEqual(id(multiples_code.components_of_multiples), id(multiples_code.multiples[0].components[0].particles_set))
        print(multiples_code.multiples[0].components)
        io.write_set_to_file(
            (
                multiples_code.singles, 
                multiples_code.singles_in_binaries, 
                multiples_code.binaries,  
                multiples_code.components_of_multiples, 
                multiples_code.multiples
            ), 
            "multiples.hdf5", 
            "hdf5",
            version = "2.0",
            names = (
                "singles",
                "singles_in_binaries",
                "binaries",
                "components_of_multiples",
                "multiples"
            )
            )
        
        multiples_code_loaded  = encounters.Multiples(
            gravity_code = Hermite(),
            handle_encounter_code = encounters.StickyHandleEncounter()
        )
        
        (
            singles,
            singles_in_binaries,
            binaries,
            components_of_multiples,
            multiples
        ) = io.read_set_from_file(
            "multiples.hdf5", 
            "hdf5",
            version = "2.0",
            names = (
                "singles",
                "singles_in_binaries",
                "binaries",
                "components_of_multiples",
                "multiples"
            )
        )
        self.assertEqual(len(multiples), 1)
        self.assertEqual(len(singles), 8)
        self.assertEqual(len(binaries), 1)
        self.assertEqual(len(singles_in_binaries), 2)
        #self.assertEquals(id(components_of_multiples), id(multiples[0].components[0].particles_set))
        
        multiples_code_loaded.singles.add_particles(singles)
        multiples_code_loaded.singles_in_binaries.add_particles(singles_in_binaries)
        multiples_code_loaded.binaries.add_particles(binaries)
        multiples_code_loaded.components_of_multiples.add_particles(components_of_multiples)
        multiples_code_loaded.multiples.add_particles(multiples)
        
        multiples_code_loaded.commit_particles()   
        
        self.assertEqual(len(multiples_code_loaded.multiples), 1)
        self.assertEqual(len(multiples_code_loaded.particles), 9)
        self.assertEqual(len(multiples_code_loaded.singles), 8)
        self.assertEqual(len(multiples_code_loaded.binaries), 1)
        self.assertEqual(len(multiples_code_loaded.singles_in_binaries), 2)
        #self.assertEquals(id(multiples_code_loaded.components_of_multiples), id(multiples_code_loaded.multiples[0].components[0].particles_set))
       
        multiples_code.evolve_model(4 | nbody_system.time)
        
        # need to use 3 here as the model_time is reset when doing a restart and we dit not set it after creating Hermite
        multiples_code_loaded.evolve_model(3.0 | nbody_system.time)
        

        print(len(multiples_code.multiples), multiples_code.particles)
        print(multiples_code.particles.position - multiples_code_loaded.particles.position)
        self.assertAlmostRelativeEquals(multiples_code.particles.position - multiples_code_loaded.particles.position, [0,0,0] | nbody_system.length)
        
        for code in [multiples_code, multiples_code_loaded]:
            self.assertEqual(len(code.multiples), 1)
            self.assertEqual(len(code.particles), 8)
            self.assertEqual(len(code.singles), 7)
            self.assertEqual(len(code.binaries), 1)
            self.assertEqual(len(code.singles_in_binaries), 2)
            self.assertEqual(len(code.components_of_multiples), 3)
            self.assertEqual(id(code.components_of_multiples), id(code.multiples[0].components[0].particles_set))
