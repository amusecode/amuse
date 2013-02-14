from amuse.test import amusetest

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
import numpy
import time
import sys

from amuse.datamodel import Particles
from amuse.datamodel import Particle
from amuse.couple import encounters

class TestAbstractHandleEncounter(amusetest.TestWithMPI):
    def new_binary(self, mass1, mass2, semi_major_axis,
                   eccentricity = 0, keyoffset = 1):
        total_mass = mass1 + mass2
        mass_fraction_particle_1 = mass1 / (total_mass)
    
        binary = Particles(keys=range(keyoffset, keyoffset+2))
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
        
    def test1(self):
        particles_in_encounter = Particles(2)
        particles_in_encounter.mass = 1 | nbody_system.mass
        particles_in_encounter[0].position = [1,0,0] | nbody_system.length
        particles_in_encounter[1].position = [0,0,0] | nbody_system.length
        particles_in_encounter.velocity = [0,0.5,0] | nbody_system.speed
        
        particles_in_field = Particles(5)
        particles_in_field.mass = [1,2,3,4,5] | nbody_system.mass
        particles_in_field[0].position = [2,0,0] | nbody_system.length
        particles_in_field[1].position = [1.5,0,0] | nbody_system.length
        particles_in_field[2].position = [0.5,1,0] | nbody_system.length
        particles_in_field[3].position = [0.5,-0.5,0] | nbody_system.length
        particles_in_field[4].position = [0,0,2] | nbody_system.length
        particles_in_field.velocity = [0,0.0,0] | nbody_system.speed
        
        x = encounters.AbstractHandleEncounter(
            particles_in_encounter,
            particles_in_field,
            G = nbody_system.G
        )
        
        x.start()
        
        self.assertAlmostRelativeEqual(x.large_scale_of_particles_in_the_encounter, 1.0 | nbody_system.length)
        self.assertAlmostRelativeEqual(x.small_scale_of_particles_in_the_encounter, 1.0 | nbody_system.length)
        
        self.assertEquals(len(x.particles_close_to_encounter), 3)
        self.assertAlmostRelativeEqual(x.particles_close_to_encounter.mass, [4,2,3] | nbody_system.mass)
        
    def test2(self):
        
        particles_in_encounter = Particles(keys=(1,2))
        particles_in_encounter.mass = 1 | nbody_system.mass
        particles_in_encounter[0].position = [1,0,0] | nbody_system.length
        particles_in_encounter[1].position = [0,0,0] | nbody_system.length
        particles_in_encounter.velocity = [0,0.5,0] | nbody_system.speed
        
        particles_in_field = Particles(keys=(3,))
        particles_in_field.mass = 1 | nbody_system.mass
        particles_in_field[0].position = [3,0,0] | nbody_system.length
        
        particles_in_multiples = Particles()
        particles_in_multiples.add_particle(particles_in_encounter[0])
        multiple = particles_in_multiples[0]
        multiple.components = Particles(2)
        multiple.components.mass = 0.5 | nbody_system.mass
        multiple.components[0].position = [0,0.1,0] | nbody_system.length
        multiple.components[1].position = [0,-0.1,0] | nbody_system.length
        multiple.components[0].velocity = [0,0,0.2] | nbody_system.speed
        multiple.components[1].velocity = [0,0.1,-0.2] | nbody_system.speed
        multiple.components.child1 = None
        multiple.components.child2 = None
        
        x = encounters.AbstractHandleEncounter(
            particles_in_encounter,
            particles_in_field,
            existing_multiples = particles_in_multiples,
            G = nbody_system.G
        )
        
        x.start()
        
        self.assertEquals(len(x.all_singles_in_encounter), 3)
        
        self.assertAlmostRelativeEqual(x.all_singles_in_encounter.mass, [0.5, 0.5, 1.0]| nbody_system.mass)
        
        child1 = x.singles_and_multiples_after_evolve[0]
        child2 = x.singles_and_multiples_after_evolve[1]
        self.assertAlmostRelativeEqual(child1.position, [1,0.1,0]| nbody_system.length)
        self.assertAlmostRelativeEqual(child1.velocity, [0,0.5,0.2]| nbody_system.speed)
        self.assertAlmostRelativeEqual(child2.position, [1,-0.1,0]| nbody_system.length)
        self.assertAlmostRelativeEqual(child2.velocity, [0,0.6,-0.2]| nbody_system.speed)
        
        
    def test3(self):
        
        particles_in_encounter = Particles(keys=(1,2))
        particles_in_encounter.mass = 1 | nbody_system.mass
        particles_in_encounter[0].position = [1,0,0] | nbody_system.length
        particles_in_encounter[1].position = [0,0,0] | nbody_system.length
        particles_in_encounter.velocity = [0,0.0,0] | nbody_system.speed
        
        particles_in_field = Particles(keys=(2,4,5,6,7,))
        particles_in_field.mass = [1,2,3,4,5] | nbody_system.mass
        particles_in_field[0].position = [2,0,0] | nbody_system.length
        particles_in_field[1].position = [1.5,0,0] | nbody_system.length
        particles_in_field[2].position = [0.5,1,0] | nbody_system.length
        particles_in_field[3].position = [0.5,-0.5,0] | nbody_system.length
        particles_in_field[4].position = [0,0,2] | nbody_system.length
        
        x = encounters.AbstractHandleEncounter(
            particles_in_encounter,
            particles_in_field,
            G = nbody_system.G
        )
        
        simple_binary = self.new_binary(
            1 | nbody_system.mass, 
            1 | nbody_system.mass, 
            0.2 | nbody_system.length
        )
                
        def evolve_singles_in_encounter_until_end_state():
            particles = x.singles_and_multiples_after_evolve
            particles.add_particles(x.all_singles_in_encounter)
            particles.child1 = None
            particles.child2 = None
            
            root_particle = particles.add_particle(Particle(
                key = 10,
                mass = 2.0 | nbody_system.mass, 
                position = particles[0:1].center_of_mass(),
                velocity = particles[0:1].center_of_mass_velocity(),
            ))
            root_particle.child1 = particles[0]
            root_particle.child2 = particles[1]
            particles[0].position = simple_binary[0].position + root_particle.position
            particles[1].position = simple_binary[1].position + root_particle.position
            
            particles[0].velocity = simple_binary[0].velocity + root_particle.velocity
            particles[1].velocity = simple_binary[1].velocity + root_particle.velocity
        
        x.evolve_singles_in_encounter_until_end_state = evolve_singles_in_encounter_until_end_state
        x.determine_structure_of_the_evolved_state = lambda : 1
        
        x.start()
        self.assertEquals(len(x.new_multiples), 1)
        multiple = x.new_multiples[0]
        self.assertEquals(len(multiple.components), 2)
        
        self.assertAlmostRelativeEqual(multiple.components[0].position, simple_binary[0].position)
        self.assertAlmostRelativeEqual(multiple.components[1].position, simple_binary[1].position) 
        
