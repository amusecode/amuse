# Run with
#
# nosetests-2.7 --nocapture --nologcapture -w test/codes_tests --tests=test_multiples 

from amuse.test.amusetest import TestWithMPI

import os
import sys
import numpy
import time
import math

from amuse.community.hermite0.interface import Hermite
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants

from amuse import datamodel
from amuse.ic import plummer
from amuse.couple import multiples

class TestSimpleMultiples(TestWithMPI):

    def new_smalln(self):
        result = SmallN()
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 2001
        return result
        
    def new_kepler_si(self):
        unit_converter = nbody_system.nbody_to_si(
            1 | units.MSun,
            1 | units.AU
        )
        kepler = Kepler(unit_converter)
        kepler.initialize_code()
        return kepler
        
    def new_kepler(self):
        kepler = Kepler()
        kepler.initialize_code()
        return kepler
        
    def new_smalln_si(self):
    
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        result = SmallN(converter)
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 2001
        return result
        
    def new_binary(self, mass1, mass2, semi_major_axis,
                   eccentricity = 0, keyoffset = 0):
        total_mass = mass1 + mass2
        mass_fraction_particle_1 = mass1 / (total_mass)
    
        binary = datamodel.Particles(keys=range(keyoffset, keyoffset+2))
        binary[0].mass = mass1
        binary[1].mass = mass2
    
        mu = nbody_system.G * total_mass
    
        velocity_perihelion = numpy.sqrt( mu / semi_major_axis  * ((1.0 + eccentricity)/(1.0 - eccentricity)))
        radius_perihelion = semi_major_axis * (1.0 - eccentricity)
        print velocity_perihelion
        
        binary[0].position = ((1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0,0.0,0.0])
        binary[1].position = -(mass_fraction_particle_1 * radius_perihelion * [1.0,0.0,0.0])
    
        binary[0].velocity = ((1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0,1.0,0.0])
        binary[1].velocity = -(mass_fraction_particle_1 * velocity_perihelion * [0.0,1.0,0.0])

        return binary
        
    def test0(self):
        code = Hermite()
        stars = datamodel.Particles(2)
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0,0,0],
            [1.2, 0, 0]
        ]|nbody_system.length
        stars.velocity = [
            [0,0,0],
            [0,0.1, 0]
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        code.particles.add_particles(stars)
        
        multiples_code = multiples.Multiples(code, self.new_smalln, self.new_kepler())
        print multiples_code.multiples_energy_correction
        total_energy0 = multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy0
        multiples_code.evolve_model(0.6|nbody_system.time)
        total_energy1 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy1
        print 
        print total_energy0
        error = abs((total_energy1 - total_energy0)/total_energy0)
        print multiples_code.multiples_energy_correction
        
        self.assertTrue(error < 1e-7)
        self.assertAlmostRelativeEquals(multiples_code.multiples_energy_correction - multiples_code.kinetic_energy, -total_energy0, 7)
    
    def test1(self):
        code = Hermite()
        stars = datamodel.Particles(3)
        stars.mass = 1 | nbody_system.mass
        stars.position = [
            [0,0,0],
            [1.2, 0, 0],
            [-10, 0, 0],
        ]|nbody_system.length
        stars.velocity = [
            [0,0,0],
            [0,0.1, 0],
            [0,0, 0],
        ]|nbody_system.speed
        stars.radius = 0.5 | nbody_system.length
        code.particles.add_particles(stars)
        
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        print converter.to_si(stars.velocity)
        print converter.to_si(0.6|nbody_system.time).as_quantity_in(units.Myr)
        
        multiples_code = multiples.Multiples(code, self.new_smalln, self.new_kepler())
        print multiples_code.multiples_energy_correction
        total_energy0 = multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy0
        multiples_code.evolve_model(0.6|nbody_system.time)
        total_energy1 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy1
        print multiples_code.multiples_energy_correction
        error = abs((total_energy1 - total_energy0)/total_energy0)
        print "error:", error
        self.assertTrue(error < 1e-7)
        
    def test2(self):
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        
        code = Hermite(converter)
        stars = datamodel.Particles(2)
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
        code.particles.add_particles(stars)
        
        multiples_code = multiples.Multiples(code, self.new_smalln_si, self.new_kepler_si(), gravity_constant = constants.G)
        print multiples_code.multiples_energy_correction
        total_energy0 = multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy0
        multiples_code.evolve_model(converter.to_si(0.6|nbody_system.time))
        total_energy1 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        print total_energy1
        print total_energy0
        
        print converter.to_nbody(total_energy1)
        print converter.to_nbody(total_energy0)
        error = abs((total_energy1 - total_energy0)/total_energy0)
        print multiples_code.multiples_energy_correction
        print converter.to_nbody(multiples_code.multiples_energy_correction)
        self.assertTrue(error < 1e-7)
        self.assertAlmostRelativeEquals(multiples_code.multiples_energy_correction - multiples_code.kinetic_energy, -total_energy0, 7)

    def test3(self):
        code = Hermite()
        stars = datamodel.Particles()
        binary1 = self.new_binary(
            1.0 | nbody_system.mass,
            0.4 | nbody_system.mass,
            0.01 |  nbody_system.length,
            keyoffset = 1
        )
        binary2 = self.new_binary(
            0.5 | nbody_system.mass,
            0.2 | nbody_system.mass,
            0.005 |  nbody_system.length,
            keyoffset = 3
        )
        binary2.position += [0.5,0,0] | nbody_system.length
        binary2.velocity += [0,0,0.03] | nbody_system.speed
        stars.add_particles(binary1)
        stars.add_particles(binary2)
        stars.radius = 0.005 | nbody_system.length
        print stars
        code.particles.add_particles(stars)

        multiples_code = multiples.Multiples(code, self.new_smalln, self.new_kepler())
        total_energy0 = multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        multiples_code.evolve_model(0.1|nbody_system.time)
        multiples_code.print_multiples()
        total_energy1 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction

        error = abs((total_energy1 - total_energy0)/total_energy0)
        
        self.assertTrue(error < 1e-7)
        multiples_code.evolve_model(0.2|nbody_system.time)
        multiples_code.print_multiples()
        total_energy2 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction

        error = abs((total_energy2 - total_energy0)/total_energy0)
        
        self.assertTrue(error < 1e-5)
        
        stars = multiples_code.stars
        self.assertEquals(len(stars), 4)
        self.assertEquals(len(code.particles), 2)
    

    def test4(self):
        code = Hermite()
        stars = datamodel.Particles()
        binary1 = self.new_binary(
            1.0 | nbody_system.mass,
            1.0 | nbody_system.mass,
            0.1 | nbody_system.length,
            keyoffset = 1
        )
        binary2 = self.new_binary(
            1.0 | nbody_system.mass,
            1.0 | nbody_system.mass,
            0.2 | nbody_system.length,
            keyoffset = 3
        )
        binary1.position += [1.0,0.0,0] | nbody_system.length
        binary2.position -= [1.0,0.0,0] | nbody_system.length
        stars.add_particles(binary1)
        stars.add_particles(binary2)
        stars.radius = 0.25 | nbody_system.length
        
        code.particles.add_particles(stars)

        multiples_code = multiples.Multiples(code, self.new_smalln, self.new_kepler())
        total_energy0 = (
            multiples_code.kinetic_energy 
            + multiples_code.potential_energy 
            - multiples_code.multiples_energy_correction
        )
        multiples_code.evolve_model(0.1|nbody_system.time)
        multiples_code.print_multiples()
        total_energy1 =  (
            multiples_code.kinetic_energy
            + multiples_code.potential_energy
            - multiples_code.multiples_energy_correction
        )

        error = abs((total_energy1 - total_energy0)/total_energy0)
        
        self.assertTrue(error < 1e-7)
        multiples_code.evolve_model(2.0|nbody_system.time)
        multiples_code.print_multiples()
        total_energy2 =(
            multiples_code.kinetic_energy
            + multiples_code.potential_energy 
            - multiples_code.multiples_energy_correction
        )

        error = abs((total_energy2 - total_energy0)/total_energy0)
        
        self.assertTrue(error < 1e-5)
        
    def test5(self):
        code1 = Hermite()
        stars = datamodel.Particles()
        binary1 = self.new_binary(
            1.0 | nbody_system.mass,
            0.2 | nbody_system.mass,
            0.01 |  nbody_system.length,
            eccentricity = 0.7,
            keyoffset = 1
        )
        binary2 = self.new_binary(
            0.8 | nbody_system.mass,
            0.3 | nbody_system.mass,
            0.01 |  nbody_system.length,
            eccentricity = 0.9,
            keyoffset = 3
        )
        binary1.position += [-0.5,-0.5,-0.5] | nbody_system.length
        binary1.velocity += [0.0,0.0,0.0] | nbody_system.speed
        binary2.position += [0.5,0.5,0.5] | nbody_system.length
        binary2.velocity += [0.0,0.0,0.0] | nbody_system.speed
        stars.add_particles(binary1)
        stars.add_particles(binary2)
        stars.radius = 0.2 | nbody_system.length
        code1.particles.add_particles(stars)
        kepler1 = self.new_kepler()
        print binary2.velocity, (binary2[0].velocity - binary2[1].velocity).length
        multiples_code1 = multiples.Multiples(code1, self.new_smalln, kepler1)
        multiples_code1.evolve_model(0.1|nbody_system.time)
        multiples_code1.print_multiples()
        stars = multiples_code1.stars
        self.assertEquals(len(stars), 4)
        code2 = Hermite()
        code2.particles.add_particles(stars)
        kepler2 = self.new_kepler()
        multiples_code2 = multiples.Multiples(code2, self.new_smalln, kepler2)
        multiples_code1.evolve_model(0.2|nbody_system.time)
        multiples_code2.evolve_model(0.1|nbody_system.time)
        multiples_code1.print_multiples()
        multiples_code2.print_multiples()
        # compare the root positions
        p01 = multiples_code1._inmemory_particles[0]
        p02 = p01.as_particle_in_set(multiples_code2._inmemory_particles)
        self.assertAlmostRelativeEquals(p01.position, p01.position)
        p11 = multiples_code1._inmemory_particles[1]
        p12 = p11.as_particle_in_set(multiples_code2._inmemory_particles)
        self.assertAlmostRelativeEquals(p11.position, p11.position)
        
        if p01.mass == 1.2 | nbody_system.mass:
            b1 = p01
            b2 = p11
        else:
            b1 = p11
            b2 = p01
        b1stars = binary1.get_intersecting_subset_in(multiples_code1.stars)
        self.assertAlmostRelativeEquals(
            b1stars.center_of_mass(),
            b1.position,
        )
        b2stars = binary2.get_intersecting_subset_in(multiples_code1.stars)
        self.assertAlmostRelativeEquals(
            b2stars.center_of_mass(),
            b2.position,
        )
        
        
        
        
