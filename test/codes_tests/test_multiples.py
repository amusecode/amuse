from amuse.test.amusetest import TestWithMPI

import os
import sys
import numpy
import time
import math

from amuse.community.hermite0.interface import Hermite
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
        
    def new_smalln_si(self):
    
        converter = nbody_system.nbody_to_si(units.MSun, units.parsec)
        result = SmallN(converter)
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 2001
        return result
        
    def new_binary(self, mass1, mass2, semi_major_axis, eccentricity = 0, keyoffset = 0):
        total_mass = mass1 + mass2
        mass_fraction_particle_1 = mass1 / (total_mass)
    
        binary = datamodel.Particles(keys=range(keyoffset, keyoffset+2))
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
        
        multiples_code = multiples.Multiples(code, self.new_smalln)
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
        
        multiples_code = multiples.Multiples(code, self.new_smalln)
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
        
        multiples_code = multiples.Multiples(code, self.new_smalln_si, gravity_constant = constants.G)
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
        
        
        multiples_code = multiples.Multiples(code, self.new_smalln)
        total_energy0 = multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction
        multiples_code.evolve_model(0.1|nbody_system.time)
        multiples_code.print_multiples()
        total_energy1 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction

        error = abs((total_energy1 - total_energy0)/total_energy0)
        
        self.assertTrue(error < 1e-7)
        multiples_code.evolve_model(0.6|nbody_system.time)
        multiples_code.print_multiples()
        total_energy2 =  multiples_code.kinetic_energy + multiples_code.potential_energy - multiples_code.multiples_energy_correction

        error = abs((total_energy2 - total_energy0)/total_energy0)
        
        print code.particles
        print error
        self.assertTrue(error < 1e-5)
        #self.assertTrue(False)
        
