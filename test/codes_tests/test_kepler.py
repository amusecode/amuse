import os
import sys
import time
from amuse.test.amusetest import TestWithMPI

import numpy
import math
from amuse.units import nbody_system
from amuse.units import units, constants
from amuse import datamodel
from amuse.ic.plummer import new_plummer_model


from amuse.community.kepler_orbiters.interface import Kepler

class TestKeplerOrbiters(TestWithMPI):
    
    def test0(self):
        instance = Kepler()
        instance.stop()
    def new_system_of_sun_and_earth(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371)     
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))
        
        return stars

    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = Kepler(convert_nbody)#, redirection="none")#, debugger="xterm")
        instance.initialize_code()

        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        instance.central_particle.add_particle(stars[0])
        earth = instance.orbiters.add_particle(stars[1])
        instance.commit_particles()

        instance.evolve_model(365 | units.day)

        instance.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

















