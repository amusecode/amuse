import unittest
import sys

from amuse.legacy.hermite0 import muse_dynamics_mpi as mpi_interface

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

import numpy

from matplotlib import pyplot


class TestMPIInterface(unittest.TestCase):
    
    def test1(self):
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        new_state = hermite.dynamics_state()
        new_state.id = 1
        new_state.mass = 11.0
        hermite.add_particle(new_state)
        retrieved_state = hermite.get_state(1)
        self.assertEquals(new_state.mass,  retrieved_state.mass)
        self.assertEquals(hermite.get_number(), 1)
        hermite.cleanup_module()
        del hermite
        
    def test2(self):
        hermite = mpi_interface.Hermite()
        hermite.eps2 = 0.101
        self.assertEquals(0.101, hermite.eps2)
        hermite.eps2 = 0.110
        self.assertEquals(0.110, hermite.eps2)
        del hermite
        
    def test3(self):
        hermite = mpi_interface.Hermite()
        hermite.flag_collision = 1
        self.assertEquals(1, hermite.flag_collision)
        hermite.flag_collision = 0
        self.assertEquals(0, hermite.flag_collision)
        del hermite
        
class TestAmuseInterface(unittest.TestCase):
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(units.MSun(1.0), units.km(149.5e6))

        hermite = mpi_interface.Hermite(convert_nbody)
        hermite.setup_module()
        hermite.dt_dia = 5000
        
        sun = core.Particle(0)
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = core.Particle(1)
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))

        hermite.add_star(sun)
        hermite.add_star(earth)

        hermite.evolve_model(365.0 | units.day)
        hermite.update_star(earth)
        
        postion_at_start = earth.position.get_value_at_time(0 | units.s)[1].in_(units.AU).number[0]
        postion_after_full_rotation = earth.position.value().in_(units.AU) .number[0]
        
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 6)
        
        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        hermite.update_star(earth)
        postion_after_half_a_rotation = earth.position.value().in_(units.AU) .number[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 2)
        
        
        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        hermite.update_star(earth)
        postion_after_half_a_rotation = earth.position.value().in_(units.AU) .number[1]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        hermite.cleanup_module()
    
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(units.MSun(1.0), units.km(149.5e6))

        instance = mpi_interface.Hermite(convert_nbody)
        instance.setup_module()
        instance.dt_dia = 5000
        
        
        sun = core.Particle(0)
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = core.Particle(1)
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))

        instance.add_star(sun)
        instance.add_star(earth)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.update_star(earth)
        
        figure = pyplot.figure(figsize = (40,40))
        plot = figure.add_subplot(1,1,1)
        
        
        for index, (time,position) in enumerate(earth.position.values):
            x_point = position.in_(units.AU).number[0]
            y_point = position.in_(units.AU).number[1]
            color = 'b'
            plot.plot([x_point],[y_point], color + 'o')
        
        figure.savefig("hermite-earth-sun.svg")    
        
        instance.cleanup_module()
        del instance

