import unittest
import sys

from amuse.legacy.bhtree import muse_dynamics_mpi as mpi_interface

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

import numpy


try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
    

class TestMPIInterface(unittest.TestCase):
    
    def test1(self):
        instance = mpi_interface.BHTree()
        instance.setup_module()
        instance.add_particle(1, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(instance.get_number(), 1)
        instance.cleanup_module()
        del instance
        
    def test2(self):
        instance = mpi_interface.BHTree()
        instance.eps2 = 0.101
        self.assertEquals(0.101, instance.eps2)
        instance.eps2 = 0.110
        self.assertEquals(0.110, instance.eps2)
        del instance
        
    def test3(self):
        instance = mpi_interface.BHTree()
        instance.flag_collision = 1
        self.assertEquals(1, instance.flag_collision)
        instance.flag_collision = 0
        self.assertEquals(0, instance.flag_collision)
        del instance
        
class TestAmuseInterface(unittest.TestCase):
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = mpi_interface.BHTree(convert_nbody)
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        instance.setup_module()
        
        stars = core.Stars(2)
        
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

        instance.add_particles(stars)

        instance.evolve_model(365.0 | units.day)
        instance.update_particles(stars)
        
        postion_at_start = earth.position.get_value_at_time(0 | units.s)[1].in_(units.AU).number[0]
        postion_after_full_rotation = earth.position.value().in_(units.AU) .number[0]
       
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.value().in_(units.AU) .number[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 2)
        
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
         
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.value().in_(units.AU) .number[1]
        
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 1)
        instance.cleanup_module()
        del instance
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(units.MSun(1.0), units.km(149.5e6))

        instance = mpi_interface.BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.timestep = 0.0001
        #instance.use_self_gravity = 0
        instance.setup_module()
        
        stars = core.Stars(2)
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

        instance.add_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure(figsize = (40,40))
            plot = figure.add_subplot(1,1,1)
            
            for index, (time,position) in enumerate(earth.position.values):
                x_point = position.in_(units.AU).number[0]
                y_point = position.in_(units.AU).number[1]
                color = 'b'
                plot.plot([x_point],[y_point], color + 'o')
            plot.set_xlim(-1.5, 1.5)
            plot.set_ylim(-1.5, 1.5)
            
            figure.savefig("bhtree-earth-sun.svg")    
        
        instance.cleanup_module()
        del instance

