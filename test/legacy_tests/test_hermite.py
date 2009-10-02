import unittest
import sys

from amuse.legacy.hermite0 import muse_dynamics_mpi as mpi_interface

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
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        hermite.add_particle(1, 11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = hermite.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
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
        
    def test4(self):
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        
        hermite.add_particle( [1,2,3,4]
            , [11.0,12.0,13.0,14.0]
            , [2.0,3.0,4.0,5.0]
            , [2.1,3.1,4.1,5.1]
            , [2.2,3.2,4.2,5.2]
            , [2.3,3.3,4.3,5.3]
            , [2.4,3.4,4.4,5.4]
            , [2.5,3.5,4.5,5.5]
            , [2.6,3.6,4.6,5.6])
        retrieved_state = hermite.get_state(1)
        print "result:", retrieved_state
        self.assertEquals(11.0,  retrieved_state['mass'])
        
        retrieved_state = hermite.get_state([2,3,4])
        self.assertEquals(12.0,  retrieved_state['mass'][0])
        self.assertEquals(hermite.get_number(), 4)
        hermite.cleanup_module()
        
    
    def test5(self):
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        n = 10000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        hermite.add_particle(ids
            , values
            , values
            , values
            , values
            , values
            , values
            , values
            , values)
        retrieved_state = hermite.get_state(3999)
        print "result:", retrieved_state
        self.assertEquals(3999.0,  retrieved_state['mass'])
        hermite.cleanup_module()
        
    def test6(self):
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        n = 4000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        for i in range(n-1):
            hermite.add_particle(ids[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i])
                
        retrieved_state = hermite.get_state(1)
        print "result:", retrieved_state
        self.assertEquals(1.0,  retrieved_state['mass'])
        hermite.cleanup_module()
        del hermite
        
class TestAmuseInterface(unittest.TestCase):
    def new_system_of_sun_and_earth(self):
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
        
        return stars
        
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = mpi_interface.Hermite(convert_nbody)
        hermite.parameters.epsilon_squared = 0.0 | units.AU**2
        hermite.setup_module()
        hermite.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        hermite.add_particles(stars)

        hermite.evolve_model(365.0 | units.day)
        hermite.update_particles(stars)
        
        postion_at_start = earth.position.get_value_at_time(0 | units.s)[1].value_in(units.AU)[0]
        postion_after_full_rotation = earth.position.value().value_in(units.AU)[0]
        
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 6)
        
        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        hermite.update_particles(stars)
        postion_after_half_a_rotation = earth.position.value().value_in(units.AU)[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 2)
        
        
        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        hermite.update_particles(stars)
        postion_after_half_a_rotation = earth.position.value().value_in(units.AU)[1]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        hermite.cleanup_module()
        del hermite

    
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = mpi_interface.Hermite(convert_nbody)
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.setup_module()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.add_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure(figsize = (40,40))
            plot = figure.add_subplot(1,1,1)
            
            
            for index, (time,position) in enumerate(earth.position.values):
                x_point = position.value_in(units.AU)[0]
                y_point = position.value_in(units.AU)[1]
                color = 'b'
                plot.plot([x_point],[y_point], color + 'o')
            
            figure.savefig("hermite-earth-sun.svg")    
        
        instance.cleanup_module()
        del instance
        
    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(units.MSun(1.0), units.km(149.5e6))

        instance = mpi_interface.Hermite(convert_nbody)
        instance.setup_module()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        instance.add_particles(stars)
        factor = instance.get_energies()[0] / (instance.get_energies()[1] / (-2 | units.none))
        print factor
        total_energy0 =  sum(instance.get_energies(), 0|units.J)
        self.assertAlmostEqual(factor.value_in(units.none), 1.000, 3)
        instance.evolve_model(100 | units.day)
        total_energy1 =  sum(instance.get_energies(), 0|units.J)
        factor = instance.get_energies()[0] / (instance.get_energies()[1] / (-2 | units.none))
        self.assertAlmostEqual(factor.value_in(units.none), 1.000, 3)
        instance.evolve_model(100 | units.day)
        total_energy2=  sum(instance.get_energies(), 0|units.J)
        self.assertAlmostEqual((total_energy2 /  total_energy0).value_in(units.none), 1.0, 7)
        self.assertAlmostEqual(total_energy2.value_in(units.J), total_energy1.value_in(units.J), 10)
        
