import unittest
import sys

from amuse.legacy.hermite0.muse_dynamics_mpi import Hermite
from amuse.legacy.bhtree.muse_dynamics_mpi import BHTree

from amuse.support.data.core import Particles
from amuse.support.units import nbody_system
from amuse.support.units import units

import numpy
import threading

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestAmuseInterface(unittest.TestCase):
    def new_system_sun_and_earth(self):
        result = Particles(2)
        sun = result[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = result[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))
        return result
        
    def evolve_model_unit_day(self, instance, particles, day):
        delta_days = 5
        for x in range(1, day + delta_days, delta_days):
            instance.evolve_model(x | units.day)
            instance.update_particles(particles)
            
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        bhtree = BHTree(convert_nbody)
        bhtree.eps2_for_gravity = 0.001
        bhtree.setup_module()
            
        bhtree_particles = self.new_system_sun_and_earth()
        bhtree.add_particles(bhtree_particles)
        
        hermite = Hermite(convert_nbody)
        hermite.dt_dia = 5000
        hermite.setup_module()
            
        hermite_particles = self.new_system_sun_and_earth()
        hermite.add_particles(hermite_particles)
        
        thread1 = threading.Thread(target = self.evolve_model_unit_day, args = (bhtree, bhtree_particles, 365))
        thread2 = threading.Thread(target = self.evolve_model_unit_day, args = (hermite, hermite_particles, 365))
        
        thread1.start()
        thread2.start()
        
        thread1.join()
        thread2.join()
        
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure(figsize = (10,10))
            plot = figure.add_subplot(1,1,1)
            
            
            earth = list(bhtree_particles)[1]
            for index, (time,position) in enumerate(earth.position.values):
                x_point = position.value_in(units.AU)[0]
                y_point = position.value_in(units.AU)[1]
                color = 'b'
                plot.plot([x_point],[y_point], color + 'o')
                
            earth = list(hermite_particles)[1]
            for index, (time,position) in enumerate(earth.position.values):
                x_point = position.value_in(units.AU)[0]
                y_point = position.value_in(units.AU)[1]
                color = 'g'
                plot.plot([x_point],[y_point], color + 'o')
                
            plot.set_xlim(-1.1, 1.1)
            plot.set_ylim(-1.1, 1.1)
            
            figure.savefig("earth-sun.svg")    
        
        bhtree.cleanup_module()
        hermite.cleanup_module()
        del bhtree
        del hermite

