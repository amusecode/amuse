from amuse.test.amusetest import TestWithMPI

import os
import sys

from amuse.community.hermite.interface import Hermite
from amuse.community.bhtree.interface import BHTree

import numpy
import threading
from amuse.units import nbody_system
from amuse.units import units
from amuse.datamodel import Particles
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestAmuseInterface(TestWithMPI):
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
            instance.particles.copy_values_of_all_attributes_to(particles)
            particles.savepoint()

    def test1(self):
        from amuse.rfi import channel
        channel.MpiChannel.ensure_mpi_initialized()
        is_threaded = channel.MpiChannel.is_threaded()
        is_multithreading_supported = channel.MpiChannel.is_multithreading_supported()
        self.assertEqual(is_threaded, is_multithreading_supported)

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        bhtree = BHTree(convert_nbody)
        bhtree.initialize_code()
        bhtree.eps2_for_gravity = 0.001
            
        bhtree_particles = self.new_system_sun_and_earth()
        bhtree.particles.add_particles(bhtree_particles)
        
        if bhtree.legacy_interface.channel_type == 'mpi':
            from mpi4py import MPI
            if not MPI.Query_thread() == MPI.THREAD_MULTIPLE:
                bhtree.stop()
                self.skip("can only test parallel with multiple thread support in mpi implementation")
            
        
        hermite = Hermite(convert_nbody)
        hermite.dt_dia = 5000
        hermite.commit_parameters()
            
        hermite_particles = self.new_system_sun_and_earth()
        hermite.particles.add_particles(hermite_particles)
        
        thread1 = threading.Thread(target = self.evolve_model_unit_day, args = (bhtree, bhtree_particles, 10))
        thread2 = threading.Thread(target = self.evolve_model_unit_day, args = (hermite, hermite_particles, 10))
        
        thread1.start()
        thread2.start()
        
        thread1.join()
        thread2.join()
        
        
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure()
            plot = figure.add_subplot(1,1,1)
            
            earth = bhtree_particles[1]
            x_points = earth.get_timeline_of_attribute("x")
            y_points = earth.get_timeline_of_attribute("y")
            x_points_in_AU = [t_x[1].value_in(units.AU) for t_x in x_points]
            y_points_in_AU = [t_x1[1].value_in(units.AU) for t_x1 in y_points]
            
            plot.scatter(x_points_in_AU,y_points_in_AU, color = "b", marker = 'o')
            
            earth = hermite_particles[1]
            x_points = earth.get_timeline_of_attribute("x")
            y_points = earth.get_timeline_of_attribute("y")
            x_points_in_AU = [t_x2[1].value_in(units.AU) for t_x2 in x_points]
            y_points_in_AU = [t_x3[1].value_in(units.AU) for t_x3 in y_points]
            
            plot.scatter(x_points_in_AU,y_points_in_AU, color = "g", marker = 'o')
            
            plot.set_xlim(-1.5, 1.5)
            plot.set_ylim(-1.5, 1.5)
            
            test_results_path = self.get_path_to_results()
            output_file = os.path.join(test_results_path, "parallel-earth-sun.svg")
            figure.savefig(output_file)
                    
        bhtree.stop()
        hermite.stop()
        bhtree.stop()

