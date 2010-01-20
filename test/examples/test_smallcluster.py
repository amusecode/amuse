import sys
import unittest
import numpy 
import random
import collections
import os

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.sse.muse_stellar_mpi import SSE
from amuse.legacy.support.core import is_mpd_running

from amuse.support.io import store

from amuse.experiments.plummer import MakePlummerModel

class SalpeterIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125 | units.MSun, alpha = -2.35):
        self.mass_min = mass_min.as_quantity_in(units.MSun)
        self.mass_max = mass_max.as_quantity_in(units.MSun)
        self.alpha = alpha
        self.random = random.Random()
        self.random.seed()
    
    def mass_mean(self):
        alpha1 = self.alpha + 1
        alpha2 = self.alpha + 2
        l1 = pow(self.mass_min.value_in(units.MSun), alpha1)
        l2 = pow(self.mass_min.value_in(units.MSun), alpha2)
        u1 = pow(self.mass_max.value_in(units.MSun), alpha1)
        u2 = pow(self.mass_max.value_in(units.MSun), alpha2)
        return ((u2 - l2) * alpha1) / ((u1 - l1) * alpha2) | units.MSun
        
    def mass(self, random_number):
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max.value_in(units.MSun) / self.mass_min.value_in(units.MSun) , alpha1) - 1.0)
        return self.mass_min.value_in(units.MSun) * (pow(1 + (factor * random_number), 1.0 / alpha1)) | units.MSun
        
    def next_mass(self):
        return self.mass(self.random.random())
        
    def next_set(self, number_of_stars):
        set_of_masses = numpy.zeros(number_of_stars)
        total_mass = 0.0 | units.MSun
        for i in range(number_of_stars):
           mass = self.next_mass()
           set_of_masses[i] = mass.value_in(units.MSun)
           total_mass += mass
        return (total_mass, units.MSun.new_quantity(set_of_masses))
        
class SalpeterIMFTests(unittest.TestCase):
    def test1(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean().value_in(units.MSun), 0.351, 3)

    def test2(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass(1.0).value_in(units.MSun), 100, 3)
        self.assertAlmostEqual(instance.mass(0).value_in(units.MSun), 0.1, 3)
       
    def test3(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        n = 10000
        total_mass, set_of_masses = instance.next_set(10000)
        mean = total_mass.value_in(units.MSun) / float(n)
        exact_mean = instance.mass_mean().value_in(units.MSun)
        self.assertTrue(abs(mean - exact_mean) < 0.1)
        
    def test4(self):
        instance = SalpeterIMF(0.1 | units.MSun, 125 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual( 1.0 / instance.mass_mean().value_in(units.MSun), 2.8253, 4)
       

def move_particles_to_center_of_mass(particles):
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()
    
    center_of_mass = particles.center_of_mass()
    center_of_mass_velocity = particles.center_of_mass_velocity()
    
    particles.position = particles.position - center_of_mass
    particles.velocity = particles.velocity - center_of_mass_velocity     
    
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()
   
     
def plot_particles(particles, name_of_the_figure):
    
    if HAS_MATPLOTLIB:
        print "plotting the data"
        
        
        figure = pyplot.figure(figsize = (40, 40))
        plots = map(lambda x : figure.add_subplot(4,4, x+1), range(4*4))
        
        index = 0
        for data in particles.history:
            if index > 15:
                break
            
            x_values = data.x.value_in(units.parsec)
            y_values = data.y.value_in(units.parsec)
            mass_values = data.mass.value_in(units.MSun)
            
            sizes = mass_values * 10.0
            plots[index].scatter(x_values,y_values, marker='o', s=sizes)
            index += 1
            
        for plot in plots:
            plot.set_xlim(-2.0, 2.0)
            plot.set_ylim(-2.0, 2.0)

        figure.savefig(name_of_the_figure)
        if False:
            from matplotlib import axes3d
            figure = pyplot.figure()
            axes_3d = axes3d.Axes3D(figure)
            positions = particles.get_values_of_attribute('position')
            xs = numpy.array([position.x.value_in(units.lightyear) for position in positions])
            ys = numpy.array([position.y.value_in(units.lightyear) for position in positions])
            zs = numpy.array([position.z.value_in(units.lightyear) for position in positions])
            #print xs, yz, zs
            
            plot = axes_3d.scatter(xs, ys, zs)
            axes_3d.set_xlim(-10.0,10.0)
            axes_3d.set_ylim(-10.0,10.0)  
            axes_3d.set_zlim(-10.0,10.0)        
            
            figure.savefig("3d_"+name_of_the_figure)

def print_log(time, gravity, particles, total_energy_at_t0, total_energy_at_this_time):
    print "Evolved model to t = " + str(time)
    print total_energy_at_t0, total_energy_at_this_time, (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
    #print "KE:" , particles.kinetic_energy().as_quantity_in(units.J)
    #print "PE:" , particles.potential_energy(gravity.parameters.epsilon_squared)
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()
    
def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    initial_mass_function = SalpeterIMF()
    total_mass, salpeter_masses = initial_mass_function.next_set(number_of_stars)
    
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)
    convert_nbody.set_as_default()
    
    particles = MakePlummerModel(number_of_stars, convert_nbody).result;
   
    gravity = BHTree()
    gravity.setup_module()
    gravity.parameters.epsilon_squared = 0.15 | units.parsec ** 2
        
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_default_parameters() 
    
    print "setting masses of the stars"
    particles.radius = 0.0 | units.RSun
    particles.mass = salpeter_masses
    
    print "initializing the particles"
    stellar_evolution.particles.add_particles(particles)
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(particles)
    
    from_stellar_evolution_to_model.copy_attributes(["mass"])
    
    print "centering the particles"
    move_particles_to_center_of_mass(particles)
    
    gravity.particles.add_particles(particles)
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
    
        
    time = 0.0 | units.Myr    
    particles.savepoint(time)
    total_energy_at_t0 = sum(gravity.get_energies(), 0.0 | units.J)
    
    print "evolving the model until t = " + str(end_time)
    while time < end_time:
        time += 0.25 | units.Myr
        
        print "gravity evolve step starting"
        gravity.evolve_model(time)
        print "gravity evolve step done"
        
        print "stellar evolution step starting"
        stellar_evolution.evolve_model(time)
        print "stellar evolution step done"

        from_gravity_to_model.copy()
        from_stellar_evolution_to_model.copy_attributes(["mass", "radius"])

        particles.savepoint(time)  
        
        from_model_to_gravity.copy_attributes(["mass"])
        
        total_energy_at_this_time = sum(gravity.get_energies(), 0.0 | units.J)    
        print_log(time, gravity, particles, total_energy_at_t0, total_energy_at_this_time)

    
    if os.path.exists('small.hdf5'):
        os.remove('small.hdf5')
    storage = store.StoreHDF("small.hdf5")
    storage.store(particles)
   
    del gravity
    del stellar_evolution
    
    plot_particles(particles, name_of_the_figure)
    
        

        
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    assert is_mpd_running()
    simulate_small_cluster(4, 4 | units.Myr)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | units.Myr, sys.argv[3])
