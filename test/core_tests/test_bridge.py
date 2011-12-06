import numpy

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import zero, AdaptingVectorQuantity, VectorQuantity
from amuse.support.exceptions import AmuseException
from amuse.datamodel import Particles
from amuse.datamodel import parameters
from amuse.ic.plummer import new_plummer_model

from amuse.test import amusetest
from amuse.couple import bridge

class TestCalculateFieldForParticles(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        
        instance = bridge.CalculateFieldForParticles(particles = particles, gravity_constant = nbody_system.G)
        
        zero = 0.0 | nbody_system.length
        print instance.get_gravity_at_point([zero], [1.0] | nbody_system.length, [zero], [zero])
        fx, fy, fz = instance.get_gravity_at_point([zero], [1.0] | nbody_system.length, [zero], [zero])
        self.assertAlmostEqual(fx, [0.0] | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fy, [0.0] | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fz, [0.0] | nbody_system.acceleration, 6)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point([zero], [x0], [zero], [zero])
            potential1 = instance.get_potential_at_point([zero], [x1], [zero], [zero])
            fx0, fy0, fz0 = instance.get_gravity_at_point([zero], [x0], [zero], [zero])
            fx1, fy1, fz1 = instance.get_gravity_at_point([zero], [x1], [zero], [zero])
            
            self.assertAlmostEqual(fy0[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz0[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fy1[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz1[0], 0.0 | nbody_system.acceleration, 6)
            
            self.assertAlmostEqual(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0[0], 5)
            self.assertAlmostEqual(potential0, potential1, 6)

class ExampleGravityCodeInterface(object):
    
    def __init__(self):
        self.particles = Particles()
        epsilon_squared_parameter = parameters.ModuleMethodParameterDefinition(
            "get_epsilon_squared",
            "set_epsilon_squared",
            "epsilon_squared",
            "gravitational softening length squared",
            default_value = 0.0 | nbody_system.length**2,
            must_set_before_get = False
        )
        self.parameters = parameters.new_parameters_instance_with_docs([epsilon_squared_parameter], self)
        self.epsilon_squared = 0.0 | nbody_system.length**2
    
    def initialize_code(self):
        pass
    
    def get_potential_at_point(self, eps, x ,y, z):
        if isinstance(x, VectorQuantity):
            return -constants.G * (self.particles.mass.reshape((-1,1)) / 
                (self.epsilon_squared + eps.reshape((1,-1))**2 + (self.particles.x.reshape((-1,1)) - 
                x.reshape((1,-1)))**2 + (self.particles.y.reshape((-1,1)) - y.reshape((1,-1)))**2 + 
                (self.particles.z.reshape((-1,1)) - z.reshape((1,-1)))**2).sqrt()).sum(axis=0)
        return -constants.G * (self.particles.mass / (self.epsilon_squared + eps**2 + (self.particles.x - x)**2 + 
            (self.particles.y - y)**2 + (self.particles.z - z)**2).sqrt()).sum()
    
    @property
    def potential_energy(self):
        return self.particles.potential_energy(smoothing_length_squared=2*self.epsilon_squared)
    
    @property
    def kinetic_energy(self):
        return self.particles.kinetic_energy()
    
    def get_epsilon_squared(self):
        return self.epsilon_squared
    
    def set_epsilon_squared(self, epsilon_squared):
        self.epsilon_squared = epsilon_squared
    

def system_from_particles(base_class, particles, eps=None):
    interface = base_class()
    interface.initialize_code()
    if eps is not None:
        interface.parameters.epsilon_squared = eps**2 
    interface.particles.add_particles(particles)
    return interface

class TestBridge(amusetest.TestWithMPI):
    
    def test1(self):
        print "Bridge potential energy with code's epsilon_squared as softening length"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 1.0e-2 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_model(100, convert_nbody=convert)
        
        cluster=test_class()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, first_half, epsilon)
        cluster2 = system_from_particles(test_class, second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,) )
        bridgesys.add_system(cluster2, (cluster1,) )
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    
    def test2(self):
        print "Bridge potential energy with radius as softening length"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 1.0e-2 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_model(100, convert_nbody=convert)
        stars.radius = epsilon
        
        cluster=test_class()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, first_half, epsilon)
        cluster2 = system_from_particles(test_class, second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,), radius_is_eps=True)
        bridgesys.add_system(cluster2, (cluster1,), radius_is_eps=True)
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    

