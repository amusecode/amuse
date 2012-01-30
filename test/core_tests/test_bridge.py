import numpy

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import zero, AdaptingVectorQuantity, VectorQuantity
from amuse.support.exceptions import AmuseException
from amuse.datamodel import Particles
from amuse.datamodel import parameters
from amuse.ic.plummer import new_plummer_sphere

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
    
    def test2(self):
        print "CalculateFieldForParticles, nbody units, no gravity_constant exceptions"
        stars = new_plummer_sphere(100)
        self.assertRaises(AmuseException, bridge.CalculateFieldForParticles, stars, 
            expected_message = "For generic units the gravity_constant must be specified")
        self.assertRaises(AmuseException, bridge.CalculateFieldForParticles, Particles(), 
            expected_message = "Particle data not yet available, so the gravity_constant must be specified")
        instance = bridge.CalculateFieldForParticles(particles = stars, gravity_constant = nbody_system.G)
    
    def test3(self):
        print "CalculateFieldForParticles get_potential_at_point, no softening"
        epsilon = 0 | units.m
        
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        
        cluster = ExampleGravityCodeInterface()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        instance = bridge.CalculateFieldForParticles(particles = stars)
        instance.smoothing_length_squared = epsilon**2
        
        zeros = numpy.zeros(9) | units.parsec
        pos_range = numpy.linspace(-1.0, 1.0, 9) | units.parsec
        self.assertAlmostRelativeEqual(
            instance.get_potential_at_point(zeros, pos_range, zeros, zeros), 
            cluster.get_potential_at_point(zeros, pos_range, zeros, zeros))
        for a_calculate_field, a_code in zip(
                instance.get_gravity_at_point(zeros, pos_range, zeros, zeros), 
                cluster.get_gravity_at_point(zeros, pos_range, zeros, zeros)):
            self.assertAlmostRelativeEqual(a_calculate_field, a_code, 12)
    
    def test4(self):
        print "CalculateFieldForParticles get_potential_at_point, with softening"
        epsilon = 0.5 | units.parsec
        
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        
        cluster = ExampleGravityCodeInterface()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        instance = bridge.CalculateFieldForParticles(particles = stars)
        instance.smoothing_length_squared = epsilon**2
        
        zeros = numpy.zeros(9) | units.parsec
        pos_range = numpy.linspace(-1.0, 1.0, 9) | units.parsec
        self.assertAlmostRelativeEqual(
            instance.get_potential_at_point(zeros, pos_range, zeros, zeros), 
            cluster.get_potential_at_point(zeros, pos_range, zeros, zeros))
        for a_calculate_field, a_code in zip(
                instance.get_gravity_at_point(zeros, pos_range, zeros, zeros), 
                cluster.get_gravity_at_point(zeros, pos_range, zeros, zeros)):
            self.assertAlmostRelativeEqual(a_calculate_field, a_code, 12)
    
    def test5(self):
        print "CalculateFieldForParticles get_potential_at_point, with individual softening"
        epsilon = 0.5 | units.parsec
        
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        stars.radius = epsilon * numpy.random.uniform(low=0.4, high=3.0, size=len(stars))
        
        cluster = ExampleGravityCodeInterface(softening_mode="individual")
        cluster.particles.add_particles(stars)
        
        instance = bridge.CalculateFieldForParticles(particles=stars, softening_mode="individual")
        
        zeros = numpy.zeros(9) | units.parsec
        pos_range = numpy.linspace(-1.0, 1.0, 9) | units.parsec
        self.assertAlmostRelativeEqual(
            instance.get_potential_at_point(zeros, pos_range, zeros, zeros), 
            cluster.get_potential_at_point(zeros, pos_range, zeros, zeros))
        for a_calculate_field, a_code in zip(
                instance.get_gravity_at_point(zeros, pos_range, zeros, zeros), 
                cluster.get_gravity_at_point(zeros, pos_range, zeros, zeros)):
            self.assertAlmostRelativeEqual(a_calculate_field, a_code, 12)
    

class ExampleGravityCodeInterface(object):
    
    def __init__(self, softening_mode="shared"):
        self.particles = Particles()
        if softening_mode == "individual":
            self.softening_mode = "individual"
            self._softening_lengths_squared = self._softening_lengths_squared_individual
            self._softening_lengths = self._softening_lengths_individual
        else:
            self.softening_mode = "shared"
            self._softening_lengths_squared = self._softening_lengths_squared_shared
            self._softening_lengths = self._softening_lengths_shared
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
    
    def _softening_lengths_squared_individual(self):
        return self.particles.radius**2
    def _softening_lengths_squared_shared(self):
        return self.epsilon_squared.as_vector_with_length(len(self.particles))
    
    def _softening_lengths_individual(self):
        return self.particles.radius
    def _softening_lengths_shared(self):
        return self.epsilon_squared.sqrt().as_vector_with_length(len(self.particles))
    
    def initialize_code(self):
        self.model_time = 0 | units.Myr
    
    def get_potential_at_point(self, eps, x ,y, z):
        if isinstance(x, VectorQuantity):
            return -constants.G * (self.particles.mass.reshape((-1,1)) / 
                (self._softening_lengths_squared().reshape((-1,1)) + eps.reshape((1,-1))**2 + (self.particles.x.reshape((-1,1)) - 
                x.reshape((1,-1)))**2 + (self.particles.y.reshape((-1,1)) - y.reshape((1,-1)))**2 + 
                (self.particles.z.reshape((-1,1)) - z.reshape((1,-1)))**2).sqrt()).sum(axis=0)
        return -constants.G * (self.particles.mass / (self._softening_lengths_squared() + eps**2 + (self.particles.x - x)**2 + 
            (self.particles.y - y)**2 + (self.particles.z - z)**2).sqrt()).sum()
    
    def get_gravity_at_point(self, eps, x ,y, z):
        if isinstance(x, VectorQuantity):
            delta_x = x.reshape((1,-1)) - self.particles.x.reshape((-1,1))
            delta_y = y.reshape((1,-1)) - self.particles.y.reshape((-1,1))
            delta_z = z.reshape((1,-1)) - self.particles.z.reshape((-1,1))
            factor = -constants.G * (self.particles.mass.reshape((-1,1)) / 
                (self._softening_lengths_squared().reshape((-1,1)) + eps.reshape((1,-1))**2 + 
                delta_x**2 + delta_y**2 + delta_z**2)**1.5)
            return (factor*delta_x).sum(axis=0), (factor*delta_y).sum(axis=0), (factor*delta_z).sum(axis=0)
        delta_x = self.particles.x - x
        delta_y = self.particles.y - y
        delta_z = self.particles.z - z
        factor = -constants.G * (self.particles.mass / 
            (self._softening_lengths_squared() + eps**2 + delta_x**2 + delta_y**2 + delta_z**2)**1.5)
        return (factor*delta_x).sum(), (factor*delta_y).sum(), (factor*delta_z).sum()
    
    @property
    def potential_energy(self):
        if self.softening_mode == "individual":
            if len(self.particles) < 2:
                return zero * constants.G
            eps_vector = self.particles.radius**2
            sum_of_energies = zero
            for i in range(len(self.particles) - 1):
                dx = self.particles[i].x - self.particles[i+1:].x
                dy = self.particles[i].y - self.particles[i+1:].y
                dz = self.particles[i].z - self.particles[i+1:].z
                dr = ((dx * dx) + (dy * dy) + (dz * dz) + eps_vector[i] + eps_vector[i+1:]).sqrt()
                sum_of_energies -= (self.particles.mass[i] * self.particles.mass[i+1:] / dr).sum()
            return constants.G * sum_of_energies
        return self.particles.potential_energy(smoothing_length_squared=2*self.epsilon_squared)
    
    @property
    def kinetic_energy(self):
        return self.particles.kinetic_energy()
    
    def get_epsilon_squared(self):
        return self.epsilon_squared
    
    def set_epsilon_squared(self, epsilon_squared):
        self.epsilon_squared = epsilon_squared
    
    def commit_particles(self):
        self.set_accelerations()
        self.set_next_timestep()
    
    def set_accelerations(self):
        accelerations = self.get_gravity_at_point(self._softening_lengths(), 
            self.particles.x, self.particles.y, self.particles.z)
        self.particles.ax = accelerations[0]
        self.particles.ay = accelerations[1]
        self.particles.az = accelerations[2]
    
    def set_next_timestep(self):
        self.next_timestep = 0.01 * (self.particles.velocity / 
            self.particles.acceleration).lengths_squared().amin().sqrt()
    
    def evolve_model(self, t_end):
        while self.model_time < t_end:
            dt = min([self.next_timestep, t_end - self.model_time])
            self.particles.position += self.particles.velocity * dt + 0.5 * self.particles.acceleration * dt**2
            old_acceleration = self.particles.acceleration
            self.set_accelerations()
            self.particles.velocity += 0.5 * (old_acceleration + self.particles.acceleration) * dt
            self.model_time += dt
            self.set_next_timestep()

def system_from_particles(base_class, kwargs, particles, eps=None):
    interface = base_class(**kwargs)
    interface.initialize_code()
    if eps is not None:
        interface.parameters.epsilon_squared = eps**2 
    interface.particles.add_particles(particles)
    interface.commit_particles()
    return interface

class TestBridge(amusetest.TestWithMPI):
    
    def test1(self):
        print "Bridge potential energy with code's epsilon_squared as softening length"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 1.0e-2 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        
        cluster=test_class()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, dict(), first_half, epsilon)
        cluster2 = system_from_particles(test_class, dict(), second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,) )
        bridgesys.add_system(cluster2, (cluster1,) )
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    
    def test2(self):
        print "Bridge potential energy with radius as softening length"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 0.1 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        stars.radius = epsilon
        
        cluster=test_class()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, dict(softening_mode="individual"), first_half)
        cluster2 = system_from_particles(test_class, dict(), second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,), radius_is_eps=True)
        bridgesys.add_system(cluster2, (cluster1,), radius_is_eps=True)
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    
    def test3(self):
        print "Bridge potential energy with radius as softening length"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 0.1 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        stars.radius = epsilon
        
        cluster=test_class()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, dict(), first_half, epsilon)
        cluster2 = system_from_particles(test_class, dict(), second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,), radius_is_eps=True)
        bridgesys.add_system(cluster2, (cluster1,))
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    
    def test4(self):
        print "Bridge evolve_model"
        convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
        epsilon = 1.0e-2 | units.parsec
        test_class=ExampleGravityCodeInterface
        
        numpy.random.seed(12345)
        stars = new_plummer_sphere(100, convert_nbody=convert)
        
        cluster=test_class()
        cluster.initialize_code()
        cluster.parameters.epsilon_squared = epsilon**2 
        cluster.particles.add_particles(stars)
        cluster.commit_particles()
        
        first_half = stars.select_array(lambda x: (x > 0 | units.m), ['x'] )
        second_half = stars - first_half
        cluster1 = system_from_particles(test_class, dict(), first_half, epsilon)
        cluster2 = system_from_particles(test_class, dict(), second_half, epsilon)
        bridgesys=bridge.Bridge()
        bridgesys.add_system(cluster1, (cluster2,) )
        bridgesys.add_system(cluster2, (cluster1,) )
        
        self.assertAlmostRelativeEqual(
            cluster1.particles.get_intersecting_subset_in(cluster.particles).position,
            cluster1.particles.position)
        
#~        old = cluster1.particles.position
        for i in range(2):
            one_timestep = cluster.next_timestep
            cluster.evolve_model(cluster.model_time + one_timestep)
            bridgesys.evolve_model(bridgesys.model_time + one_timestep, timestep=one_timestep)
        
#~        print ((old - cluster1.particles.position)/old).lengths()
        self.assertAlmostRelativeEqual(
            cluster1.particles.get_intersecting_subset_in(cluster.particles).position,
            cluster1.particles.position)
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
    

