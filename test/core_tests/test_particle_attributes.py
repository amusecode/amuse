import random
import numpy.random
import sys

from amuse.test import amusetest
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support.interface import ConvertArgumentsException

from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.datamodel import Particle, Particles, ParticlesWithUnitsConverted
from amuse.datamodel import particle_attributes

class TestParticlesAttributes(amusetest.TestCase):
    
    def test1(self):
        print("Test basic particle attributes and scale_to_standard - nbody units")
        particles = Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | nbody_system.length
        particles.velocity = [[-1, 0, 0], [1,0,0]] | nbody_system.length/nbody_system.time
        particles.mass = 0.4 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.total_mass(), 0.8 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.4 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.08 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | nbody_system.length)
        particles.scale_to_standard()
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.25 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.5 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 1.0 | nbody_system.length)
        
        particles.scale_to_standard(virial_ratio=1) # unbound
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.5 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.5 | nbody_system.energy)
        particles.scale_to_standard(virial_ratio=0) # velocities zeroed
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.5 | nbody_system.energy)
    
    def test2(self):
        print("Test basic particle attributes and scale_to_standard - SI units")
        convert_nbody = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        particles = Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | units.parsec
        particles.velocity = [[-1, 0, 0], [1,0,0]] | units.parsec / units.Myr
        particles.mass = 0.5 | units.MSun 
        
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0 | units.MSun)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 1.0 * (0.5 | units.MSun) * (1 |units.parsec / units.Myr) **2 )
        self.assertAlmostRelativeEquals(particles.potential_energy(), -constants.G *  (0.5 | units.MSun) ** 2  / ([2,0,0] | units.parsec).length() )
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | units.parsec)
        
        particles.scale_to_standard(convert_nbody)
        self.assertAlmostRelativeEquals(particles.total_mass(), convert_nbody.to_si(1.0 | nbody_system.mass))
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), convert_nbody.to_si(0.25 | nbody_system.energy))
        self.assertAlmostRelativeEquals(particles.potential_energy().as_quantity_in(units.J), convert_nbody.to_si(-0.5 | nbody_system.energy).as_quantity_in(units.J), 12)
        self.assertAlmostRelativeEquals(particles.virial_radius(), convert_nbody.to_si(1.0 | nbody_system.length))
        
        particles.scale_to_standard(convert_nbody, virial_ratio=1) # unbound
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.5 * constants.G * (1 | units.MSun**2 / units.parsec), 13)
        self.assertAlmostRelativeEquals(particles.potential_energy(), -0.5 * constants.G * (1 | units.MSun**2 / units.parsec))
        particles.scale_to_standard(convert_nbody, virial_ratio=0) # velocities zeroed
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0 | units.J)
        self.assertAlmostRelativeEquals(particles.potential_energy(), -0.5 * constants.G * (1 | units.MSun**2 / units.parsec))
    
    def test7(self):
        print("Test minimum_spanning_tree_length")
        particles = Particles(6)
        particles.position = [[-5,0,0], [-1,0,0], [0,0,0], [0,1,0], [0,-2,0], [-1,0.1,0]] | units.m
        self.assertEqual(particles[:1].minimum_spanning_tree_length(), 0 | units.m)
        self.assertEqual(particles[:2].minimum_spanning_tree_length(), 4 | units.m)
        self.assertEqual(particles[:3].minimum_spanning_tree_length(), 5 | units.m)
        self.assertEqual(particles[:4].minimum_spanning_tree_length(), 6 | units.m)
        self.assertEqual(particles[:5].minimum_spanning_tree_length(), 8 | units.m)
        self.assertEqual(particles[:6].minimum_spanning_tree_length(), 8.1 | units.m)
    
    def test8(self):
        print("Test mass_segregation_ratio")
        numpy.random.seed(123)
        random.seed(456)
        number_of_particles = 10000
        particles = new_plummer_sphere(number_of_particles)
        particles.r_squared = particles.position.lengths_squared()
        sorted = particles.sorted_by_attribute("r_squared")
        
        sorted.mass = numpy.random.uniform(1.0, 2.0, number_of_particles) | nbody_system.mass
        MSR = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=10)
        
        if sys.hexversion > 0x03000000:
            self.assertAlmostEqual(MSR, 0.7160, 3)
        else:
            self.assertAlmostEqual(MSR, 0.8877, 3)
                
        random.seed(456)
        result = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=10, also_compute_uncertainty=True)
        self.assertTrue(isinstance(result, particle_attributes.MassSegregationRatioResults))
        if sys.hexversion > 0x03000000:
            self.assertAlmostEqual(result.mass_segregation_ratio, 0.7160, 3)
            self.assertAlmostEqual(result.uncertainty, 0.2321, 3)
        else:
            self.assertAlmostEqual(result.mass_segregation_ratio, 0.8877, 3)
            self.assertAlmostEqual(result.uncertainty, 0.2482, 3)
        MSR, sigma = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=50, also_compute_uncertainty=True)
        self.assertTrue(MSR - sigma < 1.0 < MSR + sigma)
        
        sorted.mass = numpy.linspace(1.0, 2.0, number_of_particles) | nbody_system.mass
        MSR, sigma = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=20, 
            also_compute_uncertainty=True)
        self.assertTrue(MSR < 0.1)
        self.assertTrue(sigma < MSR)
        
        sorted.mass = numpy.linspace(2.0, 1.0, number_of_particles) | nbody_system.mass
        MSR, sigma = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=20, 
            also_compute_uncertainty=True)
        self.assertTrue(MSR > 10.0)
        self.assertTrue(sigma < MSR)
    
    def test9(self):
        print("Test __doc__ and help for particle attributes")
        particles = Particles(2)
        self.assertTrue("Returns the total kinetic energy of the\n    particles in the particles set." in particles.kinetic_energy.__doc__)
        self.assertEqual(particles.kinetic_energy.__class__.__name__, "BoundParticlesFunctionAttribute")
        self.assertEqual(particles.kinetic_energy.__class__._function.__name__, "kinetic_energy")
        # __doc__ must be defined on the __class__ for the Python help to work:
        self.assertTrue("Returns the total kinetic energy of the\n    particles in the particles set." in particles.kinetic_energy.__class__.__doc__)
        
        # difference between particle function and particleS function:
        self.assertTrue("Returns the specific kinetic energy of each particle in the set." in particles.specific_kinetic_energy.__doc__)
        self.assertTrue("Returns the specific kinetic energy of each particle in the set." in particles.specific_kinetic_energy.__class__.__doc__)
        self.assertTrue("Returns the specific kinetic energy of the particle." in particles[0].specific_kinetic_energy.__doc__)
        self.assertTrue("Returns the specific kinetic energy of the particle." in particles[0].specific_kinetic_energy.__class__.__doc__)
        
        self.assertTrue("Returns the potential at the position of each particle in the set." in particles.potential.__doc__)
        self.assertTrue("Returns the potential at the position of each particle in the set." in particles.potential.__class__.__doc__)
        self.assertTrue("Returns the potential at the position of the particle." in particles[0].potential.__doc__)
        self.assertTrue("Returns the potential at the position of the particle." in particles[0].potential.__class__.__doc__)
    
    def test10(self):
        particles = Particles(2)
        particles.position = [[1, 0, 0], [2,0,0]] | units.m
        particles.velocity = [[3, 0, 0], [4,0,0]] | units.m / units.s
        particles.mass = 1 | units.kg
        
        self.assertEqual(particles.total_mass(), 2 | units.kg)
        self.assertEqual(particles.total_momentum(), [7, 0, 0] | units.kg * units.m / units.s)
        self.assertEqual(particles.total_momentum(), particles.total_mass() * particles.center_of_mass_velocity())
        self.assertEqual(particles.total_radius(), 0.5 | units.m)
        
        convert_nbody = nbody_system.nbody_to_si(1000 | units.kg, 1e-6 | units.m)
        numpy.random.seed(123)
        field = new_plummer_sphere(10000, convert_nbody) # small clump of particles, can be regarded as point mass
        self.assertAlmostRelativeEquals(particles.potential_energy_in_field(field), -constants.G * (1500 | units.kg**2 / units.m), 5)
        self.assertAlmostEqual(particles.potential_energy_in_field(field), -1.001142 | 1e-7 * units.kg * units.m**2 / units.s**2, 5)
        
        field.position *= ((5 | units.m) / field.position.lengths()).reshape((-1, 1)) # spherical shell around particles
        potential_energy = particles.potential_energy_in_field(field)
        particles.position += [0, 1, 2] | units.m # as long as particles remain inside the shell, the potential doesn't change
        self.assertAlmostEqual(particles.potential_energy_in_field(field), potential_energy, 5)
        
        particles.mass = [1, 2] | units.kg
        self.assertAlmostRelativeEquals(particles.potential(), -constants.G * ([2, 1] | units.kg / units.m))
        self.assertAlmostRelativeEquals(particles.potential()[0], particles[0].potential())
        self.assertAlmostRelativeEquals(particles.potential()[1], particles[1].potential())
    
    def test11(self):
        print("Test nearest_neighbour")
        particles = Particles(21)
        particles.x = numpy.logspace(0.0, 2.0, 21) | units.m
        particles.y = 0.0 | units.m
        particles.z = 0.0 | units.m
        self.assertEqual(particles.nearest_neighbour()[0], particles[1])
        self.assertEqual(particles.nearest_neighbour()[1:].key, particles[:-1].key)
        
        neighbours = Particles(3)
        neighbours.x = [1.0, 10.0, 100.0] | units.m
        neighbours.y = 0.0 | units.m
        neighbours.z = 0.0 | units.m
        self.assertEqual(particles.nearest_neighbour(neighbours).key, neighbours.key[[0]*8 + [1]*10 + [2]*3])
        
        # A few tests to check the correct behaviour of 'max_array_length' (to prevent memory overflow)
        nearest_neighbours = particles.nearest_neighbour(max_array_length=3*21*21) # all in one go
        self.assertEqual(nearest_neighbours[0], particles[1])
        self.assertEqual(nearest_neighbours[1:].key, particles[:-1].key)
        nearest_neighbours = particles.nearest_neighbour(max_array_length=3*21*21-1) # two passes
        self.assertEqual(nearest_neighbours[0], particles[1])
        self.assertEqual(nearest_neighbours[1:].key, particles[:-1].key)
        nearest_neighbours = particles.nearest_neighbour(max_array_length=1) # 21 passes, one for each particle
        self.assertEqual(nearest_neighbours[0], particles[1])
        self.assertEqual(nearest_neighbours[1:].key, particles[:-1].key)
        self.assertEqual(particles.nearest_neighbour(neighbours, max_array_length=189).key, neighbours.key[[0]*8 + [1]*10 + [2]*3]) # all in one go
        self.assertEqual(particles.nearest_neighbour(neighbours, max_array_length=188).key, neighbours.key[[0]*8 + [1]*10 + [2]*3]) # two passes
        self.assertEqual(particles.nearest_neighbour(neighbours, max_array_length=1).key, neighbours.key[[0]*8 + [1]*10 + [2]*3]) # 21 passes, one for each particle
    
    def new_koch_star(self, level=5):
        height = numpy.sqrt(3) / 6.0
        def next_iteration(x_values, y_values):
            dx = x_values[1:] - x_values[:-1]
            dy = y_values[1:] - y_values[:-1]
            x_one_third = x_values[:-1] + dx / 3.0
            y_one_third = y_values[:-1] + dy / 3.0
            x_two_third = x_one_third + dx / 3.0
            y_two_third = y_one_third + dy / 3.0
            x_new_point = x_values[:-1] + 0.5*dx - height * dy
            y_new_point = y_values[:-1] + 0.5*dy + height * dx
            new_x = numpy.append(numpy.dstack((x_values[:-1], x_one_third, x_new_point, x_two_third)), [x_values[-1]])
            new_y = numpy.append(numpy.dstack((y_values[:-1], y_one_third, y_new_point, y_two_third)), [y_values[-1]])
            return new_x, new_y
        x, y = numpy.array([0.0, 0.5, 1.0, 0.0]), numpy.array([0.0, 3*height, 0.0, 0.0])
        for i in range(level):
            x, y = next_iteration(x, y)
        return x, y
    
    def test12(self):
        print("Test correlation_dimension")
        # Particles distributed uniformly in 3D
        particles = Particles(729)
        particles.position = numpy.mgrid[0:9.0, 0:9.0, 0:9.0].reshape(3, -1).transpose() | units.m
        dimension = particles.correlation_dimension()
        self.assertAlmostRelativeEquals(dimension, 3.0, 1)
        # Fractal dimension is scale-free
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.correlation_dimension(), 10)
        
        # Particles distributed in the x-y plane
        particles.position = numpy.concatenate((numpy.mgrid[0:27.0, 0:27.0].reshape(2, -1), numpy.zeros((1, 729)))).transpose() | units.m
        dimension = particles.correlation_dimension()
        self.assertAlmostRelativeEquals(dimension, 2.0, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.correlation_dimension(), 10)
        
        # Particles distributed along a line
        particles.position = numpy.concatenate([numpy.arange(729.0).reshape(1, -1)]*3).transpose() | units.m
        dimension = particles.correlation_dimension()
        self.assertAlmostRelativeEquals(dimension, 1.0, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.correlation_dimension(), 10)
        
        # Particles on a Koch curve
        x, y = self.new_koch_star(level=6)
        numpy.random.seed(123456)
        sel = numpy.random.randint(len(x), size=729)
        particles.position = numpy.hstack((x[sel], y[sel], numpy.zeros(729))) | units.m
        dimension = particles.correlation_dimension()
        self.assertAlmostRelativeEquals(dimension, 1.26186, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.correlation_dimension(), 10)
    
    def test13(self):
        print("Test box_counting_dimension")
        # Particles distributed uniformly in 3D
        particles = Particles(4096)
        particles.position = numpy.mgrid[0:16.0, 0:16.0, 0:16.0].reshape(3, -1).transpose() | units.m
        dimension = particles.box_counting_dimension()
        self.assertAlmostRelativeEquals(dimension, 3.0, 1)
        # Fractal dimension is scale-free
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.box_counting_dimension(), 10)
        
        # Particles distributed in the x-y plane
        particles.position = numpy.concatenate((numpy.mgrid[0:64.0, 0:64.0].reshape(2, -1), numpy.zeros((1, 4096)))).transpose() | units.m
        dimension = particles.box_counting_dimension()
        self.assertAlmostRelativeEquals(dimension, 2.0, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.box_counting_dimension(), 10)
        
        # Particles distributed along a line
        particles.position = numpy.concatenate([numpy.arange(4096.0).reshape(1, -1)]*3).transpose() | units.m
        dimension = particles.box_counting_dimension()
        self.assertAlmostRelativeEquals(dimension, 1.0, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.box_counting_dimension(), 10)
        
        # Particles on a Koch curve
        x, y = self.new_koch_star(level=7)
        numpy.random.seed(123456)
        sel = numpy.random.randint(len(x), size=4096)
        particles.position = numpy.hstack((x[sel], y[sel], numpy.zeros(4096))) | units.m
        dimension = particles.box_counting_dimension()
        self.assertAlmostRelativeEquals(dimension, 1.26186, 1)
        particles.position *= 1000
        self.assertAlmostRelativeEquals(dimension, particles.box_counting_dimension(), 10)
    
    def test14(self):
        print("Test mass_segregation_from_nearest_neighbour")
        numpy.random.seed(123)
        random.seed(4567)
        number_of_particles = 1000
        particles = new_plummer_sphere(number_of_particles)
        particles.r_squared = particles.position.lengths_squared()
        sorted = particles.sorted_by_attribute("r_squared")
        
        sorted.mass = numpy.random.uniform(1.0, 2.0, number_of_particles) | nbody_system.mass
        MSR, sigma = sorted.mass_segregation_from_nearest_neighbour(number_of_particles=10, also_compute_uncertainty=True)
        
        if sys.hexversion > 0x03000000:
            self.assertAlmostEqual(MSR, 1.72632, 3)
            self.assertAlmostEqual(sigma, 0.4127, 3)
        else:
            self.assertAlmostEqual(MSR, 1.7355, 3)
            self.assertAlmostEqual(sigma, 0.3969, 3)
                
        random.seed(456)
        MSR_of_nonsegregated_systems = []
        for i in range(10):
            sorted.mass = numpy.random.uniform(1.0, 2.0, number_of_particles) | nbody_system.mass
            MSR = sorted.mass_segregation_from_nearest_neighbour(number_of_particles=10, number_of_random_sets=50)
            MSR_of_nonsegregated_systems.append(MSR)
        self.assertAlmostEqual((MSR_of_nonsegregated_systems|units.none).mean(), 1.0, 1)
        self.assertAlmostEqual((MSR_of_nonsegregated_systems|units.none).std(), 0.3, 1)
        
        sorted.mass = numpy.linspace(2.0, 1.0, number_of_particles) | nbody_system.mass
        MSR, sigma = sorted.mass_segregation_from_nearest_neighbour(number_of_particles=10, number_of_random_sets=20, 
            also_compute_uncertainty=True)
        self.assertTrue(MSR > 5.0)
        
        if sys.hexversion > 0x03000000:
            self.assertAlmostEqual(sigma, 0.3, 1)
        else:
            self.assertAlmostEqual(sigma, 0.4, 1)

    def test15(self):
        scale_R = 1.0 | units.parsec
        scale_M = 1000.0 | units.MSun
        converter = nbody_system.nbody_to_si(scale_M,scale_R)

        for n in range(3162, 3165, 1):
            stars = new_plummer_sphere(
                    n, 
                    convert_nbody=converter,
                    )
            stars.mass=(numpy.arange(1,n+1)/(1.*n)) | units.MSun
            potential = stars.potential()

            for i,x in enumerate(stars):
              self.assertAlmostRelativeEqual(potential[i],x.potential())

    def test16(self):
        scale_R = 1.0 | units.parsec
        scale_M = 1000.0 | units.MSun
        converter = nbody_system.nbody_to_si(scale_M,scale_R)

        for n in range(0, 50):
            stars = new_plummer_sphere(
                    50, 
                    convert_nbody=converter,
                    )
            stars.mass=numpy.arange(1,51) | units.MSun
            potential = stars.potential(block_size=n)

            for i,x in enumerate(stars):
              self.assertAlmostRelativeEqual(potential[i],x.potential())


class TestParticlesDomainAttributes(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        
        self.assertAlmostRelativeEqual(particles.a.foo,  1 | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  2 | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  3 | units.s)
        
    def test2(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.bar = 2 | units.kg
        particles.foo = 3 | units.s
        self.assertEqual(
            sorted(particles.a.get_attribute_names_defined_in_store()), 
            ['foo']
        )
        self.assertEqual(
            sorted(particles.b.get_attribute_names_defined_in_store()), 
            ['bar']
        )
        
    def test3(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        particles.a.add_particle(Particle(foo = 2 | units.m))
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [1,1,2] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2,0] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3,0] | units.s)
        
        particles.add_particle(Particle(foo = 2 | units.s))
        self.assertAlmostRelativeEqual(particles.a.foo,  [1,1,2,0] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2,0,0] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3,0,2] | units.s)
        
        
    def test4(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        def set_a():
            particles.a = 1 | units.kg
        self.assertRaises(AttributeError, set_a)

    def test5(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        
        particles.a[0].foo = 3 | units.m
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [3,1] | units.m)
        
        particles[0].a.foo = 4 | units.m
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [4,1] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3] | units.s)

    def test6(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        def set_a():
            particles[0].a = 1 | units.kg
        self.assertRaises(AttributeError, set_a)
