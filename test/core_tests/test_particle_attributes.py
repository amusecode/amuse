import random
import numpy.random

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
        print "Test basic particle attributes and scale_to_standard - nbody units"
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
        print "Test basic particle attributes and scale_to_standard - SI units"
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
    
    def test3(self):
        print "Test new_particle_from_cluster_core - nbody units"
        numpy.random.seed(123)
        plummer = new_plummer_sphere(10000)
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, reuse_hop=True)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | nbody_system.length
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostEqual(result.density, 3.55015420914 | nbody_system.density)
        
        plummer.vx += 42 | nbody_system.speed
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|nbody_system.speed)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | nbody_system.length
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, reuse_hop=False)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e3 | nbody_system.density, 4)
    
    def test4(self):
        print "Test new_particle_from_cluster_core - SI units"
        numpy.random.seed(123)
        converter = nbody_system.nbody_to_si(1000.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(10000, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, reuse_hop=True)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | units.parsec
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostEqual(result.density, 3.55015420914e3 | units.MSun * units.parsec**-3)
        
        plummer.vx += 42 | units.km / units.s
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|units.km / units.s)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | units.parsec
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, reuse_hop=False)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e6 | units.MSun * units.parsec**-3, 4)
    
    def test5(self):
        print "Test new_particle_from_cluster_core - reuse_hop or not"
        converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(100, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=True)
        
        # Hop wasn't stopped, will use same Hop instance:
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=True)
        
        nbody_plummer = new_plummer_sphere(100)
        # Hop wasn't stopped, unit_converters don't match:
        self.assertRaises(AttributeError, nbody_plummer.new_particle_from_cluster_core, 
            expected_message="Cannot combine units from different systems: m and length")
        
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=False)
        
        # Hop was stopped, new instance will be made with supplied unit_converter (None in this case):
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=True)
        
        self.assertRaises(ConvertArgumentsException, plummer.new_particle_from_cluster_core, unit_converter=converter,#,
            expected_message="error while converting parameter 'mass', error: Cannot express kg in mass, the units do not have the same bases")
        
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=False)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=False)
    
    def test6(self):
        print "Test all particle attributes using Hop - each different function creates its own instance of Hop"
        numpy.random.seed(123)
        nbody_plummer = new_plummer_sphere(100)
        nbody_plummer.mass = new_salpeter_mass_distribution_nbody(100)
        
        # Each different function creates its own instance of Hop
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=True)
        result = nbody_plummer.bound_subset(G=nbody_system.G, reuse_hop=True)
        result = nbody_plummer.mass_segregation_Gini_coefficient(reuse_hop=True)
        result = nbody_plummer.LagrangianRadii(reuse_hop=True)
        result = nbody_plummer.densitycentre_coreradius_coredens(reuse_hop=True)
        
        converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0 | units.parsec)
        si_plummer = ParticlesWithUnitsConverted(nbody_plummer, converter.as_converter_from_si_to_nbody())
        functions_using_hop = [
            si_plummer.new_particle_from_cluster_core, 
            si_plummer.bound_subset, 
            si_plummer.mass_segregation_Gini_coefficient, 
            si_plummer.LagrangianRadii, 
            si_plummer.densitycentre_coreradius_coredens
        ]
        
        # Each fails since the Hop instance it tries to reuse has a different unit_converter
        for function_using_hop in functions_using_hop:
            self.assertRaises(ConvertArgumentsException, function_using_hop, unit_converter=converter,
                expected_message="error while converting parameter 'mass', error: Cannot express kg in mass, the units do not have the same bases")
        
        # Close all Hop instances:
        nbody_results = []
        nbody_results.append(nbody_plummer.new_particle_from_cluster_core(reuse_hop=False))
        nbody_results.append(nbody_plummer.bound_subset(G=nbody_system.G, reuse_hop=False))
        nbody_results.append(nbody_plummer.mass_segregation_Gini_coefficient(reuse_hop=False))
        nbody_results.append(nbody_plummer.LagrangianRadii(reuse_hop=False))
        nbody_results.append(nbody_plummer.densitycentre_coreradius_coredens(reuse_hop=False))
        
        # Now it works, because the Hop instances were closed, and new ones will be instantiated
        si_results = []
        for function_using_hop in functions_using_hop:
            si_results.append(function_using_hop(unit_converter=converter))
        
        convert = converter.as_converter_from_si_to_nbody()
        self.assertAlmostRelativeEqual(si_results[0].position, 
            ParticlesWithUnitsConverted(nbody_results[0].as_set(), convert)[0].position)
        self.assertAlmostRelativeEqual(si_results[1].position, 
            ParticlesWithUnitsConverted(nbody_results[1], convert).position)
        self.assertAlmostRelativeEqual(si_results[2], nbody_results[2], places=10)
        self.assertAlmostRelativeEqual(si_results[3][0], convert.from_target_to_source(nbody_results[3][0]), places=10)
        for in_si, in_nbody in zip(si_results[4], nbody_results[4]):
            self.assertAlmostRelativeEqual(in_si, convert.from_target_to_source(in_nbody), places=10)
    
    def test7(self):
        print "Test minimum_spanning_tree_length"
        particles = Particles(6)
        particles.position = [[-5,0,0], [-1,0,0], [0,0,0], [0,1,0], [0,-2,0], [-1,0.1,0]] | units.m
        self.assertEqual(particles[:1].minimum_spanning_tree_length(), 0 | units.m)
        self.assertEqual(particles[:2].minimum_spanning_tree_length(), 4 | units.m)
        self.assertEqual(particles[:3].minimum_spanning_tree_length(), 5 | units.m)
        self.assertEqual(particles[:4].minimum_spanning_tree_length(), 6 | units.m)
        self.assertEqual(particles[:5].minimum_spanning_tree_length(), 8 | units.m)
        self.assertEqual(particles[:6].minimum_spanning_tree_length(), 8.1 | units.m)
    
    def test8(self):
        print "Test mass_segregation_ratio"
        numpy.random.seed(123)
        random.seed(456)
        number_of_particles = 10000
        particles = new_plummer_sphere(number_of_particles)
        particles.r_squared = particles.position.lengths_squared()
        sorted = particles.sorted_by_attribute("r_squared")
        
        sorted.mass = numpy.random.uniform(1.0, 2.0, number_of_particles) | nbody_system.mass
        MSR = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=10)
        self.assertAlmostEquals(MSR, 0.8877, 3)
        random.seed(456)
        result = sorted.mass_segregation_ratio(number_of_particles=10, number_of_random_sets=10, also_compute_uncertainty=True)
        self.assertTrue(isinstance(result, particle_attributes.MassSegregationRatioResults))
        self.assertAlmostEquals(result.mass_segregation_ratio, 0.8877, 3)
        self.assertAlmostEquals(result.uncertainty, 0.2482, 3)
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
        print "Test __doc__ and help for particle attributes"
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
        
        self.assertEquals(particles.total_mass(), 2 | units.kg)
        self.assertEquals(particles.total_momentum(), [7, 0, 0] | units.kg * units.m / units.s)
        self.assertEquals(particles.total_momentum(), particles.total_mass() * particles.center_of_mass_velocity())
        self.assertEquals(particles.total_radius(), 0.5 | units.m)
        
        convert_nbody = nbody_system.nbody_to_si(1000 | units.kg, 1e-6 | units.m)
        numpy.random.seed(123)
        field = new_plummer_sphere(10000, convert_nbody) # small clump of particles, can be regarded as point mass
        self.assertAlmostRelativeEquals(particles.potential_energy_in_field(field), -constants.G * (1500 | units.kg**2 / units.m), 5)
        self.assertAlmostEquals(particles.potential_energy_in_field(field), -1.001142 | 1e-7 * units.kg * units.m**2 / units.s**2, 5)
        
        field.position *= ((5 | units.m) / field.position.lengths()).reshape((-1, 1)) # spherical shell around particles
        potential_energy = particles.potential_energy_in_field(field)
        particles.position += [0, 1, 2] | units.m # as long as particles remain inside the shell, the potential doesn't change
        self.assertAlmostEquals(particles.potential_energy_in_field(field), potential_energy, 5)
        
        particles.mass = [1, 2] | units.kg
        self.assertAlmostRelativeEquals(particles.potential(), -constants.G * ([2, 1] | units.kg / units.m))
        self.assertAlmostRelativeEquals(particles.potential()[0], particles[0].potential())
        self.assertAlmostRelativeEquals(particles.potential()[1], particles[1].potential())
    
    def test11(self):
        print "Test nearest_neighbour"
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
        print "Test correlation_dimension"
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
        print "Test box_counting_dimension"
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
        self.assertEquals(
            sorted(particles.a.get_attribute_names_defined_in_store()), 
            ['foo']
        )
        self.assertEquals(
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
