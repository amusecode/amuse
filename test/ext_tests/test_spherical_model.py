import numpy
from amuse.test.amusetest import TestCase
from amuse.support.exceptions import AmuseWarning, AmuseException
from amuse.ext.spherical_model import *
from amuse.units import units


class TestUniformSphericalDistribution(TestCase):
    
    def test1(self):
        instance = UniformSphericalDistribution(4200)
        x, y, z = instance.result
        self.assertEqual(len(x), 4200)
        r_squared = x*x + y*y + z*z
        self.assertAlmostEqual(r_squared.max(), 1.0)
        self.assertAlmostEqual(r_squared.min(), 0.0, 1)
    
    def test2(self):
        for n_i in [1003, 3210]:
            for type_i in ["cubic", "bcc", "body_centered_cubic", "random"]:
                for offset_i in [(0.406645, 0.879611, 0.573737), (0.939868, 0.796048, 0.236403)]:
                    numpy.random.seed(12345)
                    instance = UniformSphericalDistribution(n_i, type=type_i, offset=offset_i)
                    x, y, z = instance.result
                    self.assertEqual(len(x), n_i)
                    r_squared = x*x + y*y + z*z
                    self.assertAlmostEqual(r_squared.max(), 1.0, 2)
                    self.assertAlmostEqual(r_squared.min(), 0.0, 1)
    
    def test3(self):
        instance = UniformSphericalDistribution(1234, type="cubic", offset=(0.07974498,  0.77741132,  0.2993995))
        x, y, z = instance.result
        grid_spacing_x = min(x[numpy.where(x > min(x))]) - min(x)
        grid_spacing_y = min(y[numpy.where(y > min(y))]) - min(y)
        grid_spacing_z = min(z[numpy.where(z > min(z))]) - min(z)
        self.assertAlmostEqual(grid_spacing_x, grid_spacing_y)
        self.assertAlmostEqual(grid_spacing_x, grid_spacing_z)
        r_squared = x*x + y*y + z*z
        min_r_squared = sum([(grid_spacing_x * min(x,1-x))**2 for x in instance.offset])
        self.assertAlmostEqual(min(r_squared), min_r_squared)
     
    def test4(self):
        instance = UniformSphericalDistribution(1234, type="bcc", offset=(0.07974498,  0.77741132,  0.2993995))
        x, y, z = instance.result
        grid_spacing_x = min(x[numpy.where(x > min(x))]) - min(x)
        grid_spacing_y = min(y[numpy.where(y > min(y))]) - min(y)
        grid_spacing_z = min(z[numpy.where(z > min(z))]) - min(z)
        self.assertAlmostEqual(grid_spacing_x, grid_spacing_y)
        self.assertAlmostEqual(grid_spacing_x, grid_spacing_z)
        r_squared = x*x + y*y + z*z
        min_r_squared = min(
            sum([(2*grid_spacing_x * min(x,1-x))**2 for x in instance.offset]),
            sum([(2*grid_spacing_x * abs(x-.5))**2 for x in instance.offset])
        )
        self.assertAlmostEqual(min(r_squared), min_r_squared)
    
    def test5(self):
        numpy.random.seed(12345)
        theta = numpy.arccos(numpy.random.uniform(-1.0, 1.0))
        phi   = numpy.random.uniform(0.0, 2.0*numpy.pi)
        self.assertRaises(AmuseWarning, UniformSphericalDistribution, 1234, type="bcc", rotate=(theta, phi), 
            expected_message = "rotate is not yet supported")
    
    def test6(self):
        print "Test new_uniform_spherical_particle_distribution"
        particles = new_uniform_spherical_particle_distribution(14321, 1 | units.m, 1 | units.kg, type="cubic")
        self.assertEqual(len(particles), 14321)
        r_squared = particles.position.lengths_squared()
        self.assertAlmostEqual(r_squared.amax(), 1.0 | units.m**2)
        self.assertAlmostEqual(r_squared.amin(), 0.0 | units.m**2, 2)
        self.assertAlmostEqual(particles.total_mass(), 1 | units.kg)
        select_half_radius = numpy.where(r_squared < (0.5 | units.m)**2)
        self.assertAlmostEqual(particles.mass[select_half_radius].sum(), 1.0/8.0 | units.kg, places=3)
        self.assertAlmostEqual(particles.center_of_mass(), [0.0, 0.0, 0.0] | units.m, places=2)
    
    def test7(self):
        print "Test new_spherical_particle_distribution with total_mass specified"
        particles = new_spherical_particle_distribution(4200, radii = [2,4,3,1] | units.m, 
            densities = [80, 30, 50, 100] | (units.kg/units.m**3), total_mass = 10000 | units.kg)
        self.assertEqual(len(particles), 4200)
        self.assertAlmostEqual(particles.total_mass(), 10000 | units.kg)
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < (4.0 | units.m)**2))
        self.assertFalse(numpy.all( r_squared < (3.5 | units.m)**2 ))
        self.assertAlmostEqual(particles.center_of_mass(), [0.0, 0.0, 0.0] | units.m, places=2)
    
    def test8(self):
        print "Test new_spherical_particle_distribution without radii, densities tables"
        
        def my_density_func(radius):
            return (1 + radius.value_in(units.m)**2)**-2.5 | units.kg/units.m**3
        
        self.assertRaises(AmuseException, new_spherical_particle_distribution, 42, 
            radial_density_func = my_density_func, total_mass = 10000 | units.kg, expected_message = 
            "Using an arbitrary radial density function is not yet supported. Radius and density tables must be passed instead.")
    
    def test9(self):
        print "Test new_spherical_particle_distribution without total_mass"
        rad = [2,4,3,1] | units.m
        rho = [80, 30, 50, 100] | (units.kg/units.m**3)
        
        particles = new_spherical_particle_distribution(142, radii = rad, 
            densities = rho, size = 1.5 | units.m)
        self.assertEqual(len(particles), 142)
        self.assertAlmostEqual(particles.total_mass(), 
            numpy.pi * 4.0/3.0 *(100 * 1.0**3 + 80*(1.5**3-1.0**3)) | units.kg)
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < (1.5 | units.m)**2))
        self.assertFalse(numpy.all( r_squared < (1.4 | units.m)**2 ))
        
        particles = new_spherical_particle_distribution(14200, radii = rad, densities = rho)
        self.assertEqual(len(particles), 14200)
        interpolator = EnclosedMassInterpolator(rad, rho)
        self.assertAlmostEqual(particles.total_mass(), interpolator.get_enclosed_mass(max(rad)))
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < max(rad)**2 ))
        self.assertFalse(numpy.all( r_squared < (0.9*max(rad))**2 ))
        self.assertAlmostEqual(particles.center_of_mass(), [0.0, 0.0, 0.0] | units.m, places=2)
        
    def test10(self):
        print "Test new_uniform_spherical_particle_distribution, glass"
        numpy.random.seed(12345)
        # setting target_rms to 30% for test speed-up
        particles = new_uniform_spherical_particle_distribution(1421, 1|units.m, 1|units.kg, 
            type="glass", target_rms=0.3)
        self.assertEqual(len(particles), 1421)
        r_squared = particles.position.lengths_squared()
        self.assertAlmostEqual(r_squared.amax(), 1.0 | units.m**2)
        self.assertAlmostEqual(r_squared.amin(), 0.0 | units.m**2, 1)
    
    def test11(self):
        print "Test new_uniform_spherical_particle_distribution, sobol sequence"
        particles = new_uniform_spherical_particle_distribution(14321, 1 | units.m, 1 | units.kg, type="sobol")
        self.assertEqual(len(particles), 14321)
        r_squared = particles.position.lengths_squared()
        self.assertAlmostEqual(r_squared.amax(), 1.0 | units.m**2)
        self.assertAlmostEqual(r_squared.amin(), 0.0 | units.m**2, 2)
        self.assertAlmostEqual(particles.total_mass(), 1 | units.kg)
        select_half_radius = numpy.where(r_squared < (0.5 | units.m)**2)
        self.assertAlmostEqual(particles.mass[select_half_radius].sum(), 1.0/8.0 | units.kg, places=2)
        self.assertAlmostEqual(particles.center_of_mass(), [0.0, 0.0, 0.0] | units.m, places=2)
    
    def test12(self):
        print "Test new_uniform_spherical_particle_distribution, face-centered cubic"
        particles = new_uniform_spherical_particle_distribution(14321, 1|units.m, 1|units.kg, type="fcc")
        self.assertEqual(len(particles), 14321)
        r_squared = particles.position.lengths_squared()
        self.assertAlmostEqual(r_squared.amax(), 1.0 | units.m**2)
        self.assertAlmostEqual(r_squared.amin(), 0.0 | units.m**2, 2)
        self.assertAlmostEqual(particles.total_mass(), 1 | units.kg)
        select_half_radius = numpy.where(r_squared < (0.5 | units.m)**2)
        self.assertAlmostEqual(particles.mass[select_half_radius].sum(), 1.0/8.0 | units.kg, places=2)
        self.assertAlmostEqual(particles.center_of_mass(), [0.0, 0.0, 0.0] | units.m, places=2)
    

class TestEnclosedMassInterpolator(TestCase):
    
    def test1(self):
        instance = EnclosedMassInterpolator()
        self.assertFalse(instance.initialized)
        instance.initialize([2,4,3,1] | units.m, [80, 30, 50, 100] | (units.kg/units.m**3))
        self.assertTrue(instance.initialized)
        self.assertEqual(instance.radii, [0,1,2,3,4] | units.m)
        self.assertEqual(instance.densities, [100,80,50,30] | (units.kg/units.m**3))
        self.assertEqual(instance.enclosed_mass[0], 0.0 | units.kg)
        self.assertEqual(instance.enclosed_mass[1], 
            numpy.pi * 4.0/3.0 * instance.densities[0] * instance.radii[1]**3)
        self.assertEqual(instance.get_enclosed_mass(1.0 | units.m), 
            numpy.pi * 4.0/3.0 * instance.densities[0] * (1.0 | units.m)**3)
        self.assertEqual(instance.get_enclosed_mass(0.5 | units.m), 
            numpy.pi * 4.0/3.0 * instance.densities[0] * (0.5 | units.m)**3)
        self.assertEqual(instance.get_enclosed_mass(1.5 | units.m), 
            numpy.pi * 4.0/3.0 * (100*1.0**3 + 80*(1.5**3-1.0**3)) | units.kg)
        del instance
    
    def test2(self):
        interpolator = EnclosedMassInterpolator()
        self.assertRaises(AmuseException, interpolator.get_enclosed_mass, 0.0 | units.m, 
            expected_message = "Can't calculate enclosed mass: interpolator is not initialized")
        interpolator.initialize([2,4,3,1] | units.m, [80, 30, 50, 100] | (units.kg/units.m**3))
        self.assertEqual(numpy.pi * 4.0/3.0 * 100 * 0.3**3 | units.kg, 
            interpolator.get_enclosed_mass(0.3 | units.m))
        self.assertEqual([0.0, numpy.pi * 4.0/3.0 * 100 * 0.3**3] | units.kg, 
            interpolator.get_enclosed_mass([0.0, 0.3] | units.m))
        self.assertEqual(numpy.pi * 4.0/3.0 * 100 * 1.0**3 | units.kg, 
            interpolator.get_enclosed_mass(1.0 | units.m))
        self.assertEqual(numpy.pi * 4.0/3.0 *(100 * 1.0**3 + 80*(1.5**3-1.0**3)) | units.kg, 
            interpolator.get_enclosed_mass(1.5 | units.m))
        self.assertRaises(AmuseException, interpolator.get_enclosed_mass, -0.5 | units.m, expected_message = 
            "Can't find a valid index. [-0.5] m is not in the range [0.0 m, 4.0 m].")
        self.assertRaises(AmuseException, interpolator.get_enclosed_mass, 4.5 | units.m, expected_message = 
            "Can't find a valid index. [4.5] m is not in the range [0.0 m, 4.0 m].")
        self.assertRaises(AmuseException, interpolator.get_enclosed_mass, [2.5, 3.5, 4.5] | units.m, expected_message = 
            "Can't find a valid index. [4.5] m is not in the range [0.0 m, 4.0 m].")
        self.assertRaises(AmuseException, interpolator.get_enclosed_mass, [-0.5, 3.5, 4.5] | units.m, expected_message = 
            "Can't find a valid index. [-0.5, 4.5] m is not in the range [0.0 m, 4.0 m].")
        del interpolator
    
    def test3(self):
        interpolator = EnclosedMassInterpolator()
        self.assertRaises(AmuseException, interpolator.get_radius_for_enclosed_mass, 0.0 | units.kg, 
            expected_message = "Can't calculate radius for enclosed mass: interpolator is not initialized")
        interpolator.initialize([2,4,3,1] | units.m, [80, 30, 50, 100] | (units.kg/units.m**3))
        self.assertEqual(0.3 | units.m, 
            interpolator.get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 * 100 * 0.3**3 | units.kg))
        self.assertEqual([0.0, 0.3] | units.m, 
            interpolator.get_radius_for_enclosed_mass([0.0, numpy.pi * 4.0/3.0 * 100 * 0.3**3] | units.kg))
        self.assertEqual(1.0 | units.m, 
            interpolator.get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 * 100 * 1.0**3 | units.kg))
        self.assertEqual(1.5 | units.m, 
            interpolator.get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 *(100 * 1.0**3 + 80*(1.5**3-1.0**3)) | units.kg))
        self.assertRaises(AmuseException, interpolator.get_radius_for_enclosed_mass, -0.5 | units.kg, expected_message = 
            "Can't find a valid index. [-0.5] kg is not in the range [0.0 kg, 11393.509357 kg].")
        self.assertRaises(AmuseException, interpolator.get_radius_for_enclosed_mass, 12000 | units.kg, expected_message = 
            "Can't find a valid index. [12000] kg is not in the range [0.0 kg, 11393.509357 kg].")
        self.assertRaises(AmuseException, interpolator.get_radius_for_enclosed_mass, [-0.5, 1.0, 12000] | units.kg, expected_message = 
            "Can't find a valid index. [-0.5, 12000.0] kg is not in the range [0.0 kg, 11393.509357 kg].")
        del interpolator
    
