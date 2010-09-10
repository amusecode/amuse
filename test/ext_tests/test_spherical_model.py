import numpy
from amuse.test.amusetest import TestCase

from amuse.support.exceptions import AmuseWarning, AmuseException
from amuse.support.units import units
from amuse.ext.spherical_model import *


class TestUniformSphericalDistribution(TestCase):
    
    def test1(self):
        instance = UniformSphericalDistribution(42)
        x, y, z = instance.result
        self.assertEqual(len(x), 42)
        r_squared = x*x + y*y + z*z
        self.assertTrue(numpy.all( r_squared < 1.0 ))
        self.assertFalse(numpy.all( r_squared < 0.9**2 ))
    
    def test2(self):
        for n_i in [42, 103, 321]:
            for type_i in ["cubic", "bcc", "body_centered_cubic", "random"]:
                for offset_i in [(0.406645, 0.879611, 0.573737), (0.939868, 0.796048, 0.236403)]:
                    numpy.random.seed(12345)
                    instance = UniformSphericalDistribution(n_i, type=type_i, offset=offset_i)
                    x, y, z = instance.result
                    self.assertEqual(len(x), n_i)
                    r_squared = x*x + y*y + z*z
                    self.assertTrue(numpy.all( r_squared < 1.0**2 ))
                    self.assertFalse(numpy.all(r_squared < 0.9**2 ))
    
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
        particles = new_uniform_spherical_particle_distribution(1421, 1 | units.m, 1 | units.kg, type="cubic")
        self.assertEqual(len(particles), 1421)
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < (1.0 | units.m)**2))
        self.assertFalse(numpy.all( r_squared < (0.9 | units.m)**2 ))
        self.assertAlmostEqual(particles.total_mass(), 1 | units.kg)
        select_half_radius = numpy.where(r_squared < (0.5 | units.m)**2)
        self.assertAlmostEqual(particles.mass[select_half_radius].sum(), 1.0/8.0 | units.kg, places=3)
    
    def test7(self):
        print "Test new_spherical_particle_distribution with total_mass specified"
        particles = new_spherical_particle_distribution(42, radii = [2,4,3,1] | units.m, 
            densities = [80, 30, 50, 100] | (units.kg/units.m**3), total_mass = 10000 | units.kg)
        self.assertEqual(len(particles), 42)
        self.assertAlmostEqual(particles.total_mass(), 10000 | units.kg)
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < (4.0 | units.m)**2))
        self.assertFalse(numpy.all( r_squared < (3.5 | units.m)**2 ))
    
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
        
        particles = new_spherical_particle_distribution(142, radii = rad, densities = rho)
        self.assertEqual(len(particles), 142)
        self.assertAlmostEqual(particles.total_mass(), 
            get_enclosed_mass_from_tabulated(max(rad), radii = rad, densities = rho))
        r_squared = particles.position.lengths_squared()
        self.assertTrue(numpy.all( r_squared < max(rad)**2 ))
        self.assertFalse(numpy.all( r_squared < (0.9*max(rad))**2 ))
        
    

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
    
    def test2(self):
        self.assertRaises(AmuseException, get_enclosed_mass_from_tabulated, 0.0 | units.m, expected_message = 
            "Interpolator is not initialized. Radius and density tables must be passed in the first call.")
        self.assertEqual(0.0 | units.kg, get_enclosed_mass_from_tabulated(0.0 | units.m, 
            [2,4,3,1] | units.m, [80, 30, 50, 100] | (units.kg/units.m**3)))
        self.assertEqual(numpy.pi * 4.0/3.0 * 100 * 0.3**3 | units.kg, 
            get_enclosed_mass_from_tabulated(0.3 | units.m))
        self.assertEqual(numpy.pi * 4.0/3.0 * 100 * 1.0**3 | units.kg, 
            get_enclosed_mass_from_tabulated(1.0 | units.m))
        self.assertEqual(numpy.pi * 4.0/3.0 *(100 * 1.0**3 + 80*(1.5**3-1.0**3)) | units.kg, 
            get_enclosed_mass_from_tabulated(1.5 | units.m))
        self.assertRaises(AmuseException, get_enclosed_mass_from_tabulated, -0.5 | units.m, expected_message = 
            "Can't find a valid index. -0.5 m is not in the range [0.0 m, 4.0 m].")
        self.assertRaises(AmuseException, get_enclosed_mass_from_tabulated, 4.5 | units.m, expected_message = 
            "Can't find a valid index. 4.5 m is not in the range [0.0 m, 4.0 m].")
    
    def test3(self):
        self.assertRaises(AmuseException, get_radius_for_enclosed_mass, 0.0 | units.kg, expected_message = 
            "Interpolator is not initialized. Radius and density tables must be passed in the first call.")
        self.assertEqual(0.0 | units.m, get_radius_for_enclosed_mass(0.0 | units.kg, 
            [2,4,3,1] | units.m, [80, 30, 50, 100] | (units.kg/units.m**3)))
        self.assertEqual(0.3 | units.m, 
            get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 * 100 * 0.3**3 | units.kg))
        self.assertEqual(1.0 | units.m, 
            get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 * 100 * 1.0**3 | units.kg))
        self.assertEqual(1.5 | units.m, 
            get_radius_for_enclosed_mass(numpy.pi * 4.0/3.0 *(100 * 1.0**3 + 80*(1.5**3-1.0**3)) | units.kg))
        self.assertRaises(AmuseException, get_radius_for_enclosed_mass, -0.5 | units.kg, expected_message = 
            "Can't find a valid index. -0.5 kg is not in the range [0.0 kg, 11393.509357 kg].")
        self.assertRaises(AmuseException, get_radius_for_enclosed_mass, 12000 | units.kg, expected_message = 
            "Can't find a valid index. 12000 kg is not in the range [0.0 kg, 11393.509357 kg].")
    
