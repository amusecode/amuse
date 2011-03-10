import os.path
from amuse.test.amusetest import TestWithMPI
#from amuse.support.exceptions import AmuseException
from amuse.support.data.core import Particles, Particle
from amuse.support.units import units
#from amuse.legacy.mesa.interface import MESA
from amuse.community.mmams.interface import MakeMeAMassiveStarInterface, MakeMeAMassiveStar

# Change the default for some MakeMeAMassiveStar(-Interface) keyword arguments:
default_options = dict(redirection="none")

class TestMakeMeAMassiveStarInterface(TestWithMPI):
    
    def test1(self):
        print "Test 1: initialization of the interface"
        instance = MakeMeAMassiveStarInterface(**default_options)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        instance.stop()
    
    def test2(self):
        print "Test 2: define a new particle"
        instance = MakeMeAMassiveStarInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        id, error = instance.new_particle(1.0)
        self.assertEqual(error, 0)
        self.assertEqual(id, 0)
        id, error = instance.new_particle([2.0, 3.0, 4.0])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(id, [1, 2, 3])
        n_particles, error = instance.get_number_of_particles()
        self.assertEqual(error, 0)
        self.assertEqual(n_particles, 4)
        
        error = instance.add_shell([1, 1], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0], 
            [5.0, 10.0], [6.0, 12.0], [7.0, 14.0], [8.0, 16.0], [9.0, 18.0], 
            [0.4, 0.2], [0.2, 0.4], [0.15, 0.1], [0.1, 0.15], [0.05, 0.01], 
            [0.04, 0.02], [0.03, 0.03], [0.02, 0.04], [0.01, 0.05])
        self.assertEqual(error, [0, 0])
        
        number_of_shells, error = instance.get_number_of_zones(1)
        self.assertEqual(error, 0)
        self.assertEqual(number_of_shells, 2)
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, C12, N14, O16, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_stellar_model_element([0, 1], [1, 1])
        self.assertEqual(error, [0, 0])
        self.assertEqual(mass,      [2.0, 4.0])
        self.assertEqual(radius,    [3.0, 6.0])
        self.assertEqual(density,   [4.0, 8.0])
        self.assertEqual(pressure,  [5.0, 10.0])
        self.assertEqual(e_thermal, [6.0, 12.0])
        self.assertEqual(entropy,   [7.0, 14.0])
        self.assertEqual(temperature, [8.0, 16.0])
        self.assertEqual(molecular_weight, [9.0, 18.0])
        self.assertEqual(H1,   [0.4,  0.2])
        self.assertEqual(He4,  [0.2,  0.4])
        self.assertEqual(C12,  [0.15, 0.1])
        self.assertEqual(N14,  [0.1,  0.15])
        self.assertEqual(O16,  [0.05, 0.01])
        self.assertEqual(Ne20, [0.04, 0.02])
        self.assertEqual(Mg24, [0.03, 0.03])
        self.assertEqual(Si28, [0.02, 0.04])
        self.assertEqual(Fe56, [0.01, 0.05])
        instance.stop()
    
    def test3(self):
        print "Test 3: read a new particle from a usm file"
        instance = MakeMeAMassiveStarInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        usm_file = os.path.join(instance.data_directory, 'primary.usm')
        id, error = instance.read_usm(usm_file)
        self.assertEqual(error, 0)
        self.assertEqual(id, 0)
        id, error = instance.new_particle([2.0, 3.0])
        self.assertEqual(error, [0, 0])
        self.assertEqual(id, [1, 2])
        
        n_particles, error = instance.get_number_of_particles()
        self.assertEqual(error, 0)
        self.assertEqual(n_particles, 3)
        
        number_of_shells, error = instance.get_number_of_zones([0, 1, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(number_of_shells, [187, 0, 0])
        
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, C12, N14, O16, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_stellar_model_element([0, 186, 0, 0], [0, 0, 1, 3])
        self.assertEqual(error, [0, 0, -2, -1])
        self.assertAlmostEqual(mass[0],  0.0, 3)
        self.assertAlmostEqual(mass[1], 20.0, 0)
        self.assertAlmostEqual(radius[0], 0.0, 1)
        self.assertAlmostEqual(radius[1], 16.8, 1)
        self.assertAlmostEqual(temperature[0], 47318040.0, 0)
        self.assertAlmostEqual(temperature[1], 81542.0, 0)
        self.assertAlmostEqual(H1[0], 0.0121, 4)
        self.assertAlmostEqual(H1[1], 0.7, 4)
        instance.stop()
    
    def slowtest4(self):
        print "Test 4: merge particles (from usm files)"
        instance = MakeMeAMassiveStarInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        usm_file = os.path.join(instance.data_directory, 'primary.usm')
        id, error = instance.read_usm(usm_file)
        self.assertEqual(error, 0)
        self.assertEqual(id, 0)
        
        usm_file = os.path.join(instance.data_directory, 'secondary.usm')
        id, error = instance.read_usm(usm_file)
        self.assertEqual(error, 0)
        self.assertEqual(id, 1)
        
        id, error = instance.merge_two_stars(0, 1)
        self.assertEqual(error, 0)
        self.assertEqual(id, 2)
        
        n_particles, error = instance.get_number_of_particles()
        self.assertEqual(error, 0)
        self.assertEqual(n_particles, 3)
        
        number_of_shells, error = instance.get_number_of_zones([0, 1, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(number_of_shells, [187, 181, 17602])
        
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, C12, N14, O16, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_stellar_model_element([0, 10000, 17601], [2, 2, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertAlmostEqual(mass,  [0.0, 21.8369, 25.675], 3)
        self.assertAlmostEqual(radius, [0.0, 6.456, 19.458], 3)
        self.assertAlmostEqual(temperature, [39054497.9, 6788317.3, 11.8], 0)
        self.assertAlmostEqual(H1, [0.61566, 0.69942, 0.70002], 4)
        instance.stop()
    

class TestMakeMeAMassiveStar(TestWithMPI):
    
    def test1(self):
        print "Test 1: initialization of the interface"
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.recommit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: define a new particle"
        stars = Particles(4)
        stars.mass = [1.0, 2.0, 3.0, 4.0] | units.MSun
        
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particle(stars[0])
        instance.particles.add_particles(stars[1:])
        self.assertEqual(instance.number_of_particles, 4)
        
        instance.particles[1].add_shell([2.0, 4.0]|units.MSun, [3.0, 6.0]|units.RSun, 
            [4.0, 8.0]|units.g / units.cm**3, [5.0, 10.0]|units.barye, 
            [6.0, 12.0]|units.erg / units.g, [7.0, 14.0]|units.none, 
            [8.0, 16.0]|units.K, [9.0, 18.0]|units.none, [0.4, 0.2]|units.none, 
            [0.2, 0.4]|units.none, [0.15, 0.1]|units.none, [0.1, 0.15]|units.none, 
            [0.05, 0.01]|units.none, [0.04, 0.02]|units.none, [0.03, 0.03]|units.none, 
            [0.02, 0.04]|units.none, [0.01, 0.05]|units.none)
        self.assertEqual(instance.particles[1].number_of_zones, 2)
        stellar_model = instance.particles[1].internal_structure()
        self.assertTrue(set(['mass', 'radius', 'rho', 'pressure', 'e_thermal', 
            'entropy', 'temperature', 'molecular_weight', 'X_H', 'X_He', 'X_C', 
            'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']).issubset(
                stellar_model.all_attributes()
            )
        )
        self.assertEqual(stellar_model.mass,      [2.0, 4.0]|units.MSun)
        self.assertEqual(stellar_model.radius,    [3.0, 6.0]|units.RSun)
        self.assertEqual(stellar_model.rho,       [4.0, 8.0]|units.g / units.cm**3)
        self.assertEqual(stellar_model.pressure,  [5.0, 10.0]|units.barye)
        self.assertEqual(stellar_model.e_thermal, [6.0, 12.0]|units.erg / units.g)
        self.assertEqual(stellar_model.entropy,   [7.0, 14.0]|units.none)
        self.assertEqual(stellar_model.temperature, [8.0, 16.0]|units.K)
        self.assertEqual(stellar_model.molecular_weight, [9.0, 18.0]|units.none)
        self.assertEqual(stellar_model.X_H,   [0.4,  0.2]|units.none)
        self.assertEqual(stellar_model.X_He,  [0.2,  0.4]|units.none)
        self.assertEqual(stellar_model.X_C,  [0.15, 0.1]|units.none)
        self.assertEqual(stellar_model.X_N,  [0.1,  0.15]|units.none)
        self.assertEqual(stellar_model.X_O,  [0.05, 0.01]|units.none)
        self.assertEqual(stellar_model.X_Ne, [0.04, 0.02]|units.none)
        self.assertEqual(stellar_model.X_Mg, [0.03, 0.03]|units.none)
        self.assertEqual(stellar_model.X_Si, [0.02, 0.04]|units.none)
        self.assertEqual(stellar_model.X_Fe, [0.01, 0.05]|units.none)
        instance.stop()
    
    def test3(self):
        print "Test 3: read a new particle from a usm file"
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        stars = Particles(4)
        stars.usm_file = [os.path.join(instance.data_directory, filename) for 
            filename in ['primary.usm', 'secondary.usm', '', '']] | units.string
        stars[2:].mass = [3.0, 4.0] | units.MSun
        instance.imported_stars.add_particles(stars[:2])
        
        instance.particles.add_particles(stars[2:])
        self.assertEqual(instance.number_of_particles, 4)
        self.assertEqual(instance.native_stars.number_of_zones, [0, 0])
        self.assertEqual(instance.imported_stars.number_of_zones, [187, 181])
        self.assertEqual(instance.particles.number_of_zones, [0, 0, 187, 181])
        
        stellar_model = instance.imported_stars[0].internal_structure()
        self.assertAlmostEqual(stellar_model.mass[0],  0.0 | units.MSun, 3)
        self.assertAlmostEqual(stellar_model.mass[-1], 20.0 | units.MSun, 0)
        self.assertAlmostEqual(stellar_model.radius[0], 0.0 | units.RSun, 1)
        self.assertAlmostEqual(stellar_model.radius[-1], 16.8 | units.RSun, 1)
        self.assertAlmostEqual(stellar_model.temperature[0], 47318040.0 | units.K, 0)
        self.assertAlmostEqual(stellar_model.temperature[-1], 81542.0 | units.K, 0)
        self.assertAlmostEqual(stellar_model.X_H[0], 0.0121 | units.none, 4)
        self.assertAlmostEqual(stellar_model.X_H[-1], 0.7 | units.none, 4)
        instance.stop()
    
    def slowtest4(self):
        print "Test 4: merge particles (from usm files)"
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        stars = Particles(2)
        stars.usm_file = [os.path.join(instance.data_directory, filename) for 
            filename in ['primary.usm', 'secondary.usm']] | units.string
        instance.imported_stars.add_particles(stars)
        
        merge_product = Particle()
        merge_product.primary = instance.imported_stars[0]
        merge_product.secondary = instance.imported_stars[1]
        instance.merge_products.add_particle(merge_product)
        self.assertEqual(instance.number_of_particles, 3)
        self.assertEqual(instance.particles.number_of_zones, [187, 181, 17602])
        
        stellar_model = instance.merge_products[0].internal_structure()
        self.assertAlmostEqual(stellar_model.mass[[0, 10000, 17601]],  [0.0, 21.8369, 25.675] | units.MSun, 3)
        self.assertAlmostEqual(stellar_model.radius[[0, 10000, 17601]], [0.0, 6.456, 19.458] | units.RSun, 3)
        self.assertAlmostEqual(stellar_model.temperature[[0, 10000, 17601]], [39054497.9, 6788317.3, 11.8] | units.K, 0)
        self.assertAlmostEqual(stellar_model.X_H[[0, 10000, 17601]], [0.61566, 0.69942, 0.70002] | units.none, 4)
        instance.stop()
   

