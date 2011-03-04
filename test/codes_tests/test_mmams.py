import os.path
from amuse.test.amusetest import TestWithMPI
#from amuse.support.exceptions import AmuseException
#from amuse.support.data import core
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
        
        number_of_shells, error = instance.get_number_of_shells(1)
        self.assertEqual(error, 0)
        self.assertEqual(number_of_shells, 2)
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, O16, N14, C12, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_shell([1, 1], [0, 1])
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
        self.assertEqual(O16,  [0.15, 0.1])
        self.assertEqual(N14,  [0.1,  0.15])
        self.assertEqual(C12,  [0.05, 0.01])
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
        
        number_of_shells, error = instance.get_number_of_shells([0, 1, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(number_of_shells, [187, 0, 0])
        
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, O16, N14, C12, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_shell([0, 0, 1, 3], [0, 186, 0, 0])
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
        
        number_of_shells, error = instance.get_number_of_shells([0, 1, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(number_of_shells, [187, 181, 17602])
        
        mass, radius, density, pressure, e_thermal, entropy, temperature, \
            molecular_weight, H1, He4, O16, N14, C12, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_shell([2, 2, 2], [0, 10000, 17601])
        self.assertEqual(error, [0, 0, 0])
        self.assertAlmostEqual(mass,  [0.0, 21.8369, 25.675], 3)
        self.assertAlmostEqual(radius, [0.0, 6.456, 19.458], 3)
        self.assertAlmostEqual(temperature, [39054497.9, 6788317.3, 11.8], 0)
        self.assertAlmostEqual(H1, [0.61566, 0.69942, 0.70002], 4)
        instance.stop()
   

class TestMakeMeAMassiveStar(TestWithMPI):
    
    def test1(self):
        print "Testing initialization of the interface..."
        instance = MakeMeAMassiveStar(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        instance.stop()
    

