import os.path
import itertools
import pickle
from amuse.test.amusetest import get_path_to_results, TestWithMPI
from amuse.support.exceptions import AmuseException, CodeException
from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mmams.interface import MakeMeAMassiveStarInterface, MakeMeAMassiveStar
from amuse.units import units
from amuse.datamodel import Particles, Particle

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
        
        error = instance.add_shell([1, 1, 1, 1], [1.0, 2.0, 3.0, 3.0], 
            [2.0, 4.0, 6.0, 6.0], [3.0, 6.0, 9.0, 9.0], [4.0, 8.0, 12.0, 12.0], 
            [5.0, 10.0, 15.0, 15.0], [6.0, 12.0, 18.0, 18.0], 
            [-6.0, -12.0, -18.0, -18.0], [7.0, 14.0, 21.0, 21.0], 
            [0.4, 0.2, 0.4, 0.4], [0.2, 0.4, 0.2, 0.2], [0.15, 0.1, 0.1, 0.1], 
            [0.1, 0.15, 0.15, 0.15], [0.05, 0.01, 0.01, 0.01], 
            [0.04, 0.02, 0.02, 0.02], [0.03, 0.03, 0.03, 0.03], 
            [0.02, 0.04, 0.04, 0.04], [0.01, 0.05, 0.05, 0.05])
        self.assertEqual(error, [0, 0, 0, 0])
        
        number_of_shells, error = instance.get_number_of_zones(1)
        self.assertEqual(error, 0)
        self.assertEqual(number_of_shells, 3)
        d_mass, mass, radius, density, pressure, entropy, temperature, luminosity, \
            molecular_weight, H1, He4, C12, N14, O16, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_stellar_model_element([0, 1, 2], [1, 1, 1])
        self.assertEqual(error, [0, 0, 0])
        self.assertEqual(d_mass,    [1.0, 2.0, 3.0])
        self.assertEqual(mass,      [2.0, 4.0, 6.0])
        self.assertEqual(radius,    [3.0, 6.0, 9.0])
        self.assertEqual(density,   [4.0, 8.0, 12.0])
        self.assertEqual(pressure,  [5.0, 10.0, 15.0])
        self.assertEqual(temperature, [6.0, 12.0, 18.0])
        self.assertEqual(luminosity, [-6.0, -12.0, -18.0])
        self.assertEqual(molecular_weight, [7.0, 14.0, 21.0])
        self.assertEqual(H1,   [0.4,  0.2, 0.4])
        self.assertEqual(He4,  [0.2,  0.4, 0.2])
        self.assertEqual(C12,  [0.15, 0.1, 0.1])
        self.assertEqual(N14,  [0.1,  0.15, 0.15])
        self.assertEqual(O16,  [0.05, 0.01, 0.01])
        self.assertEqual(Ne20, [0.04, 0.02, 0.02])
        self.assertEqual(Mg24, [0.03, 0.03, 0.03])
        self.assertEqual(Si28, [0.02, 0.04, 0.04])
        self.assertEqual(Fe56, [0.01, 0.05, 0.05])
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
        
        d_mass, mass, radius, density, pressure, entropy, temperature, luminosity, \
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
        
        d_mass, mass, radius, density, pressure, entropy, temperature, luminosity, \
            molecular_weight, H1, He4, C12, N14, O16, Ne20, Mg24, Si28, Fe56, \
            error = instance.get_stellar_model_element([0, 10000, 17601], [2, 2, 2])
        self.assertEqual(error, [0, 0, 0])
        self.assertAlmostEqual(mass,  [0.0, 21.8369, 25.675], 3)
        self.assertAlmostEqual(radius, [0.0, 6.456, 19.458], 3)
        self.assertAlmostEqual(temperature, [39054497.9, 6788317.3, 11.8], 0)
        self.assertAlmostEqual(H1, [0.61566, 0.69942, 0.70002], 4)
        instance.stop()
    
    def test5(self):
        print "Test 5: parameters"
        instance = MakeMeAMassiveStarInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        
        dump_mixed_flag, error = instance.get_dump_mixed_flag()
        self.assertEqual(error, 0)
        self.assertEqual(dump_mixed_flag, 0)
        self.assertEqual(instance.set_dump_mixed_flag(1), 0)
        dump_mixed_flag, error = instance.get_dump_mixed_flag()
        self.assertEqual(error, 0)
        self.assertEqual(dump_mixed_flag, 1)
        
        target_n_shells_mixing, error = instance.get_target_n_shells_mixing()
        self.assertEqual(error, 0)
        self.assertEqual(target_n_shells_mixing, 200)
        self.assertEqual(instance.set_target_n_shells_mixing(300), 0)
        target_n_shells_mixing, error = instance.get_target_n_shells_mixing()
        self.assertEqual(error, 0)
        self.assertEqual(target_n_shells_mixing, 300)
        
        target_n_shells, error = instance.get_target_n_shells()
        self.assertEqual(error, 0)
        self.assertEqual(target_n_shells, 10000)
        self.assertEqual(instance.set_target_n_shells(5000), 0)
        target_n_shells, error = instance.get_target_n_shells()
        self.assertEqual(error, 0)
        self.assertEqual(target_n_shells, 5000)
        
        do_shock_heating_flag, error = instance.get_do_shock_heating_flag()
        self.assertEqual(error, 0)
        self.assertEqual(do_shock_heating_flag, 1)
        self.assertEqual(instance.set_do_shock_heating_flag(0), 0)
        do_shock_heating_flag, error = instance.get_do_shock_heating_flag()
        self.assertEqual(error, 0)
        self.assertEqual(do_shock_heating_flag, 0)
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
        
        instance.particles[1].add_shell([1.0, 2.0]|units.MSun, [2.0, 4.0]|units.MSun, [3.0, 6.0]|units.RSun, 
            [4.0, 8.0]|units.g / units.cm**3, [5.0, 10.0]|units.barye, 
            [6.0, 12.0]|units.K, [-6.0, -12.0]|units.LSun, [7.0, 14.0]|units.amu, [0.4, 0.2]|units.none, 
            [0.2, 0.4]|units.none, [0.15, 0.1]|units.none, [0.1, 0.15]|units.none, 
            [0.05, 0.01]|units.none, [0.04, 0.02]|units.none, [0.03, 0.03]|units.none, 
            [0.02, 0.04]|units.none, [0.01, 0.05]|units.none)
        self.assertEqual(instance.particles[1].number_of_zones, 2)
        stellar_model = instance.particles[1].internal_structure()
        self.assertTrue(set(['d_mass', 'mass', 'radius', 'rho', 'pressure', 'entropy', 
            'temperature', 'luminosity', 'molecular_weight', 'X_H', 'X_He', 'X_C', 
            'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']).issubset(
                stellar_model.all_attributes()
            )
        )
        self.assertEqual(stellar_model.d_mass,    [1.0, 2.0]|units.MSun)
        self.assertEqual(stellar_model.mass,      [2.0, 4.0]|units.MSun)
        self.assertEqual(stellar_model.radius,    [3.0, 6.0]|units.RSun)
        self.assertEqual(stellar_model.rho,       [4.0, 8.0]|units.g / units.cm**3)
        self.assertEqual(stellar_model.pressure,  [5.0, 10.0]|units.barye)
        self.assertAlmostEqual(stellar_model.entropy, [407038297.689, 256418059.686]|units.none, 2)
        self.assertEqual(stellar_model.temperature,   [6.0, 12.0]|units.K)
        self.assertEqual(stellar_model.luminosity,    [-6.0, -12.0]|units.LSun)
        self.assertEqual(stellar_model.molecular_weight, [7.0, 14.0]|units.amu)
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
#        instance = MakeMeAMassiveStar(debugger = 'gdb', **default_options)
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
    
    def test5(self):
        print "Test 5: parameters"
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        
        for par, value in [('dump_mixed_flag', False), ('do_shock_heating_flag', True)]:
            self.assertTrue(value is getattr(instance.parameters, par))
            setattr(instance.parameters, par, not value)
            self.assertFalse(value is getattr(instance.parameters, par))
        
        for par, value in [('target_n_shells_mixing', 200), ('target_n_shells', 10000)]:
            self.assertEquals(value | units.none, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 1 | units.none)
            self.assertEquals(1 | units.none, getattr(instance.parameters, par))
        
        instance.stop()
    
    def slowtest6(self):
        print "Test 6: MakeMeAMassiveStar with MESA particles - match composition only"
        stars = Particles(2)
        stars.mass = [20.0, 8.0] | units.MSun
        
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        self.assertEqual(instance.number_of_particles, 2)
        
        stellar_evolution = self.new_instance(MESA, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.initialize_code() 
        stellar_evolution.commit_parameters()
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(2 | units.Myr)
        
        stellar_model = [None, None]
        for i in [0, 1]:
            number_of_zones     = stellar_evolution.particles[i].get_number_of_zones()
            mass_profile        = stellar_evolution.particles[i].get_mass_profile(
                number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
            cumul_mass_profile  = stellar_evolution.particles[i].get_cumulative_mass_profile(
                number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
            density_profile     = stellar_evolution.particles[i].get_density_profile(number_of_zones = number_of_zones)
            radius_profile      = stellar_evolution.particles[i].get_radius_profile(number_of_zones = number_of_zones)
            temperature_profile = stellar_evolution.particles[i].get_temperature_profile(number_of_zones = number_of_zones)
            pressure_profile    = stellar_evolution.particles[i].get_pressure_profile(number_of_zones = number_of_zones)
            luminosity_profile  = stellar_evolution.particles[i].get_luminosity_profile(number_of_zones = number_of_zones)
            mu_profile          = stellar_evolution.particles[i].get_mu_profile(number_of_zones = number_of_zones)
            composition_profile = stellar_evolution.particles[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
            species_names       = stellar_evolution.particles[i].get_names_of_species()
            self.assertEquals(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
            
            instance.particles[i].add_shell(mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                composition_profile[1]+composition_profile[2], composition_profile[3], 
                composition_profile[4], composition_profile[5], composition_profile[6], 
                composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
            
            self.assertEqual(instance.particles[i].number_of_zones, number_of_zones)
            stellar_model[i] = instance.particles[i].internal_structure()
            self.assertTrue(set(['d_mass', 'mass', 'radius', 'rho', 'pressure', 
                'entropy', 'temperature', 'luminosity', 'molecular_weight', 'X_H', 'X_He', 'X_C', 
                'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']).issubset(
                    stellar_model[i].all_attributes()
                )
            )
        
        self.assertAlmostEqual(stellar_model[0].mass[0],  0.0 | units.MSun, 3)
        self.assertAlmostEqual(stellar_model[0].mass[-1], 20.0 | units.MSun, 0)
        self.assertAlmostEqual(stellar_model[0].radius[0], 0.0 | units.RSun, 1)
        self.assertAlmostEqual(stellar_model[0].radius[-1], 6.67 | units.RSun, 1)
        self.assertAlmostRelativeEqual(stellar_model[0].temperature[0], 35408447.7 | units.K, 2)
        self.assertAlmostRelativeEqual(stellar_model[0].temperature[-1], 33533.7 | units.K, 2)
        self.assertAlmostEqual(stellar_model[0].X_H[0], 0.58642 | units.none, 2)
        self.assertAlmostEqual(stellar_model[0].X_H[-1], 0.7 | units.none, 2)
        
        self.assertEqual(instance.number_of_particles, 2)
        merge_product = Particle()
        merge_product.primary = instance.particles[0]
        merge_product.secondary = instance.particles[1]
        instance.merge_products.add_particle(merge_product)
        self.assertEqual(instance.number_of_particles, 3)
        self.assertEqual(instance.particles.number_of_zones, [1043, 985, 18373])
        
        stellar_model = instance.merge_products[0].internal_structure()
        self.assertAlmostEqual(stellar_model.mass[[0, -1]],        [0.0, 25.7309] | units.MSun, 3)
        self.assertAlmostEqual(stellar_model.radius[[0, -1]],      [0.0,  8.353] | units.RSun, 2)
        self.assertAlmostRelativeEqual(stellar_model.temperature[[0, -1]], [41089300.4, 450625.4] | units.K, 1)
        self.assertAlmostEqual(stellar_model.X_H[[0, -1]],         [0.67437, 0.69999] | units.none, 2)
        
        merged = Particle()
        merged.mass = stellar_model.mass[-1]
        merged_in_code = stellar_evolution.native_stars.add_particle(merged)
        stellar_evolution.evolve_model(keep_synchronous = False)
        mass_profile = merged_in_code.get_cumulative_mass_profile()*merged.mass
        new_composition = instance.match_composition_to_mass_profile(stellar_model, mass_profile)
        instance.stop()
        merged_in_code.set_chemical_abundance_profiles(new_composition)
        for i in range(10):
            stellar_evolution.evolve_model(keep_synchronous = False)
            print stellar_evolution.particles
        
        stellar_evolution.stop()
    
    def slowtest7(self):
        print "Test 7: MakeMeAMassiveStar with MESA particles - import product into MESA"
        stars = Particles(2)
        stars.mass = [20.0, 8.0] | units.MSun
        
        stellar_evolution = self.new_instance(MESA, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.initialize_code() 
        stellar_evolution.commit_parameters()
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(2 | units.Myr)
        
        if os.path.exists(os.path.join(get_path_to_results(), "test_mmams_7.pkl")):
           # Remove the next line to speed-up, and only test MESA loading the merger product
           os.remove(os.path.join(get_path_to_results(), "test_mmams_7.pkl"))
        
        if os.path.exists(os.path.join(get_path_to_results(), "test_mmams_7.pkl")):
            with open(os.path.join(get_path_to_results(), "test_mmams_7.pkl"), 'r') as in_file:
                stellar_model = pickle.load(in_file)
        else:
            instance = MakeMeAMassiveStar(**default_options)
            instance.initialize_code()
            instance.parameters.target_n_shells_mixing = 1000
            instance.parameters.dump_mixed_flag = True
            instance.commit_parameters()
            instance.particles.add_particles(stars)
            self.assertEqual(instance.number_of_particles, 2)
            
            stellar_model = [None, None]
            for i in [0, 1]:
                number_of_zones     = stellar_evolution.particles[i].get_number_of_zones()
                mass_profile        = stellar_evolution.particles[i].get_mass_profile(
                    number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
                cumul_mass_profile  = stellar_evolution.particles[i].get_cumulative_mass_profile(
                    number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
                density_profile     = stellar_evolution.particles[i].get_density_profile(number_of_zones = number_of_zones)
                radius_profile      = stellar_evolution.particles[i].get_radius_profile(number_of_zones = number_of_zones)
                temperature_profile = stellar_evolution.particles[i].get_temperature_profile(number_of_zones = number_of_zones)
                pressure_profile    = stellar_evolution.particles[i].get_pressure_profile(number_of_zones = number_of_zones)
                luminosity_profile  = stellar_evolution.particles[i].get_luminosity_profile(number_of_zones = number_of_zones)
                mu_profile          = stellar_evolution.particles[i].get_mu_profile(number_of_zones = number_of_zones)
                composition_profile = stellar_evolution.particles[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
                species_names       = stellar_evolution.particles[i].get_names_of_species()
                self.assertEquals(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
                
                instance.particles[i].add_shell(mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                    pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                    composition_profile[1]+composition_profile[2], composition_profile[3], 
                    composition_profile[4], composition_profile[5], composition_profile[6], 
                    composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
                
                self.assertEqual(instance.particles[i].number_of_zones, number_of_zones)
                stellar_model[i] = instance.particles[i].internal_structure()
                self.assertTrue(set(['d_mass', 'mass', 'radius', 'rho', 'pressure', 
                    'entropy', 'temperature', 'luminosity', 'molecular_weight', 'X_H', 'X_He', 'X_C', 
                    'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']).issubset(
                        stellar_model[i].all_attributes()
                    )
                )
            
            self.assertAlmostEqual(stellar_model[0].mass[0],  0.0 | units.MSun, 3)
            self.assertAlmostEqual(stellar_model[0].mass[-1], 20.0 | units.MSun, 0)
            self.assertAlmostEqual(stellar_model[0].radius[0], 0.0 | units.RSun, 1)
            self.assertAlmostEqual(stellar_model[0].radius[-1], 6.67 | units.RSun, 1)
            self.assertAlmostRelativeEqual(stellar_model[0].temperature[0], 35408447.7 | units.K, 3)
            self.assertAlmostRelativeEqual(stellar_model[0].temperature[-1], 33533.7 | units.K, 2)
            self.assertAlmostEqual(stellar_model[0].X_H[0], 0.58642 | units.none, 2)
            self.assertAlmostEqual(stellar_model[0].X_H[-1], 0.7 | units.none, 4)
            
            self.assertEqual(instance.number_of_particles, 2)
            merge_product = Particle()
            merge_product.primary = instance.particles[0]
            merge_product.secondary = instance.particles[1]
            instance.merge_products.add_particle(merge_product)
            self.assertEqual(instance.number_of_particles, 3)
            self.assertEqual(instance.particles.number_of_zones, [1043, 985, 935])
            
            mmams_merged_model = instance.merge_products[0].internal_structure()
            self.assertAlmostEqual(mmams_merged_model.mass[[0, -1]],        [0.0259, 25.7055] | units.MSun, 3)
            self.assertAlmostEqual(mmams_merged_model.radius[[0, -1]],      [0.036,  7.574] | units.RSun, 2)
            self.assertAlmostRelativeEqual(mmams_merged_model.temperature[[0, -1]], [40844102.7, 1832217.8] | units.K, 1)
            self.assertAlmostEqual(mmams_merged_model.X_H[[0, -1]],         [0.6743, 0.6999] | units.none, 2)
            stellar_model = dict(
                mass        = mmams_merged_model.mass,
                radius      = mmams_merged_model.radius,
                rho         = mmams_merged_model.rho,
                temperature = mmams_merged_model.temperature,
                luminosity  = mmams_merged_model.luminosity,
                X_H  = mmams_merged_model.X_H,
                X_He = mmams_merged_model.X_He,
                X_C  = mmams_merged_model.X_C,
                X_N  = mmams_merged_model.X_N,
                X_O  = mmams_merged_model.X_O,
                X_Ne = mmams_merged_model.X_Ne,
                X_Mg = mmams_merged_model.X_Mg,
                X_Si = mmams_merged_model.X_Si,
                X_Fe = mmams_merged_model.X_Fe
            )
            instance.stop()
            with open(os.path.join(get_path_to_results(), "test_mmams_7.pkl"), 'w') as out_file:
                pickle.dump(stellar_model, out_file)
        
        stellar_evolution.new_particle_from_model(stellar_model, 10.0 | units.Myr)
        print stellar_evolution.particles
        for i in range(10):
            stellar_evolution.evolve_model(keep_synchronous = False)
            print stellar_evolution.particles
        stellar_evolution.stop()
    
    def slowtest8(self):
        print "Test 8: MakeMeAMassiveStar with MESA particles - multiple mergers"
        number_of_stars = 5
        stars = Particles(number_of_stars)
        stars.mass = range(20, 20+number_of_stars) | units.MSun
        
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
#        instance.parameters.target_n_shells = 20000
        instance.commit_parameters()
        
        stellar_evolution = self.new_instance(MESA, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.initialize_code() 
        stellar_evolution.commit_parameters()
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(2 | units.Myr)
        
        dummy = Particle()
        dummy.mass = 0 | units.MSun
        
        while len(stellar_evolution.particles) > 1:
            print "Selecting the two most massive stars."
            to_be_merged = [dummy, dummy]
            for particle in stellar_evolution.particles:
                if particle.mass > to_be_merged[0].mass:
                    if particle.mass > to_be_merged[1].mass:
                        to_be_merged[0] = to_be_merged[1]
                        to_be_merged[1] = particle
                    else:
                        to_be_merged[0] = particle
            print "Merging particles with mass", to_be_merged[0].mass, "and", to_be_merged[1].mass
            
            merge_product = Particle()
            for i in [0, 1]:
                number_of_zones     = to_be_merged[i].get_number_of_zones()
                mass_profile        = to_be_merged[i].get_mass_profile() * to_be_merged[i].mass
                cumul_mass_profile  = to_be_merged[i].get_cumulative_mass_profile() * to_be_merged[i].mass
                density_profile     = to_be_merged[i].get_density_profile(number_of_zones = number_of_zones)
                radius_profile      = to_be_merged[i].get_radius_profile(number_of_zones = number_of_zones)
                temperature_profile = to_be_merged[i].get_temperature_profile(number_of_zones = number_of_zones)
                pressure_profile    = to_be_merged[i].get_pressure_profile(number_of_zones = number_of_zones)
                luminosity_profile  = to_be_merged[i].get_luminosity_profile(number_of_zones = number_of_zones)
                mu_profile          = to_be_merged[i].get_mu_profile(number_of_zones = number_of_zones)
                composition_profile = to_be_merged[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
                species_names       = to_be_merged[i].get_names_of_species()
                
                new_mmams_particle = instance.native_stars.add_particle(to_be_merged[i])
                stellar_evolution.particles.remove_particle(to_be_merged[i])
                new_mmams_particle.add_shell(mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                    pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                    composition_profile[1]+composition_profile[2], composition_profile[3], 
                    composition_profile[4], composition_profile[5], composition_profile[6], 
                    composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
                self.assertEqual(new_mmams_particle.number_of_zones, number_of_zones)
                if i == 0:
                    merge_product.primary = new_mmams_particle
                else:
                    merge_product.secondary = new_mmams_particle
                
            self.assertEqual(len(instance.merge_products)+len(stellar_evolution.particles), 3)
            new_merge_product = instance.merge_products.add_particle(merge_product)
            print "Successfully merged particles. number_of_zones:", new_merge_product.number_of_zones
            self.assertEqual(len(instance.merge_products)+len(stellar_evolution.particles), 4)
            
            new_merge_product_model = new_merge_product.internal_structure()
            stellar_evolution.new_particle_from_model(new_merge_product_model, 42 | units.yr)
            print stellar_evolution.particles
            self.assertEqual(len(instance.merge_products)+len(stellar_evolution.particles), 5)
            for i in range(10):
                stellar_evolution.evolve_model(keep_synchronous = False)
                print stellar_evolution.particles

        stellar_evolution.stop()
        instance.stop()
    
    def slowtest9(self):
        print "Test 9: MakeMeAMassiveStar with MESA particles - evolved stars"
        stars = Particles(2)
        stars.mass = [20.0, 8.0] | units.MSun
        
        stellar_evolution = self.new_instance(MESA)#, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.initialize_code() 
        stellar_evolution.commit_parameters()
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        
        if os.path.exists(os.path.join(get_path_to_results(), "test_mmams_9.pkl")):
           # Remove the next line to speed-up, and only test MESA loading the merger product
           #os.remove(os.path.join(get_path_to_results(), "test_mmams_9.pkl"))
           pass
        
        if os.path.exists(os.path.join(get_path_to_results(), "test_mmams_9.pkl")):
            with open(os.path.join(get_path_to_results(), "test_mmams_9.pkl"), 'r') as in_file:
                stellar_model = pickle.load(in_file)
        else:
            instance = MakeMeAMassiveStar(**default_options)
            instance.initialize_code()
            instance.parameters.target_n_shells_mixing = 10000
            instance.parameters.dump_mixed_flag = True
            instance.commit_parameters()
            instance.particles.add_particles(stars)
            self.assertEqual(instance.number_of_particles, 2)
            
            try:
                while True:
                    stellar_evolution.particles[0].evolve_one_step()
            except AmuseException:
                stellar_evolution.particles[1].evolve_for(stellar_evolution.particles[0].age)
                print "Evolved stars to", stellar_evolution.particles.age
                print "Radius:", stellar_evolution.particles.radius
            
            stellar_model = [None, None]
            for i in [0, 1]:
                number_of_zones     = stellar_evolution.particles[i].get_number_of_zones()
                mass_profile        = stellar_evolution.particles[i].get_mass_profile(
                    number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
                cumul_mass_profile  = stellar_evolution.particles[i].get_cumulative_mass_profile(
                    number_of_zones = number_of_zones) * stellar_evolution.particles[i].mass
                density_profile     = stellar_evolution.particles[i].get_density_profile(number_of_zones = number_of_zones)
                radius_profile      = stellar_evolution.particles[i].get_radius_profile(number_of_zones = number_of_zones)
                temperature_profile = stellar_evolution.particles[i].get_temperature_profile(number_of_zones = number_of_zones)
                luminosity_profile  = stellar_evolution.particles[i].get_luminosity_profile(number_of_zones = number_of_zones)
                mu_profile          = stellar_evolution.particles[i].get_mu_profile(number_of_zones = number_of_zones)
                pressure_profile    = stellar_evolution.particles[i].get_pressure_profile(number_of_zones = number_of_zones)
                composition_profile = stellar_evolution.particles[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
                species_names       = stellar_evolution.particles[i].get_names_of_species()
                self.assertEquals(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 'si28', ''][:-(1+i)])
                
                instance.particles[i].add_shell(mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                    pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                    composition_profile[1]+composition_profile[2], composition_profile[3], 
                    composition_profile[4], composition_profile[5], composition_profile[6], 
                    composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
                
                self.assertEqual(instance.particles[i].number_of_zones, number_of_zones)
                stellar_model[i] = instance.particles[i].internal_structure()
                self.assertTrue(set(['d_mass', 'mass', 'radius', 'rho', 'pressure', 
                    'entropy', 'temperature', 'luminosity', 'molecular_weight', 'X_H', 'X_He', 'X_C', 
                    'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']).issubset(
                        stellar_model[i].all_attributes()
                    )
                )
            
            self.assertAlmostEqual(stellar_model[0].mass[0],  0.0 | units.MSun, 3)
            self.assertAlmostEqual(stellar_model[0].mass[-1], 20.0 | units.MSun)
            self.assertAlmostEqual(stellar_model[0].radius[0],     7.4820941908e-4 | units.RSun)
            self.assertAlmostEqual(stellar_model[0].radius[-1],    7.63066161669e2 | units.RSun)
            self.assertIsOfOrder(stellar_model[0].temperature[0],  3.66003952576e8 | units.K)
            self.assertIsOfOrder(stellar_model[0].temperature[-1], 3.67886047945e3 | units.K)
            self.assertAlmostEqual(stellar_model[0].X_H[0], 0.0 | units.none)
            self.assertAlmostEqual(stellar_model[0].X_H[-1], 0.63084 | units.none, 3)
            
            self.assertEqual(instance.number_of_particles, 2)
            merge_product = Particle()
            merge_product.primary = instance.particles[0]
            merge_product.secondary = instance.particles[1]
            instance.merge_products.add_particle(merge_product)
            self.assertEqual(instance.number_of_particles, 3)
            print instance.particles.number_of_zones
            #self.assertEqual(instance.particles.number_of_zones, [1043, 985, 935])
            
            mmams_merged_model = instance.merge_products[0].internal_structure()
            print mmams_merged_model.mass[[0, -1]]
            print mmams_merged_model.radius[[0, -1]]
            print mmams_merged_model.temperature[[0, -1]]
            print mmams_merged_model.X_H[[0, -1]]
#~            self.assertAlmostEqual(mmams_merged_model.mass[[0, -1]],        [0.0259, 25.7055] | units.MSun, 3)
#~            self.assertAlmostEqual(mmams_merged_model.radius[[0, -1]],      [0.036,  7.574] | units.RSun, 2)
#~            self.assertAlmostEqual(mmams_merged_model.temperature[[0, -1]], [40844102.7, 1832217.8] | units.K, 1)
#~            self.assertAlmostEqual(mmams_merged_model.X_H[[0, -1]],         [0.6743, 0.6999] | units.none, 2)
            stellar_model = dict(
                mass        = mmams_merged_model.mass,
                radius      = mmams_merged_model.radius,
                rho         = mmams_merged_model.rho,
                temperature = mmams_merged_model.temperature,
                luminosity  = mmams_merged_model.luminosity,
                X_H  = mmams_merged_model.X_H,
                X_He = mmams_merged_model.X_He,
                X_C  = mmams_merged_model.X_C,
                X_N  = mmams_merged_model.X_N,
                X_O  = mmams_merged_model.X_O,
                X_Ne = mmams_merged_model.X_Ne,
                X_Mg = mmams_merged_model.X_Mg,
                X_Si = mmams_merged_model.X_Si,
                X_Fe = mmams_merged_model.X_Fe
            )
            instance.stop()
            with open(os.path.join(get_path_to_results(), "test_mmams_9.pkl"), 'w') as out_file:
                pickle.dump(stellar_model, out_file)
        
        merged_in_code = stellar_evolution.new_particle_from_model(stellar_model, 0.0 | units.Myr)
        stellar_evolution.particles.remove_particles(stars)
        print 'number_of_backups_in_a_row', stellar_evolution.particles.get_number_of_backups_in_a_row()
        print stellar_evolution.particles
        for i in range(100):
            stellar_evolution.evolve_model((1+i) | units.yr)
            print stellar_evolution.particles
        stellar_evolution.stop()
    
    def slowtest10(self):
        print "Test 10: MakeMeAMassiveStar with EVtwin particles - import product into EVtwin (WIP)"
        stars = Particles(2)
        stars.mass = [20.0, 8.0] | units.MSun
        
        instance = MakeMeAMassiveStar(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        self.assertEqual(instance.number_of_particles, 2)
        
        stellar_evolution = EVtwin(redirection="none")
        stellar_evolution.initialize_code() 
        stellar_evolution.commit_parameters()
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(2 | units.Myr)
        
        for i in [0, 1]:
            stellar_model = stellar_evolution.particles[i].internal_structure()
            instance.particles[i].add_shell(
                stellar_model.d_mass, stellar_model.mass, stellar_model.radius, 
                stellar_model.rho, stellar_model.pressure, stellar_model.temperature, 
                stellar_model.luminosity, stellar_model.molecular_weight, stellar_model.X_H, 
                stellar_model.X_He, stellar_model.X_C, stellar_model.X_N, stellar_model.X_O, 
                stellar_model.X_Ne, stellar_model.X_Mg, stellar_model.X_Si, stellar_model.X_Fe
            )
            self.assertEqual(instance.particles[i].number_of_zones, len(stellar_model))
        
        self.assertEqual(instance.number_of_particles, 2)
        merge_product = Particle()
        merge_product.primary = instance.particles[0]
        merge_product.secondary = instance.particles[1]
        instance.merge_products.add_particle(merge_product)
        self.assertEqual(instance.number_of_particles, 3)
        self.assertEqual(instance.particles.number_of_zones, [199, 199, 16894])
        
        stellar_model = instance.merge_products[0].internal_structure()
        self.assertAlmostEqual(stellar_model.mass[[0, -1]],        [0.0, 25.6378] | units.MSun, 3)
        self.assertAlmostEqual(stellar_model.radius[[0, -1]],      [0.0,  8.1857] | units.RSun, 3)
        self.assertAlmostEqual(stellar_model.temperature[[0, -1]], [41250979.6, 445997.1] | units.K, 0)
        self.assertAlmostEqual(stellar_model.X_H[[0, -1]],         [0.68408, 0.70005] | units.none, 4)
        
        stellar_evolution.new_particle_from_model(stellar_model, 10.0 | units.Myr)
        instance.stop()
        print stellar_evolution.particles
        for i in range(10):
            stellar_evolution.evolve_model(keep_synchronous = False)
            print stellar_evolution.particles
        stellar_evolution.stop()
    
    def slowtest11(self):
        print "Test 11: MakeMeAMassiveStar with MESA particles of various masses/ages"
        masses = [1.0, 5.0, 20.0, 42.0, 80.0, 200.0] | units.MSun
        number_of_stars = len(masses)
        stars = Particles(number_of_stars)
        stars.mass = masses
        
        stellar_evolution = self.new_instance(MESA, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.parameters.metallicity = 0.0
        mesa_particles = stellar_evolution.particles.add_particles(stars)
        stellar_evolution.evolve_model(1.5 | units.Myr)
        
        mesa_models = []
        for i in range(number_of_stars):
            number_of_zones     = mesa_particles[i].get_number_of_zones()
            mass_profile        = mesa_particles[i].get_mass_profile(
                number_of_zones = number_of_zones) * mesa_particles[i].mass
            cumul_mass_profile  = mesa_particles[i].get_cumulative_mass_profile(
                number_of_zones = number_of_zones) * mesa_particles[i].mass
            density_profile     = mesa_particles[i].get_density_profile(number_of_zones = number_of_zones)
            radius_profile      = mesa_particles[i].get_radius_profile(number_of_zones = number_of_zones)
            temperature_profile = mesa_particles[i].get_temperature_profile(number_of_zones = number_of_zones)
            pressure_profile    = mesa_particles[i].get_pressure_profile(number_of_zones = number_of_zones)
            luminosity_profile  = mesa_particles[i].get_luminosity_profile(number_of_zones = number_of_zones)
            mu_profile          = mesa_particles[i].get_mu_profile(number_of_zones = number_of_zones)
            composition_profile = mesa_particles[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
            
            mesa_models.append((mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                composition_profile[1]+composition_profile[2], composition_profile[3], 
                composition_profile[4], composition_profile[5], composition_profile[6], 
                composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0))
        
        stellar_evolution.stop()
        
        instance = MakeMeAMassiveStar(**default_options)#debugger = "gdb", **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        stellar_models = []
        crashed = False
        for (index_1, index_2) in itertools.combinations(range(number_of_stars), 2):
            print
            print masses[index_1], masses[index_2]
            print
            self.assertEqual(instance.number_of_particles, len(stellar_models))
            
            colliding = instance.particles.add_particles(stars[[index_1, index_2]])
            colliding[0].add_shell(*(mesa_models[index_1]))
            colliding[1].add_shell(*(mesa_models[index_2]))
            self.assertEqual(instance.number_of_particles, 2 + len(stellar_models))
            
            merge_product = Particle()
            merge_product.primary = colliding[0]
            merge_product.secondary = colliding[1]
            try:
                instance.merge_products.add_particle(merge_product)
                self.assertEqual(instance.number_of_particles, 3 + len(stellar_models))
                instance.particles.remove_particles(stars[[index_1, index_2]])
                self.assertEqual(instance.number_of_particles, 1 + len(stellar_models))
                stellar_models.append(instance.particles[len(stellar_models)].internal_structure().copy_to_memory())
            except CodeException as ex:
                print ex
                crashed = True
                break
            
        if not crashed:
            instance.stop()
        
        for stellar_model in stellar_models:
            print "mass:", stellar_model.mass[-1]
            print "radius:", stellar_model.radius[-1]
            print "X_H:", stellar_model.X_H[[0, -1]]
        
        print "Storing the merger products in", os.path.join(get_path_to_results(), "test_mmams_11.pkl")
        mmams_merged_models = []
        for mmams_merged_model in stellar_models:
            mmams_merged_models.append(dict(
                mass        = mmams_merged_model.mass,
                radius      = mmams_merged_model.radius,
                rho         = mmams_merged_model.rho,
                temperature = mmams_merged_model.temperature,
                luminosity  = mmams_merged_model.luminosity,
                X_H  = mmams_merged_model.X_H,
                X_He = mmams_merged_model.X_He,
                X_C  = mmams_merged_model.X_C,
                X_N  = mmams_merged_model.X_N,
                X_O  = mmams_merged_model.X_O,
                X_Ne = mmams_merged_model.X_Ne,
                X_Mg = mmams_merged_model.X_Mg,
                X_Si = mmams_merged_model.X_Si,
                X_Fe = mmams_merged_model.X_Fe
            ))
        with open(os.path.join(get_path_to_results(), "test_mmams_11.pkl"), 'w') as out_file:
            pickle.dump(mmams_merged_models, out_file)
    
    def slowtest11b(self):
        print "Test 11b: Continue the stellar evolution of the products from test 11 with MESA"
        if os.path.exists(os.path.join(get_path_to_results(), "test_mmams_11.pkl")):
            with open(os.path.join(get_path_to_results(), "test_mmams_11.pkl"), 'r') as in_file:
                stellar_models = pickle.load(in_file)
        else:
            return
        print len(stellar_models)
        
        stellar_evolution = self.new_instance(MESA)#, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.initialize_code() 
        stellar_evolution.parameters.metallicity = 0.0
        stellar_evolution.commit_parameters()
        for i, stellar_model in enumerate(stellar_models):
            if i != 4:
                stellar_evolution.model_time = 0 | units.yr
                print "\n\n\n***", i, "***", stellar_model['mass'][-1]
                stellar_evolution.particles.remove_particles(stellar_evolution.particles)
                merged_in_code = stellar_evolution.new_particle_from_model(stellar_model, 0.0 | units.Myr)
                print stellar_evolution.particles
                for i in range(10):
                    stellar_evolution.evolve_model((1+i) | units.yr)
                    print stellar_evolution.particles
        
        stellar_evolution.stop()
    
    def slowtest12(self):
        print "Test 12: same as slowtest11, but new mmams instance each merger"
        masses = [1.0, 5.0, 20.0, 42.0, 80.0, 200.0] | units.MSun
#        masses = [5.0, 1000.0] | units.MSun#, 20.0, 42.0, 80.0, 200.0] | units.MSun
        number_of_stars = len(masses)
        stars = Particles(number_of_stars)
        stars.mass = masses
        
        stellar_evolution = self.new_instance(MESA, redirection="none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stellar_evolution.parameters.metallicity = 0.0
        mesa_particles = stellar_evolution.particles.add_particles(stars)
        stellar_evolution.evolve_model(1.0 | units.Myr)
        
        big_brother = Particle()
        big_brother.mass = 1000.0 | units.MSun
        stars.add_particle(big_brother)
        mesa_big_brother = stellar_evolution.particles.add_particle(big_brother)
        mesa_particles = mesa_particles + mesa_big_brother
        number_of_stars += 1
        stellar_evolution.evolve_model(1.5 | units.Myr)
        
        mesa_models = []
        for i in range(number_of_stars):
            number_of_zones     = mesa_particles[i].get_number_of_zones()
            mass_profile        = mesa_particles[i].get_mass_profile(
                number_of_zones = number_of_zones) * mesa_particles[i].mass
            cumul_mass_profile  = mesa_particles[i].get_cumulative_mass_profile(
                number_of_zones = number_of_zones) * mesa_particles[i].mass
            density_profile     = mesa_particles[i].get_density_profile(number_of_zones = number_of_zones)
            radius_profile      = mesa_particles[i].get_radius_profile(number_of_zones = number_of_zones)
            temperature_profile = mesa_particles[i].get_temperature_profile(number_of_zones = number_of_zones)
            pressure_profile    = mesa_particles[i].get_pressure_profile(number_of_zones = number_of_zones)
            luminosity_profile  = mesa_particles[i].get_luminosity_profile(number_of_zones = number_of_zones)
            mu_profile          = mesa_particles[i].get_mu_profile(number_of_zones = number_of_zones)
            composition_profile = mesa_particles[i].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
            
            mesa_models.append((mass_profile, cumul_mass_profile, radius_profile, density_profile, 
                pressure_profile, temperature_profile, luminosity_profile, mu_profile, composition_profile[0], 
                composition_profile[1]+composition_profile[2], composition_profile[3], 
                composition_profile[4], composition_profile[5], composition_profile[6], 
                composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0))
        
        stellar_evolution.stop()
        
        stellar_models = []
        for (index_1, index_2) in itertools.combinations(range(number_of_stars), 2):
        #for (index_1, index_2) in [(1, 6)]:
        #for (index_1, index_2) in [(0, 1)]:
            print
            print stars[index_1].mass, stars[index_2].mass
            print
            instance = MakeMeAMassiveStar(**default_options)
            instance.initialize_code()
            instance.commit_parameters()
            self.assertEqual(instance.number_of_particles, 0)
            
            colliding = instance.particles.add_particles(stars[[index_1, index_2]])
            colliding[0].add_shell(*(mesa_models[index_1]))
            colliding[1].add_shell(*(mesa_models[index_2]))
            self.assertEqual(instance.number_of_particles, 2)
            
            merge_product = Particle()
            merge_product.primary = colliding[0]
            merge_product.secondary = colliding[1]
            instance.merge_products.add_particle(merge_product)
            self.assertEqual(instance.number_of_particles, 3)
            instance.particles.remove_particles(stars[[index_1, index_2]])
            self.assertEqual(instance.number_of_particles, 1)
            stellar_models.append(instance.particles[0].internal_structure().copy_to_memory())
            instance.stop()
        
        for stellar_model in stellar_models:
            print "mass:", stellar_model.mass[-1]
            print "radius:", stellar_model.radius[-1]
            print "X_H:", stellar_model.X_H[[0, -1]]
        
    

