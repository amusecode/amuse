from amuse.test.amusetest import TestWithMPI
import sys
import os.path
from numpy import pi, arange

from amuse.community.evtwin.interface import EVtwin, EVtwinInterface

from amuse.support.exceptions import AmuseException
from amuse.units import nbody_system
from amuse.units import units
from amuse.units.quantities import new_quantity
from amuse import datamodel
from amuse.rfi import channel


class TestInterface(TestWithMPI):
    
    def test1(self):
        print "Testing get/set for metallicity..."
        instance = EVtwinInterface()
        (metallicity, error) = instance.get_metallicity()
        self.assertEquals(0, error)
        self.assertEquals(0.02, metallicity)
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)
        error = instance.initialize_code()
        self.assertEquals(0, error)
        error = instance.commit_parameters()
        self.assertEquals(0, error)
        instance.stop()

        for x in [0.03, 0.01, 0.004, 0.001, 0.0003, 0.0001]:
            instance = EVtwinInterface()
            error = instance.set_metallicity(x)
            self.assertEquals(0, error)
            (metallicity, error) = instance.get_metallicity()
            self.assertEquals(0, error)
            self.assertEquals(x, metallicity)
            error = instance.set_ev_path(instance.get_data_directory())
            self.assertEquals(0, error)      
            error = instance.initialize_code()
            error = instance.commit_parameters()
            if os.path.exists(os.path.join(instance.get_data_directory(),'zams','zams'+str(x)[2:]+'.mod')):
                self.assertEquals(0, error)
            else:
                self.assertEquals(-1, error)
                print "Metallicities other than solar need additional starting models!"
            instance.stop()
    
    def test2(self):
        print "Testing get/set for maximum number of stars..."
        instance = EVtwinInterface()
        
        (value, error) = instance.get_maximum_number_of_stars()
        self.assertEquals(0, error)      
        
        for x in range(10,100,10):
            error = instance.set_maximum_number_of_stars(x)
            self.assertEquals(0, error)      
            
            (value, error) = instance.get_maximum_number_of_stars()
            self.assertEquals(0, error)      
            self.assertEquals(x, value)      
        
        instance.stop()
    
    def test4(self):
        print "Testing initialization..."
        instance = EVtwinInterface()
        
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        instance.stop()
    
    def test5(self):
        print "Testing basic operations (new_particle, evolve_one_step etc.)..."
        #code/library_v2.f:602
        instance = EVtwinInterface()
        
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)
        
        error = instance.initialize_code()
        self.assertEquals(0, error)
        
        error = instance.commit_parameters()
        self.assertEquals(0, error)
        
        (index_of_the_star, error) = instance.new_particle(1.05)
        self.assertEquals(0, error)
        
        self.assertTrue(index_of_the_star >= 0)
        
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)
        self.assertEquals(1.05, mass)
        
        error = instance.evolve_one_step(index_of_the_star)
        self.assertEquals(0, error)
          
        for i in range(2):
            error = instance.evolve_one_step(index_of_the_star)
            self.assertEquals(0, error)
    
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)
        self.assertTrue(mass < 1.05)
    
        (age, error) = instance.get_age(index_of_the_star)
        self.assertEquals(0, error)
        self.assertTrue(age > 0)
        
        (x, error) = instance.get_mass_transfer_rate(index_of_the_star)
        self.assertEquals(0, error)
        self.assertEquals(0, x)
        
        (x, error) = instance.get_wind_mass_loss_rate(index_of_the_star)
        self.assertEquals(0, error)
        self.assertTrue(x < 1e-13, "mass loss should be less than 1e-13, it is {0}".format(x))
        self.assertTrue(x > 0, "mass loss rate should be more than 0, it is {0}".format(x))
        
        instance.stop()
        
    def test6(self):
        print "Testing EVtwin stop conditions..."
        instance = EVtwinInterface()
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)      
        error = instance.initialize_code()
        self.assertEquals(0, error)      
                
        (value, error) = instance.get_max_age_stop_condition()
        self.assertEquals(0, error)
        self.assertEquals(2.0e12, value)
        for x in range(10,14):
            error = instance.set_max_age_stop_condition(10 ** x)
            self.assertEquals(0, error)
            (value, error) = instance.get_max_age_stop_condition()
            self.assertEquals(0, error)      
            self.assertEquals(10 ** x, value)
            
        (value, error) = instance.get_min_timestep_stop_condition()
        self.assertEquals(0, error)
        self.assertAlmostEqual(0.031689, value, 5)
        for x in range(-3,2):
            error = instance.set_min_timestep_stop_condition(10 ** x)
            self.assertEquals(0, error)
            (value, error) = instance.get_min_timestep_stop_condition()
            self.assertEquals(0, error)
            self.assertEquals(10 ** x, value)
        
        instance.stop()

    def test7(self):
        print "Testing EVtwin parameters..."
        instance = EVtwinInterface()
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)      
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        (value, error) = instance.get_number_of_ionization_elements()
        self.assertEquals(0, error)
        self.assertEquals(2, value)
        for x in range(1,10):
            error = instance.set_number_of_ionization_elements(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_number_of_ionization_elements()
            self.assertEquals(0, error)
            self.assertEquals(x, value)

        (value, error) = instance.get_convective_overshoot_parameter()
        self.assertEquals(0, error)
        self.assertEquals(0.12, value)
        for x in [0.0, 0.1, 0.12, 0.15]:
            error = instance.set_convective_overshoot_parameter(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_convective_overshoot_parameter()
            self.assertEquals(0, error)
            self.assertEquals(x, value)

        (value, error) = instance.get_mixing_length_ratio()
        self.assertEquals(0, error)
        self.assertEquals(2.0, value)
        for x in [0.1, 1.0, 3.0]:
            error = instance.set_mixing_length_ratio(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_mixing_length_ratio()
            self.assertEquals(0, error)
            self.assertEquals(x, value)

        (value, error) = instance.get_semi_convection_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.04, value)
        for x in [0.0, 0.1]:
            error = instance.set_semi_convection_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_semi_convection_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_thermohaline_mixing_parameter()
        self.assertEquals(0, error)
        self.assertEquals(1.0, value)
        for x in [0.0, 0.5, 1.5]:
            error = instance.set_thermohaline_mixing_parameter(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_thermohaline_mixing_parameter()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_AGB_wind_setting()
        self.assertEquals(0, error)
        self.assertEquals(1, value)
        error = instance.set_AGB_wind_setting(2)
        self.assertEquals(0, error)
        (value, error) = instance.get_AGB_wind_setting()
        self.assertEquals(0, error)
        self.assertEquals(2, value)
        
        (value, error) = instance.get_RGB_wind_setting()
        self.assertEquals(0, error)
        self.assertEquals(1.0, value)
        for x in [-1.0, -0.5, 0.0, 1.0]:
            error = instance.set_RGB_wind_setting(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_RGB_wind_setting()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        instance.stop()

    def test8(self):
        print "Testing basic operations for spinning particle (new_spinning_particle, get_spin etc.)..."
        instance = EVtwinInterface()
        error = instance.set_ev_path(instance.get_data_directory())
        self.assertEquals(0, error)
        error = instance.initialize_code()
        self.assertEquals(0, error)
        error = instance.commit_parameters()
        self.assertEquals(0, error)

        
        (index_of_default_star, error) = instance.new_particle(1.05)
        self.assertEquals(0, error)
        self.assertTrue(index_of_default_star >= 0)
        (index_of_spinning_star, error) = instance.new_spinning_particle(1.05, 300.0)
        self.assertEquals(0, error)
        self.assertTrue(index_of_spinning_star >= 0)
        
        (spin, error) = instance.get_spin(index_of_default_star)
        self.assertEquals(0, error)
        self.assertAlmostEqual(54450.2652, spin,3)
        (spin, error) = instance.get_spin(index_of_spinning_star)
        self.assertEquals(0, error)
        self.assertEquals(300.0, spin)
        
        instance.stop()
    
        
    def test9(self):
        print "Testing adding and removing particles from stellar evolution code..."
        instance = EVtwinInterface()
        self.assertEquals(0, instance.set_ev_path(instance.get_data_directory()))
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        
        (indices, errors) = instance.new_particle([1.0, 1.0])
        self.assertEquals(errors, [0, 0])
        self.assertEquals(indices, [1, 2])
        
        self.assertEquals(0, instance.commit_particles())
        

        for index in indices:
            self.assertEquals(0, instance.evolve_one_step(index))
            (age_after_evolve, error) = instance.get_age(index)
            self.assertEquals(0, error)
            self.assertAlmostEqual(age_after_evolve, 422790.633054, 5)
        
        self.assertEquals(0, instance.delete_star(1))
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        (age, error) = instance.get_age(1)
        self.assertEquals(error, -1)
        
        (indices, errors) = instance.new_particle([1.0, 1.0])
        self.assertEquals(errors, [0, 0])
        self.assertEquals(indices, [1, 3])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 3)
        
        for index in indices:
            (age, error) = instance.get_age(index)
            self.assertEquals(0, error)
            self.assertTrue(age < age_after_evolve)
            self.assertEquals(0, instance.evolve_one_step(index))
            (age, error) = instance.get_age(index)
            self.assertEquals(0, error)
            self.assertAlmostEqual(age, age_after_evolve)
        
        instance.stop()
    

class TestEVtwin(TestWithMPI):
    
            
    
    def test1(self):
        print "Testing assigned default parameter values..."
        instance = EVtwin()
        
        instance.parameters.set_defaults()
        
        self.assertEquals(10.0 | units.no_unit, instance.parameters.maximum_number_of_stars)
        instance.parameters.maximum_number_of_stars = 12 | units.no_unit
        self.assertEquals(12.0 | units.no_unit, instance.parameters.maximum_number_of_stars)
        instance.stop()
    
    def test2(self):
        print "Testing basic operations (setup_particles, initialize_stars etc.)..."
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        
        path = os.path.join(instance.default_path_to_ev_database, 'run')
        path = os.path.join(path, 'muse')
        
        #instance.set_init_dat_name(os.path.join(path,'init.dat'))
        #instance.set_init_run_name(os.path.join(path,'init.run'))
        
        stars = datamodel.Stars(1)
        stars[0].mass = 10 | units.MSun
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        self.assertEquals(stars[0].mass, 10 | units.MSun)
        self.assertAlmostEquals(stars[0].luminosity.value_in(units.LSun), 5695.19757302, 6)
        instance.stop()
    
    def xtest3(self):
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        stars =  datamodel.Stars(1)
        
        star = stars[0]
        star.mass = 5.0 | units.MSun
        star.radius = 0.0 | units.RSun
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        previous_type = star.stellar_type
        results = []
        t0 = 0 | units.Myr
        current_time = t0
        
        while current_time < (115 | units.Myr):
            instance.evolve_model()
            from_code_to_model.copy()
            
            current_time = star.age
            print (star.age, star.mass, star.stellar_type)
            if not star.stellar_type == previous_type:
                results.append((star.age, star.mass, star.stellar_type))
                previous_type = star.stellar_type
        
        print results
        self.assertEqual(len(results), 6)
        
        times = ( 
            104.0 | units.Myr, 
            104.4 | units.Myr, 
            104.7 | units.Myr, 
            120.1 | units.Myr,
            120.9 | units.Myr,
            121.5 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 1)
            
        masses = ( 
            5.000 | units.MSun, 
            5.000 | units.MSun, 
            4.998 | units.MSun, 
            4.932 | units.MSun,
            4.895 | units.MSun,
            0.997 | units.MSun
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 3)
         
        types = (
            "Hertzsprung Gap",
            "First Giant Branch",
            "Core Helium Burning",
            "First Asymptotic Giant Branch",
            "Second Asymptotic Giant Branch",
            "Carbon/Oxygen White Dwarf",
        )
        
        for result, expected in zip(results, types):
            self.assertEquals(str(result[2]), expected)
        
        instance.stop()
        
    def test4(self):
        print "Testing max age stop condition..."
        #masses = [.5, 1.0, 1.5] | units.MSun # Test with fewer particles for speed-up.
        masses = [.5] | units.MSun
        max_age = 9.0 | units.Myr

        number_of_stars=len(masses)
        stars =  datamodel.Stars(number_of_stars)
        for i, star in enumerate(stars):
            star.mass = masses[i]
            star.radius = 0.0 | units.RSun

#       Initialize stellar evolution code
        instance = EVtwin() #debugger="xterm")
        instance.initialize_code()
        if instance.get_maximum_number_of_stars() < number_of_stars:
            instance.set_maximum_number_of_stars(number_of_stars)
        self.assertEqual(instance.parameters.max_age_stop_condition, 2e6 | units.Myr)
        instance.parameters.max_age_stop_condition = max_age
        self.assertEqual(instance.parameters.max_age_stop_condition, max_age)
        instance.commit_parameters()
        instance.particles.add_particles(stars)
#       Let the code perform initialization actions after all particles have been created. 
        instance.commit_particles()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        instance.evolve_model(end_time = 8.0 | units.Myr)
        from_code_to_model.copy()
        
        for i in range(number_of_stars):
            print stars[i].age.as_quantity_in(units.Myr)
            self.assertTrue(stars[i].age.value_in(units.Myr) >= 8.0)
            self.assertTrue(stars[i].age <= max_age)
            self.assertTrue(stars[i].mass <= masses[i])
            self.assertTrue(stars[i].time_step <= max_age)
                
        self.assertRaises(AmuseException, instance.evolve_model, end_time = 2*max_age, 
            expected_message = "Error when calling 'evolve_for' of a 'EVtwin', errorcode "
                "is 5, error is 'PRINTB -- age greater than limit'")

        instance.stop()
        
    def test5(self):
        print "Testing adding and removing particles from stellar evolution code..."
        
        particles = datamodel.Particles(3)
        particles.mass = 1.0 | units.MSun
        
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        stars = instance.particles
        self.assertEquals(len(stars), 0) # before creation
        stars.add_particles(particles[:-1])
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        self.assertEquals(len(stars), 2) # before remove
        self.assertAlmostEqual(stars.age, 1.0 | units.Myr)
        
        stars.remove_particle(particles[0])
        self.assertEquals(len(stars), 1)
        self.assertEquals(instance.get_number_of_particles(), 1)
        instance.evolve_model(2.0 | units.Myr)
        self.assertAlmostEqual(stars[0].age, 2.0 | units.Myr)
        
        stars.add_particles(particles[::2])
        self.assertEquals(len(stars), 3) # it's back...
        self.assertAlmostEqual(stars[0].age, 2.0 | units.Myr)
        self.assertAlmostEqual(stars[1].age, 0.0 | units.Myr)
        self.assertAlmostEqual(stars[2].age, 0.0 | units.Myr) # ... and rejuvenated.
        
        instance.evolve_model(3.0 | units.Myr) # The young stars keep their age offset from the old star
        self.assertAlmostEqual(stars.age, [3.0, 1.0, 1.0] | units.Myr)
        instance.evolve_model(4.0 | units.Myr)
        self.assertAlmostEqual(stars.age, [4.0, 2.0, 2.0] | units.Myr)
        instance.stop()

    def test6(self):
        print "Test for obtaining the stellar structure model"
        stars = datamodel.Particles(2)
        stars.mass = [1.0, 10.0] | units.MSun
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model()
        self.assertEquals(instance.particles.get_number_of_zones(), [199, 199])
        self.assertEquals(len(instance.particles[0].get_radius_profile()), 199)
        self.assertRaises(AmuseException, instance.particles.get_radius_profile, 
            expected_message = "Querying radius profiles of more than one particle at a time is not supported.")
        self.assertEquals(len(instance.particles[1].get_density_profile()), 199)
        self.assertIsOfOrder(instance.particles[0].get_radius_profile()[-1],          1.0 | units.RSun)
        self.assertIsOfOrder(instance.particles[0].get_temperature_profile()[0],  1.0e7 | units.K)
        self.assertIsOfOrder(instance.particles[0].get_temperature_profile()[-1],  5.0e3 | units.K)
        radius1 = instance.particles[0].get_radius_profile()
        radius2 = radius1[:-1]
        radius2.prepend(0|units.m)
        delta_radius_cubed = (radius1**3 - radius2**3)
        total_mass = (4./3. * pi * instance.particles[0].get_density_profile() * delta_radius_cubed).sum()
        self.assertAlmostRelativeEqual(total_mass, stars[0].mass, places = 1)
        self.assertAlmostEquals(instance.particles[0].get_mu_profile(), [0.62]*199 | units.amu, places=1)
        instance.stop()
        del instance
    
    def test7(self):
        print "Test for obtaining the stellar composition structure"
        stars = datamodel.Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model()
        instance.evolve_model()
        number_of_zones   = instance.particles.get_number_of_zones()[0]
        number_of_species = instance.particles.get_number_of_species()[0]
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        species_names     = instance.particles[0].get_names_of_species()
        self.assertEquals(number_of_zones,    199)
        self.assertEquals(number_of_species,    9)
        self.assertEquals(len(species_names),  number_of_species)
        self.assertEquals(len(composition),    number_of_species)
        self.assertEquals(len(composition[0]), number_of_zones)
        self.assertEquals(species_names, ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 'si28', 'fe56'])
        self.assertAlmostEquals(composition[0, -1],        0.7, 4)
        self.assertAlmostEquals(composition[1, -1],        0.3 - instance.parameters.metallicity, 4)
        self.assertAlmostEquals(composition[2:,-1].sum(),  instance.parameters.metallicity, 4)
        self.assertAlmostEquals(composition.sum(axis=0), [1.0]*number_of_zones)
        instance.stop()
        del instance
    
    def slowtest8(self):
        print "Test for obtaining the stellar composition structure - evolved star"
        stars = datamodel.Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(11.7 | units.Gyr)
        self.assertTrue(instance.particles[0].age >= 11.7 | units.Gyr)
        self.assertTrue(str(instance.particles[0].stellar_type) == "First Giant Branch")
        number_of_zones   = instance.particles.get_number_of_zones()[0]
        number_of_species = instance.particles.get_number_of_species()[0]
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        species_names     = instance.particles[0].get_names_of_species()
        self.assertEquals(number_of_zones,    199)
        self.assertEquals(number_of_species,    9)
        self.assertEquals(len(species_names),  number_of_species)
        self.assertEquals(len(composition),    number_of_species)
        self.assertEquals(len(composition[0]), number_of_zones)
        self.assertEquals(species_names, ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 'si28', 'fe56'])
        self.assertAlmostRelativeEquals(composition[0, -1],        0.7 | units.none, 1)
        self.assertAlmostRelativeEquals(composition[1, -1],        0.3 - instance.parameters.metallicity, 1)
        self.assertAlmostRelativeEquals(composition[2:,-1].sum(),  instance.parameters.metallicity, 1)
        self.assertAlmostEquals(composition.sum(axis=0), [1.0]*number_of_zones | units.none)
        self.assertAlmostEquals(composition[0, 0],        0.00 | units.none)
        self.assertAlmostEquals(composition[1, 0],        1.00 - instance.parameters.metallicity, 3)
        self.assertAlmostEquals(composition[2:,0].sum(),  instance.parameters.metallicity, 3)
        instance.stop()
        del instance
    
    def xtest9(self):
        print "Test for changing the stellar structure model (not yet implemented)"
        star = datamodel.Particles(1)
        star.mass = 1.0 | units.MSun
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model()

        density_profile = instance.particles[0].get_density_profile()
        
        self.assertRaises(AmuseException, instance.particles[0].set_density_profile, density_profile[2:], 
            expected_message = "The length of the supplied vector (197) does not match the number of "
            "mesh zones of the star (199).")
        
        mass_factor = 1.1
        instance.particles[0].set_density_profile(mass_factor*density_profile)
        self.assertAlmostRelativeEqual(instance.particles[0].get_density_profile(), density_profile*mass_factor, places=10)
        instance.particles.mass *= mass_factor
        instance.evolve_model()
        
        outer_radius = instance.particles[0].get_radius_profile()
        inner_radius = outer_radius[:-1]
        inner_radius.prepend(0|units.m)
        delta_radius_cubed = (outer_radius**3 - inner_radius**3)
        integrated_mass = (4./3.*pi*delta_radius_cubed*instance.particles[0].get_density_profile()).sum()
        self.assertAlmostRelativeEqual(integrated_mass, star.mass*mass_factor, places = 3)
        instance.stop()
        del instance
    
    def xtest10(self):
        print "Test for changing the stellar composition (not yet implemented)"
        star = datamodel.Particles(1)
        star.mass = 1.0 | units.MSun
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model()
        
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        h1_profile = composition[0] * 1
        he4_profile = composition[1] * 1
        k_surface = -1 # index to the outer mesh cell (surface)
        
        self.assertAlmostEquals(composition[0, k_surface],  0.7 | units.none, 4)
        self.assertAlmostEquals(composition[1, k_surface],  (0.3 | units.none) - instance.parameters.metallicity, 4)
        self.assertAlmostEquals(composition[2: , k_surface].sum(),  instance.parameters.metallicity, 4)
        
        composition[0] = he4_profile
        composition[1] = h1_profile
        instance.particles[0].set_chemical_abundance_profiles(composition)
        instance.evolve_model()
        
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        self.assertAlmostEquals(composition[0, k_surface],  (0.3 | units.none) - instance.parameters.metallicity, 4)
        self.assertAlmostEquals(composition[1, k_surface],  0.7 | units.none, 4)
        self.assertAlmostEquals(composition[2: , k_surface].sum(),  instance.parameters.metallicity, 4)
        self.assertAlmostEquals(composition.sum(axis=0), 1.0 | units.none)
        
        self.assertRaises(AmuseException, instance.particles[0].set_chemical_abundance_profiles, composition[:7], 
            expected_message = "The length of the supplied vector (7) does not match the number of "
            "chemical species of the star (8).")
        instance.stop()
        del instance
    
    def slowtest11(self):
        print "Test 11: Continue the stellar evolution of a 'merger product' - WIP"
        instance = EVtwin()
        instance.initialize_code()
        instance.commit_parameters()
        
        instance.parameters.min_timestep_stop_condition = 1.0 | units.s
        
        stars = datamodel.Particles(3)
        stars.mass = [1.0, 2.0, 1.0] | units.MSun
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        stellar_models = instance.native_stars.internal_structure()
        
        self.assertEqual(len(stellar_models), 3)
        self.assertEqual(len(stellar_models[0]), 199)
        self.assertEqual(len(stellar_models[1]), 199)
        self.assertAlmostEqual(stellar_models[0].mass[198], 1.0 | units.MSun, 2)
        self.assertAlmostEqual(stellar_models[1].mass[198], 2.0 | units.MSun, 2)
        self.assertAlmostEqual(stellar_models[0].mass[0], 0.0 | units.MSun, 2)
        
        instance.new_particle_from_model(stellar_models[0], instance.particles[0].age)
        self.assertEqual(len(instance.particles), 4)
        self.assertEqual(len(instance.imported_stars), 1)
        imported_stellar_model = instance.imported_stars[0].internal_structure()
        self.assertEqual(len(imported_stellar_model), 199)
        self.assertAlmostEqual(imported_stellar_model.mass[198], 1.0 | units.MSun, 2)
        self.assertAlmostEqual(imported_stellar_model.mass[0], 0.0 | units.MSun, 2)
        self.assertAlmostRelativeEqual(imported_stellar_model.X_H, stellar_models[0].X_H, 5)
        self.assertAlmostRelativeEqual(imported_stellar_model.X_He, stellar_models[0].X_He, 5)
        self.assertAlmostRelativeEqual(imported_stellar_model.mass, stellar_models[0].mass, 2)
        self.assertAlmostRelativeEqual(imported_stellar_model.radius[1:], stellar_models[0].radius[1:], 2)
#        instance.evolve_model(2.0 | units.Myr)
        print instance.particles
        instance.stop()
        del instance
    
    def test12(self):
        print "Testing basic operations: evolve_one_step and evolve_for"
        stars = datamodel.Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = EVtwin()
        se_stars = instance.particles.add_particles(stars)
        self.assertAlmostEqual(se_stars.age, [0.0, 0.0] | units.yr)
        
        for i in range(3):
            se_stars[0].evolve_one_step()
        self.assertAlmostEqual(se_stars.age, [1352930.0257, 0.0] | units.yr, 3)
        number_of_steps = 10
        step_size = se_stars[0].age / number_of_steps
        for i in range(1, number_of_steps + 1):
            se_stars[1].evolve_for(step_size)
            self.assertAlmostEqual(se_stars.age, [number_of_steps, i] * step_size)
        print se_stars
        self.assertAlmostRelativeEqual(se_stars[0].age,         se_stars[1].age)
        self.assertAlmostRelativeEqual(se_stars[0].luminosity,  se_stars[1].luminosity, 3)
        self.assertAlmostRelativeEqual(se_stars[0].radius,      se_stars[1].radius, 3)
        self.assertAlmostRelativeEqual(se_stars[0].temperature, se_stars[1].temperature, 3)
        instance.stop()
    
    def test13(self):
        print "Test evolve_model optional arguments: end_time and keep_synchronous"
        stars = datamodel.Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        instance = EVtwin()
        instance.particles.add_particles(stars)
        
        self.assertAlmostEqual(instance.particles.age, [0.0, 0.0, 0.0] | units.yr)
        self.assertAlmostEqual(instance.particles.time_step, [422790.6330, 36382.1271, 11259.1953] | units.yr, 3)
        
        print "evolve_model without arguments: use shared timestep = 0.99*min(particles.time_step)"
        instance.evolve_model()
        self.assertAlmostEqual(instance.particles.age, 0.99*([11259.1953, 11259.1953, 11259.1953] | units.yr), 3)
        self.assertAlmostEqual(instance.particles.time_step, [422790.6330, 36382.1271, 11259.1953] | units.yr, 3)
        self.assertAlmostEqual(instance.model_time, 0.99*11259.1953 | units.yr, 3)
        
        print "evolve_model with end_time: take timesteps, until end_time is reached exactly"
        instance.evolve_model(15000 | units.yr)
        self.assertAlmostEqual(instance.particles.age, [15000.0, 15000.0, 15000.0] | units.yr, 3)
        self.assertAlmostEqual(instance.particles.time_step, [422790.6330, 36382.1271, 11259.1953] | units.yr, 3)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3)
        
        print "evolve_model with keep_synchronous: use non-shared timestep, particle ages will typically diverge"
        instance.evolve_model(keep_synchronous = False)
        self.assertAlmostEqual(instance.particles.age, (15000 | units.yr) + ([422790.6330, 36382.1271, 11259.1953] | units.yr), 3)
        self.assertAlmostEqual(instance.particles.time_step, [507348.7596, 43180.1460, 13511.0343] | units.yr, 3)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3) # Unchanged!
        instance.stop()
    
    def test14(self):
        print "Testing EVtwin states"
        stars = datamodel.Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = EVtwin()
        
        print "First do everything manually:",
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particle(stars[0])
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()
        print "ok"

        print "initialize_code(), commit_parameters(), (re)commit_particles(), " \
            "and cleanup_code() should be called automatically:",
        instance = EVtwin()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.RGB_wind_setting = -0.5
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particle(stars[0])
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.particles.add_particle(stars[1])
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
        print "ok"
