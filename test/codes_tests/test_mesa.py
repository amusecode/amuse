from amuse.test.amusetest import TestWithMPI, get_path_to_results
import sys
import os.path
from subprocess import PIPE, Popen
import numpy
from math import ceil

from amuse.community.mesa.interface import MESA, MESAInterface

from amuse.support.exceptions import AmuseException
from amuse.units import units
from amuse.datamodel import Particles
from amuse.datamodel import Particle
from amuse.rfi import channel
from amuse.ext.spherical_model import EnclosedMassInterpolator

class TestMESAInterface(TestWithMPI):
    
    def test1(self):
        print("Testing initialization of the interface...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        inlist_path = instance.default_path_to_inlist
        #print "Path to inlist: ", inlist_path
        MESA_data_path = instance.default_path_to_MESA_data
        #print "Path to MESA data directory: ", MESA_data_path
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        instance.stop()
        
    def slowtest2(self):
        print("Testing get/set of metallicity (tests new ZAMS model implicitly)...")
        print("The first time this test will take quite some time" \
            " to generate new starting models.")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        (metallicity, error) = instance.get_metallicity()
        self.assertEqual(0, error)
        self.assertEqual(0.02, metallicity)
        for x in [0.04, 0.02, 0.01, 0.00]:
            error = instance.set_metallicity(x)
            self.assertEqual(0, error)
            (metallicity, error) = instance.get_metallicity()
            self.assertEqual(0, error)
            self.assertEqual(x, metallicity)
        instance.stop()
    
    def test3(self):
        print("Testing basic operations: new_particle...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        (maximum_number_of_stars, error) = instance.get_maximum_number_of_stars()
        self.assertEqual(0, error)
        self.assertEqual(maximum_number_of_stars,1000)
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        number_of_stars = 1
        for i in range(number_of_stars):
            #print i
            (index_of_the_star, error) = instance.new_particle(0.5+i*1.0/number_of_stars)
            self.assertEqual(0, error)
            self.assertEqual(index_of_the_star,i+1)
        #import time    # To check the amount of memory is used ...
        #time.sleep(10) # when a large number of stars is created.
        instance.stop()
        
    def test4(self):
        print("Testing basic operations: evolve...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        status = instance.initialize_code()
        (index_of_the_star, error) = instance.new_particle(1.0)
        self.assertEqual(0, error)
        self.assertEqual(index_of_the_star, 1)
        self.assertEqual(0, instance.commit_particles())
        
        initial_dt = 1.0e5
        dt_factor = 1.2
        self.assertEqual([initial_dt, 0], list(instance.get_time_step(index_of_the_star).values()))
        self.assertEqual(0, instance.evolve_one_step(index_of_the_star))
        self.assertEqual([initial_dt, 0], list(instance.get_age(index_of_the_star).values()))
        
        target_end_time = 3.0e5 # (years)
        self.assertEqual(0, instance.evolve_for(index_of_the_star, target_end_time-initial_dt))
        self.assertEqual([initial_dt*(1 + dt_factor + dt_factor**2), 0], list(instance.get_age(index_of_the_star).values()))
        self.assertEqual([round(initial_dt*dt_factor**3), 0], list(instance.get_time_step(index_of_the_star).values()))
        self.assertTrue(instance.get_age(index_of_the_star)['age'] >= target_end_time)
        
        (L_of_the_star, error) = instance.get_luminosity(index_of_the_star)
        self.assertEqual(0, error)
        self.assertAlmostEqual(L_of_the_star,0.725,1)
        (M_of_the_star, error) = instance.get_mass(index_of_the_star)
        self.assertEqual(0, error)
        self.assertAlmostEqual(M_of_the_star,1.000,3)
        (T_of_the_star, error) = instance.get_temperature(index_of_the_star)
        self.assertEqual(0, error)
        self.assertAlmostEqual(T_of_the_star,5650.998,-2)
        instance.stop()
    
    def slowtest5(self):
        print("Testing evolve with varying Z (tests new ZAMS model implicitly)...")
        print("If the required starting models do not exist, this test will " \
            "take quite some time to generate them.")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        metallicities = [0.00, 0.01, 0.02, 0.04]
        luminosities = [1.717, 0.938, 0.725, 0.592]
        for (i,(Z,L)) in enumerate(zip(metallicities, luminosities)):
            status = instance.set_metallicity(Z)
            self.assertEqual(0, status)
            (index_of_the_star, status) = instance.new_particle(1.0)
            self.assertEqual(0, status)
            self.assertEqual(index_of_the_star,i+1)
            instance.evolve_for(index_of_the_star, 5.0e5)
            (L_of_the_star, status) = instance.get_luminosity(index_of_the_star)
            self.assertEqual(0, status)
            self.assertAlmostEqual(L_of_the_star,L,1)
        instance.stop()
    
    def test6(self):
        print("Testing MESA stop conditions...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        (value, error) = instance.get_max_age_stop_condition()
        self.assertEqual(0, error) 
        self.assertEqual(1.0e36, value)
        for x in range(10,14):
            error = instance.set_max_age_stop_condition(10 ** x)
            self.assertEqual(0, error)
            (value, error) = instance.get_max_age_stop_condition()
            self.assertEqual(0, error)
            self.assertEqual(10 ** x, value)
        (value, error) = instance.get_min_timestep_stop_condition()
        self.assertEqual(0, error)
        self.assertEqual(1.0e-6, value)
        for x in range(-9,-2):
            error = instance.set_min_timestep_stop_condition(10.0 ** x)
            self.assertEqual(0, error)
            (value, error) = instance.get_min_timestep_stop_condition()
            self.assertEqual(0, error)
            self.assertEqual(10.0 ** x, value)
        instance.stop()
    
    def test7(self):
        print("Testing MESA parameters...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        (value, error) = instance.get_mixing_length_ratio()
        self.assertEqual(0, error) 
        self.assertEqual(2.0, value)
        for x in [1.0, 1.5, 3.0]:
            error = instance.set_mixing_length_ratio(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_mixing_length_ratio()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        (value, error) = instance.get_semi_convection_efficiency()
        self.assertEqual(0, error)
        self.assertEqual(0.0, value)
        for x in [0.1, 0.04, 0.001]:
            error = instance.set_semi_convection_efficiency(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_semi_convection_efficiency()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        instance.stop()
    
    def test8(self):
        print("Testing MESA wind parameters...")
        instance = self.new_instance_of_an_optional_code(MESAInterface)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        (value, error) = instance.get_RGB_wind_scheme()
        self.assertEqual(0, error) 
        self.assertEqual(1, value)
        for x in range(6):
            error = instance.set_RGB_wind_scheme(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_RGB_wind_scheme()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        
        (value, error) = instance.get_AGB_wind_scheme()
        self.assertEqual(0, error) 
        self.assertEqual(1, value)
        for x in range(6):
            error = instance.set_AGB_wind_scheme(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_AGB_wind_scheme()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
            
        (value, error) = instance.get_reimers_wind_efficiency()
        self.assertEqual(0, error)
        self.assertEqual(0.5, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_reimers_wind_efficiency(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_reimers_wind_efficiency()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        
        (value, error) = instance.get_blocker_wind_efficiency()
        self.assertEqual(0, error)
        self.assertEqual(0.1, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_blocker_wind_efficiency(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_blocker_wind_efficiency()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        
        (value, error) = instance.get_de_jager_wind_efficiency()
        self.assertEqual(0, error)
        self.assertEqual(0.8, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_de_jager_wind_efficiency(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_de_jager_wind_efficiency()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        
        (value, error) = instance.get_dutch_wind_efficiency()
        self.assertEqual(0, error)
        self.assertEqual(0.8, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_dutch_wind_efficiency(x)
            self.assertEqual(0, error)
            (value, error) = instance.get_dutch_wind_efficiency()
            self.assertEqual(0, error)
            self.assertEqual(x, value)
        
        instance.stop()


class TestMESA(TestWithMPI):
    
    def test1(self):
        print("Testing initialization and default MESA parameters...")
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_output_directory())
        instance.initialize_code()
        self.assertEqual(0.02 | units.no_unit, instance.parameters.metallicity)
        self.assertEqual(1.0e36 | units.yr, instance.parameters.max_age_stop_condition)
        instance.parameters.max_age_stop_condition = 1.0e2 | units.Myr
        self.assertEqual(1.0e2 | units.Myr, instance.parameters.max_age_stop_condition)
        instance.stop()
    
    def test2(self):
        print("Testing basic operations: evolve and get_...")
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code() 
        instance.commit_parameters() 
        index_of_the_star = instance.new_particle(1.0 | units.MSun)
        self.assertEqual(index_of_the_star,1)
        
        initial_dt = 1.0e5 | units.yr
        dt_factor = 1.2
        time_step = instance.get_time_step(index_of_the_star)
        self.assertAlmostEqual(time_step, initial_dt)
        instance.evolve_one_step(index_of_the_star)
        age_of_the_star = instance.get_age(index_of_the_star)
        self.assertAlmostEqual(age_of_the_star, initial_dt)
        
        target_end_time = 3.0e5 | units.yr
        instance.evolve_for(index_of_the_star, target_end_time - age_of_the_star)
        self.assertAlmostEqual(initial_dt*(1 + dt_factor + dt_factor**2), instance.get_age(index_of_the_star))
        self.assertAlmostEqual(initial_dt*dt_factor**3, instance.get_time_step(index_of_the_star))
        self.assertTrue(instance.get_age(index_of_the_star) >= target_end_time)
        
        L_of_the_star = instance.get_luminosity(index_of_the_star)
        self.assertAlmostEqual(L_of_the_star,0.725 | units.LSun,1)
        M_of_the_star = instance.get_mass(index_of_the_star)
        self.assertAlmostEqual(M_of_the_star,1.000 | units.MSun,3)
        T_of_the_star = instance.get_temperature(index_of_the_star)
        self.assertAlmostEqual(T_of_the_star,5650.998 | units.K,-2)
        instance.stop()
    
    def test3(self):
        print("Testing basic operations: evolve_model and channels...")
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code() 
        instance.commit_parameters() 
        stars = Particles(1)
        mass = 10. | units.MSun
        stars.mass = mass
        instance.particles.add_particles(stars)
        instance.commit_particles()
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        #print stars
        #instance.evolve_model(end_time = 0.03 | units.Myr) # speeding test up:
        self.assertEqual(stars[0].mass, mass)
        self.assertAlmostRelativeEqual(stars[0].luminosity, 5841. | units.LSun, 1)
        instance.evolve_model()
        from_code_to_model.copy()
        self.assertAlmostEqual(stars[0].mass, mass, 5)
        self.assertAlmostRelativeEqual(stars[0].luminosity, 5820.85 | units.LSun, 1)
        instance.stop()
    
    def slowtest4(self):
        print("Testing stellar type...")
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code() 
        instance.commit_parameters() 
        stars =  Particles(1)
        
        star = stars[0]
        star.mass = 5.0 | units.MSun
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        previous_type = 15 | units.stellar_type
        results = []
        current_time = 0 | units.Myr
        
        while current_time < (101 | units.Myr):
            if not star.stellar_type == previous_type:
                print((star.age, star.mass, star.stellar_type))
                results.append((star.age, star.mass, star.stellar_type))
                previous_type = star.stellar_type
            instance.evolve_model()
            from_code_to_model.copy()
            current_time = star.age
        
        self.assertEqual(len(results), 4)
        
        times = ( 
            0.0 | units.Myr, 
            81.6 | units.Myr, 
            99.9 | units.Myr, 
            100.3 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 1)
            
        masses = ( 
            5.000 | units.MSun, 
            5.000 | units.MSun, 
            5.000 | units.MSun, 
            5.000 | units.MSun
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 3)
         
        types = (
            "Main Sequence star",
            "First Giant Branch",
            "Core Helium Burning",
            "First Asymptotic Giant Branch"
        )
        
        for result, expected in zip(results, types):
            self.assertEqual(str(result[2]), expected)
        
        instance.stop()
    
    def test5(self):
        print("Testing evolve_model for particle set...")
        instance = self.new_instance_of_an_optional_code(MESA)
        masses = [0.5, 1.0] | units.MSun
        max_age = 0.6 | units.Myr
        number_of_stars=len(masses)
        stars =  Particles(number_of_stars)
        stars.mass = masses
        instance.initialize_code()
        self.assertEqual(instance.parameters.max_age_stop_condition, 1e30 | units.Myr)
        instance.parameters.max_age_stop_condition = max_age
        self.assertEqual(instance.parameters.max_age_stop_condition, max_age)
        instance.particles.add_particles(stars)
        instance.commit_particles()
        from_code_to_model = instance.particles.new_channel_to(stars)
        instance.evolve_model(end_time = 0.5 | units.Myr)
        from_code_to_model.copy()
        #print stars
        for i in range(number_of_stars):
            self.assertTrue(stars[i].age.value_in(units.Myr) >= 0.5)
            self.assertTrue(stars[i].age <= max_age)
            self.assertTrue(stars[i].mass <= masses[i])
            self.assertTrue(stars[i].age+stars[i].time_step <= max_age)
        
        self.assertRaises(AmuseException, instance.evolve_model, end_time = 2*max_age, 
            expected_message = "Error when calling 'evolve_for' of a 'MESA', "
                "errorcode is -12, error is 'Evolve terminated: Maximum age reached.'")
        instance.stop()
    
    def test6(self):
        print("Test for obtaining the stellar structure model")
        stars = Particles(2)
        stars.mass = [1.0, 10.0] | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model()
        self.assertEqual(instance.particles.get_number_of_zones(), [575, 2262])
        self.assertEqual(len(instance.particles[0].get_mass_profile()), 575)
        self.assertAlmostEqual(instance.particles[0].get_mass_profile().sum(), 1.0)
        self.assertRaises(AmuseException, instance.particles.get_mass_profile, 
            expected_message = "Querying mass profiles of more than one particle at a time is not supported.")
        print(instance.particles)
        self.assertEqual(len(instance.particles[1].get_density_profile()), 2262)
        self.assertIsOfOrder(instance.particles[0].get_radius_profile()[-1],          1.0 | units.RSun)
        self.assertIsOfOrder(instance.particles[0].get_temperature_profile()[0],  1.0e7 | units.K)
        self.assertIsOfOrder(instance.particles[0].get_luminosity_profile()[-1],      1.0 | units.LSun)
        self.assertIsOfOrder(instance.particles[0].get_pressure_profile()[0],      1.0e17 | units.barye)
        delta_mass = instance.particles[0].get_mass_profile() * instance.particles[0].mass
        radius1 = instance.particles[0].get_radius_profile()
        radius2 = radius1[:-1]
        radius2.prepend(0|units.m)
        delta_radius_cubed = (radius1**3 - radius2**3)
        self.assertAlmostEqual(instance.particles[0].get_density_profile() / 
            (delta_mass/(4./3.*numpy.pi*delta_radius_cubed)), [1]*575, places=3)
        self.assertAlmostEqual(instance.particles[1].get_mu_profile(), [0.62]*2262 | units.amu, places=1)
        instance.stop()
    
    def test7(self):
        print("Test for obtaining the stellar composition structure")
        stars = Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model()
        number_of_zones   = instance.particles.get_number_of_zones()[0]
        number_of_species = instance.particles.get_number_of_species()[0]
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        species_names     = instance.particles[0].get_names_of_species()
        species_IDs       = instance.particles[0].get_IDs_of_species()
        species_masses    = instance.particles[0].get_masses_of_species()
        self.assertEqual(number_of_zones,    575)
        self.assertEqual(number_of_species,    8)
        self.assertEqual(len(species_names),  number_of_species)
        self.assertEqual(len(composition),    number_of_species)
        self.assertEqual(len(composition[0]), number_of_zones)
        self.assertEqual(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
        self.assertEqual(species_IDs,   [2,    5,     6,     38,    51,    69,    114,    168])
        self.assertAlmostEqual(species_masses, [1.0078250, 3.0160293, 4.0026032, 12.0, 
                                14.0030740, 15.9949146, 19.9924401, 23.9850417] | units.amu, places=5)
        self.assertAlmostEqual(composition[ :1,-1].sum(),  0.7)
        self.assertAlmostEqual(composition[1:3,-1].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEqual(composition[3: ,-1].sum(),  instance.parameters.metallicity)
        self.assertAlmostEqual(composition.sum(axis=0), [1.0]*number_of_zones)
        instance.stop()
    
    def slowtest8(self):
        print("Test for obtaining the stellar composition structure - evolved star with zero metalicity")
        stars = Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.parameters.metallicity = 0.0
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(5.85 | units.Gyr)
        self.assertTrue(instance.particles[0].age >= 5.85 | units.Gyr)
        self.assertTrue(str(instance.particles[0].stellar_type) == "First Giant Branch")
        number_of_zones   = instance.particles.get_number_of_zones()[0]
        number_of_species = instance.particles.get_number_of_species()[0]
        composition       = instance.particles[0].get_chemical_abundance_profiles()
        species_names     = instance.particles[0].get_names_of_species()
        species_IDs       = instance.particles[0].get_IDs_of_species()
        self.assertEqual(number_of_zones,    578)
        self.assertEqual(number_of_species,    8)
        self.assertEqual(len(species_names),  number_of_species)
        self.assertEqual(len(composition),    number_of_species)
        self.assertEqual(len(composition[0]), number_of_zones)
        self.assertEqual(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
        self.assertEqual(species_IDs,   [2,    5,     6,     38,    51,    69,    114,    168])
        self.assertAlmostEqual(composition[ :1,number_of_zones-1].sum(),  0.76 | units.none)
        self.assertAlmostEqual(composition[1:3,number_of_zones-1].sum(),  0.24 | units.none)
        self.assertAlmostEqual(composition[3: ,number_of_zones-1].sum(),  0.00 | units.none)
        self.assertAlmostEqual(composition.sum(axis=0), [1.0]*number_of_zones | units.none)
        self.assertAlmostEqual(composition[ :1,0].sum(),  0.00 | units.none)
        self.assertAlmostEqual(composition[1:3,0].sum(),  1.00 | units.none)
        self.assertAlmostEqual(composition[3: ,0].sum(),  0.00 | units.none)
        instance.stop()
    
    def test9(self):
        print("Test for changing the stellar structure model")
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model()

        density_profile = instance.particles[0].get_density_profile()
        
        self.assertRaises(AmuseException, instance.particles[0].set_density_profile, density_profile[2:], 
            expected_message = "The length of the supplied vector (573) does not match the number of "
            "mesh zones of the star (575).")
        
        mass_factor = 1.1
        instance.particles[0].set_density_profile(mass_factor*density_profile)
        self.assertAlmostRelativeEqual(instance.particles[0].get_density_profile(), density_profile*mass_factor, places=10)
        instance.particles.mass *= mass_factor
        instance.evolve_model()
        
        outer_radius = instance.particles[0].get_radius_profile()
        inner_radius = outer_radius[:-1]
        inner_radius.prepend(0|units.m)
        delta_radius_cubed = (outer_radius**3 - inner_radius**3)
        integrated_mass = (4./3.*numpy.pi*delta_radius_cubed*instance.particles[0].get_density_profile()).sum()
        self.assertAlmostRelativeEqual(integrated_mass, star.mass*mass_factor, places = 3)
        instance.stop()
    
    def test10(self):
        print("Test for changing the stellar composition")
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)#, redirection = 'none')
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model(0.3 | units.Myr)
        
        composition = instance.particles[0].get_chemical_abundance_profiles()
        k_surface = -1 # index to the outer mesh cell (surface)
        
        self.assertAlmostEqual(composition[ :1, k_surface].sum(),  0.7)
        self.assertAlmostEqual(composition[1:3, k_surface].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEqual(composition[3: , k_surface].sum(),  instance.parameters.metallicity)
        
        # Gradually and consistently increase helium and decrease hydrogen abundances until reversed
        for alpha in [0.3, 1.0, -0.5, -0.125]:
            h1_profile = composition[0] * 1
            he4_profile = composition[2] * 1
            composition[0] = alpha * he4_profile + (1-alpha) * h1_profile
            composition[2] = (1-alpha) * he4_profile + alpha * h1_profile
            instance.particles[0].set_chemical_abundance_profiles(composition)
            instance.evolve_model()
            instance.evolve_model()
            composition = instance.particles[0].get_chemical_abundance_profiles()
        
        self.assertAlmostEqual(composition[ :2, k_surface].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEqual(composition[2:3, k_surface].sum(),  0.7)
        self.assertAlmostEqual(composition[3: , k_surface].sum(),  instance.parameters.metallicity)
        self.assertAlmostEqual(composition.sum(axis=0), 1.0)
        
        self.assertRaises(AmuseException, instance.particles[0].set_chemical_abundance_profiles, composition[:7], 
            expected_message = "The length of the supplied vector (7) does not match the number of "
            "chemical species of the star (8).")
        instance.stop()
    
    def test11(self):
        print("Test evolve_model optional arguments: end_time and keep_synchronous")
        stars = Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        self.assertAlmostEqual(instance.particles.age, [0.0, 0.0, 0.0] | units.yr)
        self.assertAlmostEqual(instance.particles.time_step, [100000.0, 17677.6695, 6415.0029] | units.yr, 3)
        
        print("evolve_model without arguments: use shared timestep = 0.99*min(particles.time_step)")
        instance.evolve_model()
        self.assertAlmostEqual(instance.particles.age, [100000.0, 17677.6695, 6415.0029] | units.yr, 3)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.2*([100000.0, 17677.6695, 6415.0029] | units.yr), 6)
        self.assertAlmostEqual(instance.model_time, 0.99 * 6415.0029 | units.yr, 3)
        
        print("evolve_model with end_time: take timesteps, until end_time is reached exactly")
        instance.evolve_model(15000 | units.yr)
        self.assertAlmostEqual(instance.particles.age, [100000.0, 17677.6695, 6415.0029*(1+1.2+1.44)] | units.yr, 3)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.2*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 4)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3)
        
        print("evolve_model with keep_synchronous: use non-shared timestep, particle ages will typically diverge")
        instance.evolve_model(keep_synchronous = False)
        self.assertAlmostRelativeEquals(instance.particles.age, ([100000.0, 17677.6695, 6415.0029*(1+1.2+1.44)] | units.yr)
            + 1.2*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 5)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.44*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 4)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3) # Unchanged!
        instance.stop()
    
    def test12(self):
        print("Test for importing new stellar models")
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.parameters.stabilize_new_stellar_model_flag = False
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model()
        
        number_of_zones = instance.particles[0].get_number_of_zones()
        composition     = instance.particles[0].get_chemical_abundance_profiles(number_of_zones = number_of_zones)
        instance.new_particle_from_model(dict(
            mass = instance.particles[0].get_cumulative_mass_profile(number_of_zones = number_of_zones) * instance.particles[0].mass,
            radius = instance.particles[0].get_radius_profile(number_of_zones = number_of_zones),
            rho    = instance.particles[0].get_density_profile(number_of_zones = number_of_zones),
            temperature = instance.particles[0].get_temperature_profile(number_of_zones = number_of_zones),
            luminosity  = instance.particles[0].get_luminosity_profile(number_of_zones = number_of_zones),
            X_H  = composition[0],
            X_He = composition[1] + composition[2],
            X_C  = composition[3],
            X_N  = composition[4],
            X_O  = composition[5],
            X_Ne = composition[6],
            X_Mg = composition[7],
            X_Si = composition[7]*0.0,
            X_Fe = composition[7]*0.0), 10.0 | units.Myr)
        self.assertEqual(len(instance.particles), 2)
        self.assertEqual(len(instance.imported_stars), 1)
        self.assertEqual(instance.imported_stars[0].get_number_of_zones(), number_of_zones)
        self.assertIsOfOrder(instance.imported_stars[0].get_radius_profile()[-1],          1.0 | units.RSun)
        self.assertIsOfOrder(instance.imported_stars[0].get_temperature_profile()[0],  1.0e7 | units.K)
        self.assertIsOfOrder(instance.imported_stars[0].get_pressure_profile()[0],      1.0e17 | units.barye)
        self.assertAlmostEqual(instance.imported_stars[0].get_mass_profile(), 
                               instance.native_stars[0].get_mass_profile())
        self.assertAlmostRelativeEqual(instance.imported_stars[0].get_pressure_profile(), 
                               instance.native_stars[0].get_pressure_profile(),7)
        self.assertAlmostEqual(instance.imported_stars[0].get_radius_profile(), 
                               instance.native_stars[0].get_radius_profile())
        self.assertAlmostEqual(instance.imported_stars[0].get_temperature_profile(), 
                               instance.native_stars[0].get_temperature_profile())
        
        print(instance.particles)
        instance.evolve_model(keep_synchronous = False)
        print(instance.particles)
        instance.stop()
    
    def slowtest13(self):
        print("Testing MESA wind parameters...")
        stars = Particles(9)
        stars.mass = 10.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.parameters.reimers_wind_efficiency = 0.5
        instance.parameters.blocker_wind_efficiency = 0.1
        instance.parameters.de_jager_wind_efficiency = 0.8
        instance.parameters.dutch_wind_efficiency = 0.8
        instance.commit_parameters()
        for wind_scheme in [0, 1, 2, 3, 4]:
            instance.parameters.RGB_wind_scheme = wind_scheme
            instance.recommit_parameters() 
            instance.particles.add_particle(stars[wind_scheme])
        instance.parameters.reimers_wind_efficiency *= 2.0
        instance.parameters.blocker_wind_efficiency *= 2.0
        instance.parameters.de_jager_wind_efficiency *= 2.0
        instance.parameters.dutch_wind_efficiency *= 2.0
        for wind_scheme in [1, 2, 3, 4]:
            instance.parameters.RGB_wind_scheme = wind_scheme
            instance.recommit_parameters() 
            instance.particles.add_particle(stars[wind_scheme+4])
        instance.commit_particles()
        instance.evolve_model(keep_synchronous = False)
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        print(stars)
        self.assertAlmostEqual(stars[0].wind, 0.0 | units.MSun / units.yr)
        self.assertAlmostRelativeEqual(stars[1:5].wind, 
            [4.59318475897e-10, 5.20742729636e-11, 1.05565558121e-09, 3.62519254311e-09] | units.MSun / units.yr, places = 7)
        self.assertAlmostRelativeEqual(stars[5:].wind, 2.0 * stars[1:5].wind, places = 7)
        instance.stop()
    
    def test14(self):
        print("Testing MESA wind parameters... (short version of slowtest13)")
        stars = Particles(3)
        stars.mass = 10.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        instance.initialize_code()
        instance.parameters.RGB_wind_scheme = 0
        instance.commit_parameters()
        instance.particles.add_particle(stars[0])
        instance.parameters.RGB_wind_scheme = 1
        for i, wind_efficiency in enumerate([0.5, 1.0]):
            instance.parameters.reimers_wind_efficiency = wind_efficiency
            instance.recommit_parameters()
            instance.particles.add_particle(stars[i+1])
        instance.commit_particles()
        instance.evolve_model(keep_synchronous = False)
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        self.assertAlmostEqual(stars[0].wind, 0.0 | units.MSun / units.yr)
        self.assertAlmostRelativeEqual(stars[1].wind, 4.59318475897e-10 | units.MSun / units.yr, places = 1)
        self.assertAlmostRelativeEqual(stars[2].wind, 2.0 * stars[1].wind, places = 7)
        instance.stop()
    
    def test15(self):
        print("Testing MESA states")
        stars = Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        
        print("First do everything manually:", end=' ')
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particle(stars[0])
        instance.commit_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEqual(instance.get_name_of_current_state(), 'END')
        instance.stop()
        print("ok")

        print("initialize_code(), commit_parameters(), (re)commit_particles(), " \
            "and cleanup_code() should be called automatically:", end=' ')
        instance = MESA()
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.RGB_wind_scheme = 1
        instance.parameters.reimers_wind_efficiency = 0.5
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particle(stars[0])
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.particles.add_particle(stars[1])
        self.assertEqual(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.particles[0].mass
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
        print("ok")
    
    def test16(self):
        print("Testing basic operations: evolve_one_step and evolve_for")
        stars = Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance_of_an_optional_code(MESA)#, redirection = 'none')
        if instance is None:
            print("MESA was not built. Skipping test.")
            return
        
        se_stars = instance.particles.add_particles(stars)
        
        for i in range(3):
            se_stars[0].evolve_one_step()
        self.assertAlmostEqual(se_stars.age, [364000.0, 0] | units.yr)
        
        se_stars[1].evolve_for(se_stars[0].age)
        self.assertAlmostRelativeEqual(se_stars[0].age,         se_stars[1].age)
        self.assertAlmostRelativeEqual(se_stars[0].luminosity,  se_stars[1].luminosity, 2)
        self.assertAlmostRelativeEqual(se_stars[0].radius,      se_stars[1].radius, 2)
        self.assertAlmostRelativeEqual(se_stars[0].temperature, se_stars[1].temperature, 2)
        instance.stop()
    
    def test17(self):
        print("MESA validation")
        
        mesa_src_path = os.path.join(os.path.dirname(sys.modules[MESA.__module__].__file__), 'src', 'mesa')
        mesa_star_path = os.path.join(mesa_src_path, 'star', 'test', 'star')
        
        if not os.path.exists(mesa_star_path) or not os.access(mesa_star_path, os.X_OK):
            self.skip("no mesa executable in binary distribution, test cannot run")
        
        number_of_steps = 11#4
        star = Particle()
        star.mass = 1.0 | units.MSun
        
        testpath = get_path_to_results()
        outputfile_name = os.path.join(testpath, "mesa_output")
       
        instance = self.new_instance_of_an_optional_code(MESA, redirection = 'file', redirect_file = outputfile_name)
        se_star = instance.particles.add_particle(star)
        for i in range(number_of_steps):
            se_star.evolve_one_step()
        instance.stop()
      
        rfile = open(outputfile_name, 'r')
        amuse_output = rfile.read()
        rfile.close()
        
        
        #generate the inlist for the (stand-alone) MESA star run:
        instance = self.new_instance_of_an_optional_code(MESA, redirection = 'null')
        with open(instance.default_path_to_inlist, 'r') as default_inlist:
            with open(os.path.join(testpath, 'inlist'), 'w') as test_inlist:
                for one_line in default_inlist.readlines():
                    if ("max_model_number" in one_line):
                        test_inlist.write("         max_model_number = "+str(number_of_steps+1)+"\n")
                    elif ("mesa_data_dir" in one_line):
                        test_inlist.write("      mesa_data_dir = '"+os.path.join(mesa_src_path, 'data')+"'\n")
                    elif ("zams_filename" in one_line):
                        test_inlist.write("      zams_filename = '"+os.path.join(instance.get_output_directory(), "star_data", "starting_models", "zams_z20m3.data")+"'\n")
                    else:
                        test_inlist.write(one_line)
        instance.stop()
        
        (stdout, stderr) = Popen([mesa_star_path], cwd = testpath, stdin = PIPE, stdout = PIPE, stderr = PIPE).communicate()
        self.assertEqual(stderr, "")
        for i, line in enumerate(stdout.splitlines()):
            #print i, line, line in amuse_output
            if i == 52 + 4 * number_of_steps + 8 * (number_of_steps/10) + (3 if number_of_steps > 55 else 0):
                self.assertEqual(line, "stop because model_number >= max_model_number")
                break
            else:
                self.assertTrue(line in amuse_output)
    
    def test18(self):
        print("Testing MESA mass_change (User-specified wind/accretion)")
        instance = self.new_instance_of_an_optional_code(MESA)
        instance.parameters.RGB_wind_scheme = 0 # must be turned off for user-specified rates
        instance.parameters.AGB_wind_scheme = 0 # must be turned off for user-specified rates
        
        star = instance.particles.add_particle(Particle(mass=1|units.MSun))
        star.mass_change = 1.0e-8 | units.MSun / units.yr # positive -> accretion
        star.evolve_one_step()
        
        self.assertAlmostRelativeEqual(star.mass_change, 1.0e-8 | units.MSun / units.yr)
        self.assertAlmostRelativeEqual(star.wind, -1.0e-8 | units.MSun / units.yr, 3)
        self.assertAlmostRelativeEqual(star.age, 1.0e5 | units.yr)
        self.assertAlmostRelativeEqual(star.mass, 1.0010 | units.MSun)
        
        star.mass_change = -1.0e-8 | units.MSun / units.yr # negative -> wind
        star.evolve_one_step()
        
        self.assertAlmostRelativeEqual(star.mass_change, -1.0e-8 | units.MSun / units.yr)
        self.assertAlmostRelativeEqual(star.wind, 1.0e-8 | units.MSun / units.yr, 3)
        self.assertAlmostRelativeEqual(star.age, 1.8e5 | units.yr)
        self.assertAlmostRelativeEqual(star.mass, 1.0002 | units.MSun)
        print(star.as_set())
        instance.stop()
    
    def slowtest19a(self):
        print("Testing MESA core mass")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.particles.add_particle(Particle(mass=3|units.MSun))
        star.evolve_for(330 | units.Myr)
        for i in range(10):
            star.evolve_for(10 | units.Myr)
            index = numpy.searchsorted(star.get_chemical_abundance_profiles(number_of_species=1)[0], 1.0e-4)
            h_poor_mass = EnclosedMassInterpolator(radii=star.get_radius_profile(), densities=star.get_density_profile()).enclosed_mass[index].as_quantity_in(units.MSun)
            print(h_poor_mass, star.core_mass)
            self.assertAlmostEqual(star.core_mass, h_poor_mass, 2)
        instance.stop()
    
    def test19b(self):
        print("Testing MESA core mass (short version of slowtest19a)")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.particles.add_particle(Particle(mass=3|units.MSun))
        star.evolve_one_step()
        index = numpy.searchsorted(star.get_chemical_abundance_profiles(number_of_species=1)[0], 1.0e-4)
        h_poor_mass = EnclosedMassInterpolator(radii=star.get_radius_profile(), densities=star.get_density_profile()).enclosed_mass[index].as_quantity_in(units.MSun)
        self.assertAlmostEqual(star.core_mass, h_poor_mass, 2)
        instance.stop()
    
    def test20(self):
        print("Testing MESA pre-main-sequence star")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.pre_ms_stars.add_particle(Particle(mass=1.0|units.MSun))
        
        self.assertAlmostEqual(star.time_step, 1.0e-3 | units.yr)
        star.evolve_one_step()
        self.assertAlmostEqual(star.age, 1.0e-3 | units.yr)
        
        instance.evolve_model(1.0 | units.yr)
        self.assertEqual(star.stellar_type, 17 | units.stellar_type)
        self.assertEqual(str(star.stellar_type), "Pre-main-sequence Star")
        self.assertTrue(star.age > 1.0 | units.yr)
        self.assertTrue(star.temperature < 4500 | units.K)
        self.assertTrue(star.luminosity > 10 | units.LSun)
        self.assertTrue(star.radius > 2 | units.RSun)
        instance.stop()
    
    def slowtest21(self):
        print("Testing MESA calculate_core_mass")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.particles.add_particle(Particle(mass=40|units.MSun))
        instance.evolve_model(4.56|units.Myr)
        total_core_mass = star.calculate_core_mass()
        self.assertAlmostRelativeEqual(star.core_mass, total_core_mass, 2)
        self.assertTrue(star.calculate_core_mass(core_H_abundance_limit=1.0e-2) > total_core_mass)
        self.assertEqual(star.calculate_core_mass(core_H_abundance_limit=1.0e-4), total_core_mass)
        self.assertTrue(star.calculate_core_mass(core_H_abundance_limit=1.0e-6) < total_core_mass)
        self.assertAlmostRelativeEqual(star.calculate_core_mass(core_H_abundance_limit=0.8), star.mass, 2)
        
        h1_core_mass = star.calculate_core_mass(species=["h1"])
        he3_core_mass = star.calculate_core_mass(species=["he3"])
        he4_core_mass = star.calculate_core_mass(species=["he4"])
        c12_core_mass = star.calculate_core_mass(species=["c12"])
        n14_core_mass = star.calculate_core_mass(species=["n14"])
        o16_core_mass = star.calculate_core_mass(species=["o16"])
        ne20_core_mass = star.calculate_core_mass(species=["ne20"])
        mg24_core_mass = star.calculate_core_mass(species=["mg24"])
        metal_core_mass = star.calculate_core_mass(species=["c12", "n14", "o16", "ne20", "mg24"])
        print(h1_core_mass, he3_core_mass, he4_core_mass, c12_core_mass, n14_core_mass, o16_core_mass)
        print(ne20_core_mass, mg24_core_mass, metal_core_mass)
        self.assertAlmostRelativeEqual(star.core_mass, total_core_mass, 2)
        instance.stop()
        self.assertAlmostRelativeEqual(total_core_mass, he4_core_mass, 1)
        self.assertAlmostRelativeEqual(total_core_mass, he4_core_mass + metal_core_mass, 4)
        self.assertAlmostRelativeEqual(total_core_mass, he4_core_mass + metal_core_mass + h1_core_mass, 7)
        self.assertAlmostRelativeEqual(metal_core_mass, 
            c12_core_mass + n14_core_mass + o16_core_mass + ne20_core_mass + mg24_core_mass, 7)
        self.assertAlmostEqual(he3_core_mass, 0 | units.MSun)
    
    def test22(self):
        print("Testing MESA calculate_core_mass (short version of slowtest21)")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.particles.add_particle(Particle(mass=1|units.MSun))
        instance.evolve_model(0.3|units.Gyr) # VERY short, for test speed up
        central_hydrogen_abundance = star.get_chemical_abundance_profiles()[0][0]
        self.assertTrue(central_hydrogen_abundance < 0.68) # some hydrogen is burned
        self.assertTrue(central_hydrogen_abundance > 0.67) # ... but not that much yet
        self.assertEqual(star.calculate_core_mass(core_H_abundance_limit=0.67), 0 | units.MSun)
        self.assertAlmostEqual(star.calculate_core_mass(core_H_abundance_limit=0.70), 1 | units.MSun, 3)
        
        # For test speed up, we use a weird core_H_abundance_limit to define the "hydrogen exhausted core"
        limit = 0.68
        expected_core_mass = 0.01786033709 | units.MSun
        self.assertAlmostEqual(star.calculate_core_mass(core_H_abundance_limit=limit), expected_core_mass, 3)
        
        h1_core_mass = star.calculate_core_mass(species=["h1"], core_H_abundance_limit=limit)
        he3_core_mass = star.calculate_core_mass(species=["he3"], core_H_abundance_limit=limit)
        he4_core_mass = star.calculate_core_mass(species=["he4"], core_H_abundance_limit=limit)
        c12_core_mass = star.calculate_core_mass(species=["c12"], core_H_abundance_limit=limit)
        n14_core_mass = star.calculate_core_mass(species=["n14"], core_H_abundance_limit=limit)
        o16_core_mass = star.calculate_core_mass(species=["o16"], core_H_abundance_limit=limit)
        ne20_core_mass = star.calculate_core_mass(species=["ne20"], core_H_abundance_limit=limit)
        mg24_core_mass = star.calculate_core_mass(species=["mg24"], core_H_abundance_limit=limit)
        metal_core_mass = star.calculate_core_mass(species=["c12", "n14", "o16", "ne20", "mg24"], core_H_abundance_limit=limit)
        instance.stop()
        self.assertAlmostRelativeEqual(h1_core_mass, expected_core_mass*0.68, 2)
        self.assertAlmostRelativeEqual(he4_core_mass, expected_core_mass*0.30, 2)
        self.assertAlmostRelativeEqual(metal_core_mass, expected_core_mass*0.02, 1)
        self.assertAlmostRelativeEqual(expected_core_mass, he4_core_mass + he3_core_mass + metal_core_mass + h1_core_mass, 7)
        self.assertAlmostRelativeEqual(metal_core_mass, 
            c12_core_mass + n14_core_mass + o16_core_mass + ne20_core_mass + mg24_core_mass, 7)
        self.assertAlmostEqual(he3_core_mass, 0 | units.MSun, 5)
    
    def test23(self):
        print("Testing MESA central_temperature and central_density")
        instance = self.new_instance_of_an_optional_code(MESA)
        stars = instance.particles.add_particles(Particles(mass=[0.1, 1, 10]|units.MSun))
        self.assertIsOfOrder(stars.central_temperature, [4e6, 13e6, 31e6] | units.K)
        self.assertIsOfOrder(stars.central_density, [400, 77, 9] | units.g * units.cm**-3)
        instance.stop()
    
    def test24(self):
        print("Testing MESA calculate_helium_exhausted_core_mass")
        instance = self.new_instance_of_an_optional_code(MESA)
        star = instance.particles.add_particle(Particle(mass=2|units.MSun))
        
        composition = star.get_chemical_abundance_profiles()
        # Mimic hydrogen exhausted core:
        composition[2, :100] = composition[2, :100] + composition[0, :100]
        composition[0, :100] = 0
        # Mimic helium exhausted core:
        carbon_oxygen = composition[:6, :50].sum(axis=0)
        composition[3, :50] = carbon_oxygen * 0.6
        composition[5, :50] = carbon_oxygen * 0.4
        composition[1:3, :50] = 0
        composition[4, :50] = 0
        star.set_chemical_abundance_profiles(composition)
        
        self.assertAlmostRelativeEqual(star.calculate_core_mass(), 
            star.mass * star.get_cumulative_mass_profile()[100], 1)
        self.assertAlmostRelativeEqual(star.calculate_helium_exhausted_core_mass(), 
            star.mass * star.get_cumulative_mass_profile()[50], 1)
        
        core_mass = star.calculate_helium_exhausted_core_mass(split_species=False)
        core_mass_by_species = star.calculate_helium_exhausted_core_mass(split_species=True)
        carbon_mass_in_core, oxygen_mass_in_core = star.calculate_helium_exhausted_core_mass(
            split_species=True, species=["c12", "o16"])
        
        self.assertEqual(len(core_mass_by_species), len(star.get_names_of_species()))
        instance.stop()
        self.assertAlmostRelativeEqual(core_mass, core_mass_by_species.sum())
        self.assertEqual(core_mass_by_species[0:3].sum(), 0 | units.MSun)
        self.assertEqual(core_mass_by_species[4], 0 | units.MSun)
        self.assertEqual(core_mass_by_species[3], carbon_mass_in_core)
        self.assertEqual(core_mass_by_species[5], oxygen_mass_in_core)
    
    def test25(self):
        print("Testing MESA accretion")
        instance = self.new_instance_of_an_optional_code(MESA)
        instance.parameters.RGB_wind_scheme = 0
        instance.parameters.AGB_wind_scheme = 0
        star = instance.particles.add_particle(Particle(mass=2|units.MSun))
        
        self.assertEqual(star.get_accrete_same_as_surface(), 1)
        star.set_accrete_same_as_surface(0)
        self.assertEqual(star.get_accrete_same_as_surface(), 0)
        
        self.assertEqual(star.get_accrete_composition_non_metals(), [-1.0, -1.0, -1.0, -1.0])
        self.assertEqual(star.get_accrete_composition_metals_identifier(), -1)
        self.assertEqual(star.get_accrete_composition_metals(), [-1.0]*28)
        print("Accreting 75% deuterium", end=' ')
        composition_light = [0, 0.75, 0, 0]
        print("and 25% iron")
        composition_metals = [0]*23 + [1.0] + [0]*4
        star.set_accrete_composition_non_metals(*composition_light)
        star.set_accrete_composition_metals_identifier(0) # i.e. specified below:
        star.set_accrete_composition_metals(*composition_metals)
        self.assertEqual(star.get_accrete_composition_non_metals(), composition_light)
        self.assertEqual(star.get_accrete_composition_metals_identifier(), 0)
        self.assertEqual(star.get_accrete_composition_metals(), composition_metals)
        
        star.mass_change = 1.0e-8 | units.MSun / units.yr
        star.time_step = 0.1 | units.yr
        instance.evolve_model(1 | units.yr)
        composition = star.get_chemical_abundance_profiles()
        species = star.get_names_of_species()
        print("Both deuterium and iron are not in the current net,")
        print("so have been added to {0} and {1}".format(species[0], species[-1]))
        self.assertEqual(composition[:, -1], [0.75, 0, 0, 0, 0, 0, 0, 0.25])
        instance.stop()
    


