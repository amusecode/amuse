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
        print "Testing initialization of the interface..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        inlist_path = instance.default_path_to_inlist
        #print "Path to inlist: ", inlist_path
        MESA_data_path = instance.default_path_to_MESA_data
        #print "Path to MESA data directory: ", MESA_data_path
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        instance.stop()
        
    def slowtest2(self):
        print "Testing get/set of metallicity (tests new ZAMS model implicitly)..."
        print "The first time this test will take quite some time" \
            " to generate new starting models."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        (metallicity, error) = instance.get_metallicity()
        self.assertEquals(0, error)
        self.assertEquals(0.02, metallicity)
        for x in [0.04, 0.02, 0.01, 0.00]:
            error = instance.set_metallicity(x)
            self.assertEquals(0, error)
            (metallicity, error) = instance.get_metallicity()
            self.assertEquals(0, error)
            self.assertEquals(x, metallicity)
        instance.stop()
    
    def test3(self):
        print "Testing basic operations: new_particle..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        (maximum_number_of_stars, error) = instance.get_maximum_number_of_stars()
        self.assertEquals(0, error)
        self.assertEqual(maximum_number_of_stars,1000)
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        number_of_stars = 1
        for i in range(number_of_stars):
            #print i
            (index_of_the_star, error) = instance.new_particle(0.5+i*1.0/number_of_stars)
            self.assertEquals(0, error)
            self.assertEqual(index_of_the_star,i+1)
        #import time    # To check the amount of memory is used ...
        #time.sleep(10) # when a large number of stars is created.
        instance.stop()
        
    def test4(self):
        print "Testing basic operations: evolve..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        status = instance.initialize_code()
        (index_of_the_star, error) = instance.new_particle(1.0)
        self.assertEquals(0, error)
        self.assertEqual(index_of_the_star, 1)
        self.assertEquals(0, instance.commit_particles())
        
        initial_dt = 1.0e5
        dt_factor = 1.2
        self.assertEqual([initial_dt, 0], instance.get_time_step(index_of_the_star).values())
        self.assertEquals(0, instance.evolve_one_step(index_of_the_star))
        self.assertEqual([initial_dt, 0], instance.get_age(index_of_the_star).values())
        
        self.assertEquals(0, instance.evolve_for(index_of_the_star, initial_dt))
        # nothing happens; the target_time of the star is increased with initial_dt, but star already has this age
        self.assertEqual([initial_dt, 0], instance.get_age(index_of_the_star).values())
        
        target_end_time = 3.0e5 # (years)
        self.assertEquals(0, instance.evolve_for(index_of_the_star, target_end_time-initial_dt))
        self.assertEqual([initial_dt*(1 + dt_factor + dt_factor**2), 0], instance.get_age(index_of_the_star).values())
        self.assertEqual([round(initial_dt*dt_factor**3), 0], instance.get_time_step(index_of_the_star).values())
        self.assertTrue(instance.get_age(index_of_the_star)['age'] >= target_end_time)
        
        (L_of_the_star, error) = instance.get_luminosity(index_of_the_star)
        self.assertEquals(0, error)
        self.assertAlmostEqual(L_of_the_star,0.725,1)
        (M_of_the_star, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)
        self.assertAlmostEqual(M_of_the_star,1.000,3)
        (T_of_the_star, error) = instance.get_temperature(index_of_the_star)
        self.assertEquals(0, error)
        self.assertAlmostEqual(T_of_the_star,5650.998,-2)
        instance.stop()
    
    def slowtest5(self):
        print "Testing evolve with varying Z (tests new ZAMS model implicitly)..."
        print "If the required starting models do not exist, this test will " \
            "take quite some time to generate them."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        status = instance.initialize_code()
        self.assertEqual(status,0)
        metallicities = [0.00, 0.01, 0.02, 0.04]
        luminosities = [1.717, 0.938, 0.725, 0.592]
        for (i,(Z,L)) in enumerate(zip(metallicities, luminosities)):
            status = instance.set_metallicity(Z)
            self.assertEquals(0, status)
            (index_of_the_star, status) = instance.new_particle(1.0)
            self.assertEquals(0, status)
            self.assertEqual(index_of_the_star,i+1)
            instance.evolve_for(index_of_the_star, 5.0e5)
            (L_of_the_star, status) = instance.get_luminosity(index_of_the_star)
            self.assertEquals(0, status)
            self.assertAlmostEqual(L_of_the_star,L,1)
        instance.stop()
    
    def test6(self):
        print "Testing MESA stop conditions..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        (value, error) = instance.get_max_age_stop_condition()
        self.assertEquals(0, error) 
        self.assertEquals(1.0e12, value)
        for x in range(10,14):
            error = instance.set_max_age_stop_condition(10 ** x)
            self.assertEquals(0, error)
            (value, error) = instance.get_max_age_stop_condition()
            self.assertEquals(0, error)
            self.assertEquals(10 ** x, value)
        (value, error) = instance.get_min_timestep_stop_condition()
        self.assertEquals(0, error)
        self.assertEquals(1.0e-6, value)
        for x in range(-9,-2):
            error = instance.set_min_timestep_stop_condition(10.0 ** x)
            self.assertEquals(0, error)
            (value, error) = instance.get_min_timestep_stop_condition()
            self.assertEquals(0, error)
            self.assertEquals(10.0 ** x, value)
        instance.stop()
    
    def test7(self):
        print "Testing MESA parameters..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        (value, error) = instance.get_mixing_length_ratio()
        self.assertEquals(0, error) 
        self.assertEquals(2.0, value)
        for x in [1.0, 1.5, 3.0]:
            error = instance.set_mixing_length_ratio(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_mixing_length_ratio()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        (value, error) = instance.get_semi_convection_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.0, value)
        for x in [0.1, 0.04, 0.001]:
            error = instance.set_semi_convection_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_semi_convection_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        instance.stop()
    
    def test8(self):
        print "Testing MESA wind parameters..."
        instance = self.new_instance(MESAInterface)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        (value, error) = instance.get_RGB_wind_scheme()
        self.assertEquals(0, error) 
        self.assertEquals(1, value)
        for x in range(6):
            error = instance.set_RGB_wind_scheme(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_RGB_wind_scheme()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_AGB_wind_scheme()
        self.assertEquals(0, error) 
        self.assertEquals(1, value)
        for x in range(6):
            error = instance.set_AGB_wind_scheme(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_AGB_wind_scheme()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
            
        (value, error) = instance.get_reimers_wind_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.5, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_reimers_wind_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_reimers_wind_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_blocker_wind_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.1, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_blocker_wind_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_blocker_wind_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_de_jager_wind_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.8, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_de_jager_wind_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_de_jager_wind_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        (value, error) = instance.get_dutch_wind_efficiency()
        self.assertEquals(0, error)
        self.assertEquals(0.8, value)
        for x in [0.0, 0.1, 0.5, 1.0]:
            error = instance.set_dutch_wind_efficiency(x)
            self.assertEquals(0, error)
            (value, error) = instance.get_dutch_wind_efficiency()
            self.assertEquals(0, error)
            self.assertEquals(x, value)
        
        instance.stop()


class TestMESA(TestWithMPI):
    
    def test1(self):
        print "Testing initialization and default MESA parameters..."
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.set_MESA_paths(instance.default_path_to_inlist, 
            instance.default_path_to_MESA_data, instance.get_data_directory())
        instance.initialize_code()
        self.assertEquals(0.02 | units.no_unit, instance.parameters.metallicity)
        self.assertEquals(1.0e12 | units.yr, instance.parameters.max_age_stop_condition)
        instance.parameters.max_age_stop_condition = 1.0e2 | units.Myr
        self.assertEquals(1.0e2 | units.Myr, instance.parameters.max_age_stop_condition)
        instance.stop()
    
    def test2(self):
        print "Testing basic operations: evolve and get_..."
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        instance.evolve_for(index_of_the_star, age_of_the_star)
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
        print "Testing basic operations: evolve_model and channels..."
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        self.assertEquals(stars[0].mass, mass)
        self.assertAlmostEquals(stars[0].luminosity, 5841. | units.LSun, 0)
        instance.evolve_model()
        from_code_to_model.copy()
        self.assertAlmostEquals(stars[0].mass, mass, 5)
        self.assertAlmostEquals(stars[0].luminosity, 5820.85 | units.LSun, 0)
        instance.stop()
    
    def slowtest4(self):
        print "Testing stellar type..."
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
                print (star.age, star.mass, star.stellar_type)
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
            self.assertEquals(str(result[2]), expected)
        
        instance.stop()
    
    def test5(self):
        print "Testing evolve_model for particle set..."
        instance = self.new_instance(MESA)
        masses = [0.5, 1.0] | units.MSun
        max_age = 0.6 | units.Myr
        number_of_stars=len(masses)
        stars =  Particles(number_of_stars)
        stars.mass = masses
        instance.initialize_code()
        self.assertEqual(instance.parameters.max_age_stop_condition, 1e6 | units.Myr)
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
        print "Test for obtaining the stellar structure model"
        stars = Particles(2)
        stars.mass = [1.0, 10.0] | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model()
        self.assertEquals(instance.particles.get_number_of_zones(), [479, 985])
        self.assertEquals(len(instance.particles[0].get_mass_profile()), 479)
        self.assertAlmostEquals(instance.particles[0].get_mass_profile().sum(), 1.0)
        self.assertRaises(AmuseException, instance.particles.get_mass_profile, 
            expected_message = "Querying mass profiles of more than one particle at a time is not supported.")
        print instance.particles
        self.assertEquals(len(instance.particles[1].get_density_profile()), 985)
        self.assertIsOfOrder(instance.particles[0].get_radius_profile()[-1],          1.0 | units.RSun)
        self.assertIsOfOrder(instance.particles[0].get_temperature_profile()[0],  1.0e7 | units.K)
        self.assertIsOfOrder(instance.particles[0].get_luminosity_profile()[-1],      1.0 | units.LSun)
        self.assertIsOfOrder(instance.particles[0].get_pressure_profile()[0],      1.0e17 | units.barye)
        delta_mass = instance.particles[0].get_mass_profile() * instance.particles[0].mass
        radius1 = instance.particles[0].get_radius_profile()
        radius2 = radius1[:-1]
        radius2.prepend(0|units.m)
        delta_radius_cubed = (radius1**3 - radius2**3)
        self.assertAlmostEquals(instance.particles[0].get_density_profile() / 
            (delta_mass/(4./3.*numpy.pi*delta_radius_cubed)), [1]*479, places=3)
        self.assertAlmostEquals(instance.particles[1].get_mu_profile(), [0.62]*985 | units.amu, places=1)
        instance.stop()
    
    def test7(self):
        print "Test for obtaining the stellar composition structure"
        stars = Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        self.assertEquals(number_of_zones,    479)
        self.assertEquals(number_of_species,    8)
        self.assertEquals(len(species_names),  number_of_species)
        self.assertEquals(len(composition),    number_of_species)
        self.assertEquals(len(composition[0]), number_of_zones)
        self.assertEquals(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
        self.assertEquals(species_IDs,   [2,    5,     6,     38,    51,    69,    114,    168])
        self.assertAlmostEquals(species_masses, [1.0078250, 3.0160293, 4.0026032, 12.0, 
                                14.0030740, 15.9949146, 19.9924401, 23.9850417] | units.amu, places=5)
        self.assertAlmostEquals(composition[ :1,-1].sum(),  0.7)
        self.assertAlmostEquals(composition[1:3,-1].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEquals(composition[3: ,-1].sum(),  instance.parameters.metallicity)
        self.assertAlmostEquals(composition.sum(axis=0), [1.0]*number_of_zones)
        instance.stop()
    
    def slowtest8(self):
        print "Test for obtaining the stellar composition structure - evolved star with zero metalicity"
        stars = Particles(1)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        self.assertEquals(number_of_zones,    578)
        self.assertEquals(number_of_species,    8)
        self.assertEquals(len(species_names),  number_of_species)
        self.assertEquals(len(composition),    number_of_species)
        self.assertEquals(len(composition[0]), number_of_zones)
        self.assertEquals(species_names, ['h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'])
        self.assertEquals(species_IDs,   [2,    5,     6,     38,    51,    69,    114,    168])
        self.assertAlmostEquals(composition[ :1,number_of_zones-1].sum(),  0.76 | units.none)
        self.assertAlmostEquals(composition[1:3,number_of_zones-1].sum(),  0.24 | units.none)
        self.assertAlmostEquals(composition[3: ,number_of_zones-1].sum(),  0.00 | units.none)
        self.assertAlmostEquals(composition.sum(axis=0), [1.0]*number_of_zones | units.none)
        self.assertAlmostEquals(composition[ :1,0].sum(),  0.00 | units.none)
        self.assertAlmostEquals(composition[1:3,0].sum(),  1.00 | units.none)
        self.assertAlmostEquals(composition[3: ,0].sum(),  0.00 | units.none)
        instance.stop()
    
    def test9(self):
        print "Test for changing the stellar structure model"
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model()

        density_profile = instance.particles[0].get_density_profile()
        
        self.assertRaises(AmuseException, instance.particles[0].set_density_profile, density_profile[2:], 
            expected_message = "The length of the supplied vector (477) does not match the number of "
            "mesh zones of the star (479).")
        
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
        print "Test for changing the stellar composition"
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)#, redirection = 'none')
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(star)
        instance.commit_particles()
        instance.evolve_model(0.3 | units.Myr)
        
        composition = instance.particles[0].get_chemical_abundance_profiles()
        k_surface = -1 # index to the outer mesh cell (surface)
        
        self.assertAlmostEquals(composition[ :1, k_surface].sum(),  0.7)
        self.assertAlmostEquals(composition[1:3, k_surface].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEquals(composition[3: , k_surface].sum(),  instance.parameters.metallicity)
        
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
        
        self.assertAlmostEquals(composition[ :2, k_surface].sum(),  (0.3) - instance.parameters.metallicity)
        self.assertAlmostEquals(composition[2:3, k_surface].sum(),  0.7)
        self.assertAlmostEquals(composition[3: , k_surface].sum(),  instance.parameters.metallicity)
        self.assertAlmostEquals(composition.sum(axis=0), 1.0)
        
        self.assertRaises(AmuseException, instance.particles[0].set_chemical_abundance_profiles, composition[:7], 
            expected_message = "The length of the supplied vector (7) does not match the number of "
            "chemical species of the star (8).")
        instance.stop()
    
    def test11(self):
        print "Test evolve_model optional arguments: end_time and keep_synchronous"
        stars = Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        instance.initialize_code()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        self.assertAlmostEqual(instance.particles.age, [0.0, 0.0, 0.0] | units.yr)
        self.assertAlmostEqual(instance.particles.time_step, [100000.0, 17677.6695, 6415.0029] | units.yr, 3)
        
        print "evolve_model without arguments: use shared timestep = 0.99*min(particles.time_step)"
        instance.evolve_model()
        self.assertAlmostEqual(instance.particles.age, [100000.0, 17677.6695, 6415.0029] | units.yr, 3)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.2*([100000.0, 17677.6695, 6415.0029] | units.yr), 6)
        self.assertAlmostEqual(instance.model_time, 0.99 * 6415.0029 | units.yr, 3)
        
        print "evolve_model with end_time: take timesteps, until end_time is reached exactly"
        instance.evolve_model(15000 | units.yr)
        self.assertAlmostEqual(instance.particles.age, [100000.0, 17677.6695, 6415.0029*(1+1.2+1.44)] | units.yr, 3)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.2*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 4)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3)
        
        print "evolve_model with keep_synchronous: use non-shared timestep, particle ages will typically diverge"
        instance.evolve_model(keep_synchronous = False)
        self.assertAlmostRelativeEquals(instance.particles.age, ([100000.0, 17677.6695, 6415.0029*(1+1.2+1.44)] | units.yr)
            + 1.2*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 5)
        self.assertAlmostRelativeEquals(instance.particles.time_step, 1.44*([100000.0, 17677.6695, 1.44*6415.0029] | units.yr), 4)
        self.assertAlmostEqual(instance.model_time, 15000.0 | units.yr, 3) # Unchanged!
        instance.stop()
    
    def test12(self):
        print "Test for importing new stellar models"
        star = Particles(1)
        star.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        
        print instance.particles
        instance.evolve_model(keep_synchronous = False)
        print instance.particles
        instance.stop()
    
    def slowtest13(self):
        print "Testing MESA wind parameters..."
        stars = Particles(9)
        stars.mass = 10.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        print stars
        self.assertAlmostEqual(stars[0].wind, 0.0 | units.MSun / units.yr)
        self.assertAlmostRelativeEqual(stars[1:5].wind, 
            [4.59318475897e-10, 5.20742729636e-11, 1.05565558121e-09, 3.62519254311e-09] | units.MSun / units.yr, places = 7)
        self.assertAlmostRelativeEqual(stars[5:].wind, 2.0 * stars[1:5].wind, places = 7)
        instance.stop()
    
    def test14(self):
        print "Testing MESA wind parameters... (short version of slowtest13)"
        stars = Particles(3)
        stars.mass = 10.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
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
        self.assertAlmostRelativeEqual(stars[1].wind, 4.59318475897e-10 | units.MSun / units.yr, places = 7)
        self.assertAlmostRelativeEqual(stars[2].wind, 2.0 * stars[1].wind, places = 7)
        instance.stop()
    
    def test15(self):
        print "Testing MESA states"
        stars = Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)
        if instance is None:
            print "MESA was not built. Skipping test."
            return
        
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
        instance = MESA()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.RGB_wind_scheme = 1
        instance.parameters.reimers_wind_efficiency = 0.5
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
    
    def test16(self):
        print "Testing basic operations: evolve_one_step and evolve_for"
        stars = Particles(2)
        stars.mass = 1.0 | units.MSun
        instance = self.new_instance(MESA)#, redirection = 'none')
        if instance is None:
            print "MESA was not built. Skipping test."
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
        print "MESA validation"
        
        mesa_src_path = os.path.join(os.path.dirname(sys.modules[MESA.__module__].__file__), 'src', 'mesa')
        mesa_star_path = os.path.join(mesa_src_path, 'star', 'test', 'star')
        
        if not os.path.exists(mesa_star_path) or not os.access(mesa_star_path, os.X_OK):
            self.skip("no mesa executable in binary distribution, test cannot run")
        
        number_of_steps = 11#4
        star = Particle()
        star.mass = 1.0 | units.MSun
        
        testpath = get_path_to_results()
        outputfile_name = os.path.join(testpath, "mesa_output")
       
        instance = self.new_instance(MESA, redirection = 'file', redirect_file = outputfile_name)
        se_star = instance.particles.add_particle(star)
        for i in range(number_of_steps):
            se_star.evolve_one_step()
        instance.stop()
      
        rfile = open(outputfile_name, 'r')
        amuse_output = rfile.read()
        rfile.close()
        
        
        #generate the inlist for the (stand-alone) MESA star run:
        instance = self.new_instance(MESA, redirection = 'null')
        with open(instance.default_path_to_inlist, 'r') as default_inlist:
            with open(os.path.join(testpath, 'inlist'), 'w') as test_inlist:
                for one_line in default_inlist.readlines():
                    if ("max_model_number" in one_line):
                        test_inlist.write("         max_model_number = "+str(number_of_steps+1)+"\n")
                    elif ("mesa_data_dir" in one_line):
                        test_inlist.write("      mesa_data_dir = '"+os.path.join(mesa_src_path, 'data')+"'\n")
                    else:
                        test_inlist.write(one_line)
        instance.stop()
        
        (stdout, stderr) = Popen([mesa_star_path], cwd = testpath, stdin = PIPE, stdout = PIPE, stderr = PIPE).communicate()
        self.assertEquals(stderr, "")
        for i, line in enumerate(stdout.splitlines()):
            #print i, line, line in amuse_output
            if i == 52 + 4 * number_of_steps + 8 * (number_of_steps/10) + (3 if number_of_steps > 55 else 0):
                self.assertEqual(line, "stop because model_number >= max_model_number")
                break
            else:
                self.assertTrue(line in amuse_output)
    
    def test18(self):
        print "Testing MESA mass_change (User-specified wind/accretion)"
        instance = self.new_instance(MESA)
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
        print star.as_set()
        instance.stop()
    
    def slowtest19a(self):
        print "Testing MESA core mass"
        instance = self.new_instance(MESA)
        star = instance.particles.add_particle(Particle(mass=3|units.MSun))
        star.evolve_for(330 | units.Myr)
        for i in range(10):
            star.evolve_for(10 | units.Myr)
            index = numpy.searchsorted(star.get_chemical_abundance_profiles(number_of_species=1)[0], 1.0e-4)
            h_poor_mass = EnclosedMassInterpolator(radii=star.get_radius_profile(), densities=star.get_density_profile()).enclosed_mass[index].as_quantity_in(units.MSun)
            print h_poor_mass, star.core_mass
            self.assertAlmostEqual(star.core_mass, h_poor_mass, 2)
        instance.stop()
    
    def test19b(self):
        print "Testing MESA core mass (short version of slowtest19a)"
        instance = self.new_instance(MESA)
        star = instance.particles.add_particle(Particle(mass=3|units.MSun))
        star.evolve_one_step()
        index = numpy.searchsorted(star.get_chemical_abundance_profiles(number_of_species=1)[0], 1.0e-4)
        h_poor_mass = EnclosedMassInterpolator(radii=star.get_radius_profile(), densities=star.get_density_profile()).enclosed_mass[index].as_quantity_in(units.MSun)
        self.assertAlmostEqual(star.core_mass, h_poor_mass, 2)
        instance.stop()
    
    def test20(self):
        print "Testing MESA pre-main-sequence star"
        instance = self.new_instance(MESA)
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
    


