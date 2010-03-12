from legacy_support import TestWithMPI
import sys
import os.path

from amuse.legacy.mesa.interface import MESA, MESAInterface

from amuse.support.data import core
from amuse.support.units import units
from amuse.legacy.support import channel

class TestMESAInterface(TestWithMPI):
    
#    def worker_found(self):
#       self.name_of_the_worker
       
    def test1(self):
        print "Testing initialization of the interface..."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            inlist_path = instance.default_path_to_inlist
            #print "Path to inlist: ", inlist_path
            MESA_data_path = instance.default_path_to_MESA_data
            #print "Path to MESA data directory: ", MESA_data_path
            status = instance.initialize(inlist_path, MESA_data_path)
            self.assertEqual(status,0)
        else:
            print "MESA was not built. Skipping test."
        del instance
        
    def xtest2(self):
        print "Testing get/set of interface parameters..."
        print "The first time this test will take quite some time"
        print "to generate new starting models."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
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
        else:
            print "MESA was not built. Skipping test."
        del instance
    
    def test3(self):
        print "Testing basic operations: new_particle..."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            (maximum_number_of_stars, error) = instance.get_maximum_number_of_stars()
            self.assertEquals(0, error)
            self.assertEqual(maximum_number_of_stars,1000)
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
            self.assertEqual(status,0)
            number_of_stars = 10
            for i in range(number_of_stars):
                #print i
                (index_of_the_star, error) = instance.new_particle(0.5+i*1.0/number_of_stars)
                self.assertEquals(0, error)
                self.assertEqual(index_of_the_star,i+1)
            #import time    # To check the amount of memory is used ...
            #time.sleep(10) # when a large number of stars is created.
        else:
            print "MESA was not built. Skipping test."
        del instance
        
    def test4(self):
        print "Testing basic operations: evolve..."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
            self.assertEqual(status,0)
            (index_of_the_star, error) = instance.new_particle(1.0)
            self.assertEquals(0, error)
            self.assertEqual(index_of_the_star,1)
            (time_step, error) = instance.get_time_step(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(time_step,1.0e5,1)
            error = instance.evolve(index_of_the_star)
            self.assertEquals(0, error)
            end_time = 5.0e5 # (years)
            instance.evolve_to(index_of_the_star, end_time)
            (age_of_the_star, error) = instance.get_age(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(age_of_the_star,end_time,3)
            (L_of_the_star, error) = instance.get_luminosity(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(L_of_the_star,0.725,3)
            (M_of_the_star, error) = instance.get_mass(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(M_of_the_star,1.000,3)
            (T_of_the_star, error) = instance.get_temperature(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(T_of_the_star,5650.998,3)
            (time_step, error) = instance.get_time_step(index_of_the_star)
            self.assertEquals(0, error)
            self.assertAlmostEqual(time_step,163200.0,1)
        else:
            print "MESA was not built. Skipping test."
        del instance
            
    def xtest5(self):
        print "Testing new ZAMS model explicitly..."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
            self.assertEqual(status,0)
            status = instance.new_zams_model()
            self.assertEqual(status,0)
        else:
            print "MESA was not built. Skipping test."
        del instance     
    
    def xtest6(self):
        print "Testing new ZAMS model implicitly..."
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESAInterface()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
            self.assertEqual(status,0)
            metallicities = [0.00, 0.01, 0.02, 0.04]
            luminosities = [1.717, 0.938, 0.725, 0.592]
            for (i,(Z,L)) in enumerate(zip(metallicities, luminosities)):
                status = instance.set_metallicity(Z)
                self.assertEquals(0, status)
                (index_of_the_star, status) = instance.new_particle(1.0)
                self.assertEquals(0, status)
                self.assertEqual(index_of_the_star,i+1)
                instance.evolve_to(index_of_the_star, 5.0e5)
                (L_of_the_star, status) = instance.get_luminosity(index_of_the_star)
                self.assertEquals(0, status)
                self.assertAlmostEqual(L_of_the_star,L,2)
        else:
            print "MESA was not built. Skipping test."
        del instance     
    
    def test7(self):
        print "Testing MESA stop conditions..."
        instance = MESAInterface()
        if instance.MESA_exists:
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
        else:
            print "MESA was not built. Skipping test."
        del instance
    
    def test8(self):
        print "Testing MESA parameters..."
        instance = MESAInterface()
        if instance.MESA_exists:
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
        else:
            print "MESA was not built. Skipping test."
        del instance


class TestMESA(TestWithMPI):
    
    def test1(self):
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESA()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            status = instance.initialize(instance.default_path_to_inlist, 
                instance.default_path_to_MESA_data)
            self.assertEqual(status,0)
            instance.parameters.set_defaults()
            self.assertEquals(0.02 | units.no_unit , instance.parameters.metallicity)
            instance.parameters.metallicity = 0.01 | units.no_unit
            self.assertEquals(0.01| units.no_unit , instance.parameters.metallicity)
        else:
            print "MESA was not built. Skipping test."
        del instance
    
    def test2(self):
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = MESA()
        #channel.MessageChannel.DEBUGGER = None
        if instance.MESA_exists:
            instance.initialize_module_with_default_parameters()
            stars = core.Stars(1)
            mass = 1.0 | units.MSun # Getting weird errors only (?) if mass == 10
            stars[0].mass = mass
            instance.setup_particles(stars)
            instance.initialize_stars()
            from_code_to_model = instance.particles.new_channel_to(stars)
            from_code_to_model.copy()
            print instance.particles._get_attributes()
#        instance.evolve_model(end_time = 0.003 | units.Myr)
            print instance.particles.time_step
            print instance.get_luminosity(1)
            print instance.particles.stellar_type
            print stars
            self.assertEquals(stars[0].mass, mass)
            #self.assertAlmostEquals(stars[0].luminosity, 5751. | units.LSun, 0)
        else:
            print "MESA was not built. Skipping test."
        del instance
    
    def xtest3(self):
        instance = MESA()
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 5.0 | units.MSun
        star.radius = 0.0 | units.RSun
        
        instance.particles.add_particles(stars)
        instance.initialize_stars()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        previous_type = star.stellar_type
        results = []
        current_time = 0 | units.Myr
        
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
        
        del instance
        
    def xtest4(self):
#       Test whether a set of stars evolve synchronously
#       Create an array of stars with a range in stellar mass
#        masses = [.5, 1., 2., 5., 10., 30.] | units.MSun
        # no high mass stars for now.. (problems with 2nd Asymp. Giant Branch)
        masses = [.5, 1., 1.5] | units.MSun
        max_age = 12 | units.Myr

        number_of_stars=len(masses)
        stars =  core.Stars(number_of_stars)
        for i, star in enumerate(stars):
            star.mass = masses[i]
            star.radius = 0.0 | units.RSun

#       Initialize stellar evolution code
        instance = EVtwin()
        instance.initialize_module_with_default_parameters() 
        if instance.get_maximum_number_of_stars() < number_of_stars:
            instance.set_maximum_number_of_stars(number_of_stars)
        self.assertEqual(instance.parameters.max_age_stop_condition, 2e6 | units.Myr)
        instance.parameters.max_age_stop_condition = max_age
        self.assertEqual(instance.parameters.max_age_stop_condition, max_age)
        instance.setup_particles(stars)
#       Let the code perform initialization actions after all particles have been created. 
        instance.initialize_stars()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        instance.evolve_model(end_time = 10 | units.Myr)
        from_code_to_model.copy()
        
        for i in range(number_of_stars):
            self.assertTrue(stars[i].age.value_in(units.Myr) > 10)
            self.assertTrue(stars[i].age < max_age)
            self.assertTrue(stars[i].mass < masses[i])
            self.assertTrue(stars[i].time_step < max_age)
                
        try:
            instance.evolve_model(end_time = 2*max_age)
            self.fail("Should not be able to evolve beyond maximum age.")
        except Exception as ex:
            self.assertEquals("  2 -- BACKUP -- tstep reduced below limit; quit", str(ex))

        del instance
        

