from legacy_support import TestWithMPI
import sys
import os.path

from amuse.legacy.evtwin.interface import EVtwin, EVtwinInterface

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy.support import channel

class TestInterface(TestWithMPI):
    
    def test1(self):
        instance = EVtwinInterface()
        
        (metallicity, error) = instance.get_metallicity()
        self.assertEquals(0, error)
        
        for x in [0.1, 0.01, 0.001, 0.0001]:
            error = instance.set_metallicity(x)
            self.assertEquals(0, error)      
            
            (metallicity, error) = instance.get_metallicity()
            self.assertEquals(0, error)      
            self.assertEquals(x, metallicity)      
        
        del instance
        
    
    def test2(self):
        instance = EVtwinInterface()
        
        (value, error) = instance.get_maximum_number_of_stars()
        self.assertEquals(0, error)      
        
        for x in range(10,100,10):
            error = instance.set_maximum_number_of_stars(x)
            self.assertEquals(0, error)      
            
            (value, error) = instance.get_maximum_number_of_stars()
            self.assertEquals(0, error)      
            self.assertEquals(x, value)      
        
        del instance

    
    
    def test3(self):
        instance = EVtwinInterface()
        dir = os.path.dirname(sys.modules[instance.__module__].__file__)
        path_to_ev_database = os.path.join(dir, 'src')
        error = instance.set_ev_path(path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        del instance     
        
    
    def test4(self):
        instance = EVtwinInterface()
        
        error = instance.set_ev_path(instance.default_path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        del instance     
    
    def test5(self):
        #code/library_v2.f:602
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = EVtwinInterface()
        #channel.MessageChannel.DEBUGGER = None
        
        error = instance.set_ev_path(instance.default_path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        (index_of_the_star, error) = instance.new_particle(1.05)
        self.assertEquals(0, error)       
        
        self.assertTrue(index_of_the_star >= 0)      
        
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)      
        self.assertEquals(1.05, mass)    
        
        error = instance.evolve(index_of_the_star)
        self.assertEquals(0, error)      
          
        for i in range(2):
            error = instance.evolve(index_of_the_star)
            self.assertEquals(0, error)      
    
        
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)      
        self.assertTrue(mass < 1.05)  
        
        
    
        (age, error) = instance.get_age(index_of_the_star)
        self.assertEquals(0, error) 
        self.assertTrue(age > 0)      
        
        del instance   
        
class TestInterfaceBinding(TestWithMPI):
    
            
    
    def test1(self):
        instance = EVtwin()
        
        instance.parameters.set_defaults()
        
        self.assertEquals(10.0 | units.no_unit , instance.parameters.maximum_number_of_stars)
        instance.parameters.maximum_number_of_stars = 12 | units.no_unit
        self.assertEquals(12.0 | units.no_unit , instance.parameters.maximum_number_of_stars)
        del instance
    
    def test2(self):
        instance = EVtwin()
        instance.initialize_module_with_default_parameters()
        
        path = os.path.join(instance.default_path_to_ev_database, 'run')
        path = os.path.join(path, 'muse')
        
        #instance.set_init_dat_name(os.path.join(path,'init.dat'))
        #instance.set_init_run_name(os.path.join(path,'init.run'))
        
        stars = core.Stars(1)
        stars[0].mass = 10 | units.MSun
        
        instance.setup_particles(stars)
        instance.initialize_stars()
        print instance.particles[0]
        
        #instance.evolve_particles(stars, 2 | units.Myr)
        instance.update_particles(stars)
        
        self.assertEquals(stars[0].mass, 10 | units.MSun)
        self.assertAlmostEquals(stars[0].luminosity.value_in(units.LSun), 5695.19757302 , 6)
    
    def xtest3(self):
        channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD
        instance = EVtwin()
        channel.MessageChannel.DEBUGGER = None
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 5.0 | units.MSun
        star.radius = 0.0 | units.RSun
        
        instance.particles.add_particles(stars)
        instance.initialize_stars()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        previous_type = star.type
        results = []
        t0 = 0 | units.Myr
        current_time = t0
        
        while current_time < (115 | units.Myr):
            instance.evolve_model()
            from_code_to_model.copy()
            
            current_time = star.age
            print (star.age, star.mass, star.type)
            if not star.type == previous_type:
                results.append((star.age, star.mass, star.type))
                previous_type = star.type
        
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
            
        
        
        
