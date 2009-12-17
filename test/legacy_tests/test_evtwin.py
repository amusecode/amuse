from legacy_support import TestWithMPI
import sys
import os.path

from amuse.legacy.evtwin.interface import EVtwin, EVtwinBinding

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

class TestInterface(TestWithMPI):
    
    def test1(self):
        instance = EVtwin()
        
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
        instance = EVtwin()
        
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
        instance = EVtwin()
        dir = os.path.dirname(sys.modules[instance.__module__].__file__)
        path_to_ev_database = os.path.join(dir, 'src')
        error = instance.set_ev_path(path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        del instance     
        
    
    def test4(self):
        instance = EVtwin()
        
        error = instance.set_ev_path(instance.default_path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        del instance     
    
    def test5(self):
        instance = EVtwin()
        
        error = instance.set_ev_path(instance.default_path_to_ev_database)
        self.assertEquals(0, error)      
        
        error = instance.initialize_code()
        self.assertEquals(0, error)      
        
        (index_of_the_star, error) = instance.new_particle(10)
        self.assertEquals(0, error)       
        
        self.assertTrue(index_of_the_star >= 0)      
        
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)      
        self.assertEquals(10, mass)    
        
        error = instance.evolve(index_of_the_star)
        self.assertEquals(0, error)      
          
        
        (mass, error) = instance.get_mass(index_of_the_star)
        self.assertEquals(0, error)      
        self.assertEquals(10, mass)  
        
        
        (age, error) = instance.get_age(index_of_the_star)
        self.assertEquals(0, error)      
        self.assertTrue(age > 0)      
        
        del instance   
        
class TestInterfaceBinding(TestWithMPI):
    
    class EVtwinWithBinding(EVtwin, EVtwinBinding):
        """
        """
        def __init__(self):
            EVtwin.__init__(self)
            EVtwinBinding.__init__(self)
            
    
    def test1(self):
        instance = self.EVtwinWithBinding()
        
        instance.parameters.set_defaults()
        
        self.assertEquals(10.0 | units.no_unit , instance.parameters.maximum_number_of_stars)
        instance.parameters.maximum_number_of_stars = 12 | units.no_unit
        self.assertEquals(12.0 | units.no_unit , instance.parameters.maximum_number_of_stars)
        del instance
    
    def test2(self):
        instance = self.EVtwinWithBinding()
        instance.parameters.set_defaults()
        instance.set_ev_path(instance.default_path_to_ev_database)
        
        path = os.path.join(instance.default_path_to_ev_database, 'run')
        path = os.path.join(path, 'muse')
        
        #instance.set_init_dat_name(os.path.join(path,'init.dat'))
        #instance.set_init_run_name(os.path.join(path,'init.run'))
        
        stars = core.Stars(1)
        stars[0].mass = 10 | units.MSun
        
        instance.initialize_code()
        
        instance.setup_particles(stars)
        instance.initialize_stars()
        
        #instance.evolve_particles(stars, 2 | units.Myr)
        instance.update_particles(stars)
        
        self.assertEquals(stars[0].mass, 10 | units.MSun)
        self.assertAlmostEquals(stars[0].luminosity.value_in(units.LSun), 5695.19757302 , 6)
        
    
        
        
        
