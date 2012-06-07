from amuse.test.amusetest import TestWithMPI

from amuse.community.seba.interface import SeBaInterface, SeBa

from amuse.units import units
from amuse.datamodel import Particle

class TestSeBaInterface(TestWithMPI):

    def test1(self):
        instance = self.new_instance_of_an_optional_code(SeBaInterface)
            
        endtime, mass, radius, luminosity, temperature, time_step, stellar_type, error = instance.evolve_star(1, 4600, 0.02)
        self.assertEquals(error, 0)
        self.assertTrue( endtime <= 4600.0)
        self.assertAlmostRelativeEqual(endtime, 4600.0, 4)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        self.assertAlmostRelativeEqual(radius, 0.9856, 4)
        self.assertAlmostRelativeEqual(luminosity, 0.9585, 4)
        self.assertAlmostRelativeEqual(temperature, 5751, 4)
        self.assertAlmostRelativeEqual(time_step, 1089.3, 4)
        self.assertEqual(stellar_type, 1)
        
        instance.stop()
        
    def test2(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle(1.)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 0.88824945029751212, 6)
    
        stellar_type, error = instance.get_stellar_type(index)
        self.assertEquals(error, 0)
        self.assertEquals(stellar_type, 1)
        
        instance.stop()
        
    def test3(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle(1.)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        error = instance.evolve_model(4600)
        self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 0.9856, 4)
        value, error = instance.get_temperature(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 5751, 4)
        value, error = instance.get_time_step(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 1089.3, 4)
    
        stellar_type, error = instance.get_stellar_type(index)
        self.assertEquals(error, 0)
        self.assertEquals(stellar_type, 1)
        
        instance.stop()
                
    def test4(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle(1.)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        for t in range(46):
            error = instance.evolve_model((t+1) * 100)
            self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 0.9856, 4)
        value, error = instance.get_temperature(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 5751, 4)
        value, error = instance.get_time_step(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(value, 1089.3, 4)
    
        stellar_type, error = instance.get_stellar_type(index)
        self.assertEquals(error, 0)
        self.assertEquals(stellar_type, 1)
        
        instance.stop()      
        
    def test5(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle([1., 2., 3.])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2,3])
        
        mass, error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 2 , 6)
        
        mass, error = instance.get_mass(3)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)
        
        
        error = instance.evolve_model(4600)
        self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(index)
        print mass
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass[0], 1.0, 6)
        self.assertAlmostRelativeEqual(mass[1], 0.63151129, 6)
        self.assertAlmostRelativeEqual(mass[2], 0.74066592, 6)
        
        
        instance.stop()  
        
    def test6(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle([1., 2., 3.])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2,3])
        
        for t in range(46):
            error = instance.evolve_model((t+1) * 100)
            self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(index)
        print mass
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, [1.0, 0.63151, 0.74066], 4)
        
        instance.stop()
    
    def test7(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle([1., 2., 3.])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2,3])
        
        mass, error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 2 , 6)
        
        mass, error = instance.get_mass(3)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)
        
        mass, error = instance.get_mass(4)
        self.assertEquals(error, -1)
        
        error = instance.delete_star(2)
        self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(2)
        self.assertEquals(error, -1)
        
        mass, error = instance.get_mass(3)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)
        
        index, error = instance.new_particle(4.)
        self.assertEquals(error, 0)
        self.assertEquals(index, 4)
        
        instance.stop()
    
    def test8(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        index,error = instance.new_particle([3.0,1.0,2.0])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2,3])
        
        error = instance.delete_star(1)
        self.assertEquals(error, 0)
        
        error = instance.evolve_model(4600);
        
        mass, error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 1, 6)
        
        error = instance.delete_star(3)
        self.assertEquals(error, 0)
        
        
        index,error = instance.new_particle([5.0])
        self.assertEquals(error, 0)
        
        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 5.0, 6)
        error = instance.evolve_model(5000);
        
        
        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.99057, 4)
        
        error = instance.delete_star(2)
        self.assertEquals(error, 0)
        error = instance.delete_star(index)
        self.assertEquals(error, 0)
        
        for i in range(4):
            mass, error = instance.get_mass(index+1)
            self.assertEquals(error, -1)
        
        index,error = instance.new_particle([5.0])
        self.assertEquals(error, 0)
        
        error = instance.evolve_model(5000);

        mass, error = instance.get_mass(index)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.99057, 4)
        
        instance.stop()
    
class TestSeBa(TestWithMPI):

    def test1(self):
        
        instance = self.new_instance_of_an_optional_code(SeBa)
            
        endtime, mass, radius, luminosity, temperature, time_step, stellar_type = instance.evolve_star(1 | units.MSun, 4600 | units.Myr, 0.02)
        
        self.assertTrue( endtime <= 4600 | units.Myr)
        self.assertAlmostRelativeEqual(mass, 1.0 | units.MSun, 4)
        self.assertAlmostRelativeEqual(radius, 0.9856 | units.RSun, 4)
        self.assertAlmostRelativeEqual(luminosity, 0.9585 | units.LSun, 4)
        self.assertAlmostRelativeEqual(temperature, 5751 | units.K, 4)
        self.assertAlmostRelativeEqual(time_step, 1089.3 | units.Myr, 4)
        self.assertEqual(stellar_type, 1 | units.stellar_type)

        
    def test2(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        
        p = Particle()
        p.mass = 5 | units.MSun
        p.metallicity = 0.02
        
        p = instance.particles.add_particle(p)
        instance.evolve_model(130 | units.Myr)

        self.assertAlmostRelativeEqual(p.mass, 0.9906 | units.MSun, 4)
        
