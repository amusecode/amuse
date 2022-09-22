from amuse.test import amusetest
from amuse.units import constants
from amuse.units.nbody_system import *

class TestNbodyUnits(amusetest.TestCase):
    def test1(self):
        convert_nbody = nbody_to_si(1 | units.parsec, 20 |units.MSun)
        y = 1 | mass
        self.assertEqual(str(y), '1 mass')
        y_in_si = convert_nbody.to_si(y)
        y_in_msun = y_in_si.as_quantity_in(units.MSun)
        self.assertEqual(str(y_in_msun), '20.0 MSun')
        y_in_nbody = convert_nbody.to_nbody(y_in_msun) 
        self.assertEqual(str(y_in_nbody), '1.0 mass')
       
    def test2(self):
        convert_nbody = nbody_to_si(1 | units.MSun, 149.5e6 | units.km)
        y = 29800 | units.m / units.s
        y_in_nbody = convert_nbody.to_nbody(y) 
        self.assertEqual(str(y_in_nbody.unit), 'length * time**-1')
        self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
       
    def test3(self):
        convert_nbody = nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        y = 1 | length * (time**-1)
        y_in_si = convert_nbody.to_si(y) 
        #self.assertEqual(str(y_in_nbody.unit), 'length * time**-1')
        self.assertAlmostEqual(y_in_si.number, 29795.4, -1)
       
    def test4(self):
        convert_nbody = nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        y_in_nbody = convert_nbody.to_nbody(constants.G) 
        self.assertEqual(str(y_in_nbody.unit), 'length**3 * mass**-1 * time**-2' )
        self.assertAlmostEqual(y_in_nbody.number, 1.0, 9)
       
    def test5(self):
        convert_nbody = nbody_to_si(149.5e6 | units.km, 1 | units.MSun)
        y = 29800 | units.m / units.s
        y_in_nbody = convert_nbody.to_nbody(y) 
        self.assertEqual(str(y_in_nbody.unit), 'length * time**-1')
        self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
       

    def test6(self):
        convert_nbody = nbody_to_si(10 | units.kg, 5 | units.m / units. s )
        y = 5 | units.m / units.s
        y_in_nbody = convert_nbody.to_nbody(y) 
        self.assertEqual(str(y_in_nbody.unit), 'length * time**-1')
        self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
        y_in_si = convert_nbody.to_si(y_in_nbody) 
        self.assertAlmostEqual(y_in_si.number, 5.0, 3)
       
    def test7(self):
        convert_nbody = nbody_to_si(1 | units.kg, 1 | units.m / units. s )
        y = 1 | time
        y_in_si = convert_nbody.to_si(y) 
        self.assertEqual(str(y_in_si.unit), 's')
        self.assertAlmostEqual(y_in_si.number, 6.6730000000000296e-11, 3)
       
    def test8(self):
        self.assertTrue(is_nbody_unit(time / length))
        self.assertFalse(is_nbody_unit(units.s / units.m))
        
        
    def test9(self):
        convert_nbody = nbody_to_si(1 | units.kg, 1 | units.m / units. s )
        y = 1.0 | units.none
        y_in_nbody = convert_nbody.to_nbody(y) 
        y_in_si = convert_nbody.to_si(y) 
        self.assertEqual(y_in_nbody, 1.0 | units.none)
        self.assertEqual(y_in_si, 1.0 | units.none)
    
    
        
        
    def test10(self):
        self.assertEqual((time / length).to_array_of_floats(), [1,2,-1, 1, 0, 0, 0, 0, 0])
        self.assertEqual((time / (length * mass ** 2)).to_array_of_floats(), [1,2,-1, 1, -2, 0, 0, 0, 0])
        
        

