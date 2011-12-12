from amuse.test import amusetest

from amuse.units import constants
from amuse.units.generic_unit_converter import *
class TestGenericUnits(amusetest.TestCase):
    def test1(self):
        L = 1 | length
        T = 1 | time
        V = 1 | speed
        self.assertTrue(L/T == V)

    def test2(self):
        #natural units
        convert_generic = ConvertBetweenGenericAndSiUnits(constants.c, constants.hbar, constants.G, constants.kB)
        M = 1 | mass
        T = 1 | time
        L = 1 | length
        H = 1 | temperature
    
        M_in_si = convert_generic.to_si(M)
        T_in_si = convert_generic.to_si(T)
        L_in_si = convert_generic.to_si(L)
        H_in_si = convert_generic.to_si(H)
    
        self.assertAlmostRelativeEqual(M_in_si, 2.1764411e-8|units.kg, 3)
        self.assertAlmostRelativeEqual(T_in_si, 5.39124e-44|units.s, 3)
        self.assertAlmostRelativeEqual(L_in_si, 1.616252e-35|units.m, 3)
        self.assertAlmostRelativeEqual(H_in_si, 1.416785e32|units.K, 3)

    def test3(self):
        #Gadget units
        UnitLength_in_cm = 3.085678e21 | units.cm# 1.0 kpc
        UnitMass_in_g = 1.989e43 | units.g    # 1.0e10 solar masses
        UnitVelocity_in_cm_per_s = 1e5 | units.cm / units.s# 1 km/sec
        
        convert_generic = ConvertBetweenGenericAndSiUnits(UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s)
        M = 1 | mass
        T = 1 | time
        L = 1 | length
    
        M_in_si = convert_generic.to_si(M)
        T_in_si = convert_generic.to_si(T)
        L_in_si = convert_generic.to_si(L)
    
        self.assertAlmostRelativeEqual(M_in_si, 1.989e40    | units.kg, 3)
        self.assertAlmostRelativeEqual(T_in_si, 3.085678e16 | units.s, 3)
        self.assertAlmostRelativeEqual(L_in_si, 3.085678e19 | units.m, 3)
    
    def test4(self):
        print "Generic units and vector quantities"
        UnitLength_in_cm = 3.085678e21 | units.cm# 1.0 kpc
        UnitMass_in_g = 1.989e43 | units.g    # 1.0e10 solar masses
        UnitVelocity_in_cm_per_s = 1e5 | units.cm / units.s# 1 km/sec
        
        convert_generic = ConvertBetweenGenericAndSiUnits(UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s)
        M = [1, 2] | mass
        T = [1, 2] | time
        L = [1, 2] | length
    
        M_in_si = convert_generic.to_si(M)
        T_in_si = convert_generic.to_si(T)
        L_in_si = convert_generic.to_si(L)
    
        self.assertAlmostEqual(M_in_si, [1.989e40, 2*1.989e40   ] | units.kg, 3, in_units=1.0e10*units.MSun)
        self.assertAlmostEqual(T_in_si, [3.085678e16, 2*3.085678e16] | units.s, 3, in_units=units.s*units.kpc/units.km)
        self.assertAlmostEqual(L_in_si, [3.085678e19, 2*3.085678e19] | units.m, 3, in_units=3.085678e19*units.m)
    

    def test5(self):
       
        self.assertRaises(UnitsNotOrtogonalException, ConvertBetweenGenericAndSiUnits, 1 | units.km, 2 | units.s, 
            3 | units.km / units.s, expected_message = "The number of orthoganal units is incorrect, expected 3 "
            "but found 2. To convert between S.I. units and another system of units a set of quantities with "
            "orthogonal units is needed. These can be quantities with a single unit (such as length or time) "
            "or quantities with a derived units (such as velocity or force)")
            
        
    def test6(self):
        self.assertRaises(NotAQuantityException, ConvertBetweenGenericAndSiUnits, 1)
        self.assertRaises(NotAScalarException, ConvertBetweenGenericAndSiUnits,  1 | units.km, [1,2] | units.s, )
       
    def test7(self):
        converter = ConvertBetweenGenericAndSiUnits(1.0 | units.MSun, 1.0 | units.AU, 1.0 | units.yr)
        G_in_generic = converter.to_generic(constants.G)
        self.assertEqual(str(G_in_generic.unit), 'length**3 * mass**-1 * time**-2')
        self.assertAlmostEqual(G_in_generic.number, 39.4852492465)
        
        time_unit = ((1.0 | units.AU)**3 / (constants.G * (1.0 | units.MSun))).sqrt().as_quantity_in(units.yr)
        self.assertAlmostEqual(time_unit, 0.159141174195 | units.yr)
        converter = ConvertBetweenGenericAndSiUnits(1.0 | units.MSun, 1.0 | units.AU, time_unit)
        G_in_generic = converter.to_generic(constants.G)
        self.assertEqual(str(G_in_generic.unit), 'length**3 * mass**-1 * time**-2')
        self.assertAlmostEqual(G_in_generic.number, 1.0)
    

