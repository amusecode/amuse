from amuse.test import amusetest

from amuse.support.exceptions import AmuseException
from amuse.units import units
from amuse.units import nbody_system
from amuse.units import constants
from amuse.support import console
from amuse.support.console import set_printing_strategy
from amuse.support.console import get_current_printing_strategy

class TestPrintingStrategy(amusetest.TestCase):

    def tearDown(self):
        set_printing_strategy("default")
 
    def test1(self):
        print "Testing get/set of printing strategy"
        self.assertEqual(get_current_printing_strategy(), console.DefaultPrintingStrategy)
        set_printing_strategy("no_unit")
        self.assertEqual(get_current_printing_strategy(), console.NoUnitsPrintingStrategy)
        set_printing_strategy("default")
        self.assertEqual(get_current_printing_strategy(), console.DefaultPrintingStrategy)
        set_printing_strategy("no_units")
        self.assertEqual(get_current_printing_strategy(), console.NoUnitsPrintingStrategy)
        set_printing_strategy("with_units")
        self.assertEqual(get_current_printing_strategy(), console.DefaultPrintingStrategy)
        set_printing_strategy("formal")
        self.assertEqual(get_current_printing_strategy(), console.FormalPrintingStrategy)
        set_printing_strategy("nbody")
        self.assertEqual(get_current_printing_strategy(), console.NBodyPrintingStrategy)
        set_printing_strategy(console.NoUnitsPrintingStrategy)
        self.assertEqual(get_current_printing_strategy(), console.NoUnitsPrintingStrategy)
        set_printing_strategy(console.DefaultPrintingStrategy)
        self.assertEqual(get_current_printing_strategy(), console.DefaultPrintingStrategy)
    
    def test2(self):
        print "Testing no units printing strategy with SI quantities"
        mass     = 1.0 | units.kg
        acc      = 9.8 | units.m / units.s**2
        position = [1, 2.0, 3] | units.m
        mass_in_g = mass.as_quantity_in(units.g)
        pi       = 3.14 | units.none
        
        self.assertEqual(str(mass),      "1.0 kg")
        self.assertEqual(str(acc),       "9.8 m / (s**2)")
        self.assertEqual(str(position),  "[1.0, 2.0, 3.0] m")
        self.assertEqual(str(mass_in_g), "1000.0 g")
        self.assertEqual(str(pi),        "3.14 none")
        
        set_printing_strategy("no_unit")
        self.assertEqual(get_current_printing_strategy(), console.NoUnitsPrintingStrategy)
        self.assertEqual(str(mass),      "1.0")
        self.assertEqual(str(acc),       "9.8")
        self.assertEqual(str(position),  "[1.0, 2.0, 3.0]")
        self.assertEqual(str(mass_in_g), "1000.0")
        self.assertEqual(str(pi),        "3.14")
        set_printing_strategy("default")
    
    def test3(self):
        print "Testing no units printing strategy with N-body quantities"
        mass     = 1.0 | nbody_system.mass
        acc      = 9.8 | nbody_system.acceleration
        position = [1, 2, 3] | nbody_system.length
        
        self.assertEqual(str(mass),      "1.0 mass")
        self.assertEqual(str(acc),       "9.8 length / (time**2)")
        self.assertEqual(str(position),  "[1, 2, 3] length")
        
        set_printing_strategy("no_unit")
        self.assertEqual(str(mass),      "1.0")
        self.assertEqual(str(acc),       "9.8")
        self.assertEqual(str(position),  "[1, 2, 3]")
        set_printing_strategy("default")
    
    def test4(self):
        print "Testing formal printing strategy"
        mass     = 1.0 | units.kg
        acc      = 9.8 | units.m / units.s**2
        position = [1, 2.0, 3] | nbody_system.length
        mass_in_g = mass.as_quantity_in(units.g) * 1.0
        pi       = 3.14 | units.none
        
        set_printing_strategy("formal")
        self.assertEqual(str(mass),      "<quantity 1.0 | kg>")
        self.assertEqual(str(acc),       "<quantity 9.8 | m / (s**2)>")
        self.assertEqual(str(position),  "<quantity [1.0, 2.0, 3.0] | length>")
        self.assertEqual(str(mass_in_g), "<quantity 1000.0 | 0.001 * kg>")#<quantity 1000.0 | g>")
        self.assertEqual(str(pi),        "<quantity 3.14 | none>")
        set_printing_strategy("default")
    
    def test5(self):
        print "Testing nbody printing strategy"
        mass     = 1.0 | nbody_system.mass
        acc      = 9.8 | nbody_system.length / units.s**2
        position = [1, 2, 3] | units.m
        pi       = 3.14 | units.none
        converter = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.m)
        
        set_printing_strategy("nbody")
        self.assertEqual(str(mass), "1.0")
        self.assertRaises(AmuseException, str, acc, expected_message = 
            "Unable to convert length / (s**2) to N-body units. No nbody_converter given")
        self.assertEqual(str(converter.to_nbody(acc * constants.G.number)), "9.8")
        self.assertEqual(str(converter.to_nbody(position)), "[1, 2, 3]")
        self.assertEqual(str(pi), "3.14")
        
        set_printing_strategy("nbody", nbody_converter = converter)
        self.assertEqual(str(mass), "1.0")
        self.assertEqual(str(acc * constants.G.number), "9.8")
        self.assertEqual(str(position), "[1, 2, 3]")
        set_printing_strategy("default")
    
    def test6(self):
        print "Testing astro printing strategy without units printed"
        mass     = 2.0 | 0.5 * units.MSun
        acc      = (0.0098 | nbody_system.length) * (1 | units.Myr**-2).as_quantity_in(units.s**-2)
        position = [0.1, 0.2, 0.3] | nbody_system.length
        energy   = 1e8 | units.erg
        temperature = 5000 | units.K
        pi       = 3.14 | units.none
        converter = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.kpc)
        
        set_printing_strategy("astro", print_units = False)
        self.assertEqual(str(mass), "1.0")
        self.assertRaises(AmuseException, str, acc, expected_message = 
            "Unable to convert length * s**-2 to SI units. No nbody_converter given")
        self.assertEqual(str(converter.to_si(acc)), "9.8")
        self.assertEqual(str(converter.to_si(position)), "[100.0, 200.0, 300.0]")
        self.assertEqual(str(energy), "10.0")
        self.assertEqual(str(constants.G), "0.00449945056135")
        self.assertEqual(str(temperature), "5000")
        self.assertEqual(str(pi), "3.14")
        
        set_printing_strategy("astro", nbody_converter = converter, print_units = False)
        self.assertEqual(str(acc), "9.8")
        set_printing_strategy("default")
    
    def test7(self):
        print "Testing astro printing strategy with units printed"
        mass     = 2.0 | 0.5 * units.MSun
        acc      = (0.0098 | nbody_system.length) * (1 | units.Myr**-2).as_quantity_in(units.s**-2)
        position = [0.1, 0.2, 0.3] | nbody_system.length
        energy   = 1e8 | units.erg
        temperature = 5000 | units.K
        pi       = 3.14 | units.none
        converter = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.kpc)
        
        set_printing_strategy("astro")
        self.assertEqual(str(mass), "1.0 MSun")
        self.assertEqual(str(acc), "0.0098 length * Myr**-2")
        self.assertEqual(str(converter.to_si(acc)), "9.8 parsec * Myr**-2")
        self.assertEqual(str(converter.to_si(position)), "[100.0, 200.0, 300.0] parsec")
        self.assertEqual(str(energy), "10.0 J")
        self.assertEqual(str(constants.G), "0.00449945056135 parsec**3 * MSun**-1 * Myr**-2")
        self.assertEqual(str(temperature), "5000 K")
        self.assertEqual(str(pi), "3.14 none")
        
        set_printing_strategy("astro", nbody_converter = converter)
        self.assertEqual(str(acc), "9.8 parsec * Myr**-2")
        set_printing_strategy("astro", ignore_converter_exceptions = False)
        self.assertRaises(AmuseException, str, acc, expected_message = 
            "Unable to convert length * s**-2 to SI units. No nbody_converter given")
        set_printing_strategy("default")
    
    def test8(self):
        print "Testing SI printing strategy"
        mass     = 2.0 | 0.5 * units.MSun
        acc      = 0.0098 | nbody_system.length / units.Myr**2
        position = [0.1, 0.2, 0.3] | nbody_system.length
        energy   = 1e8 | units.erg
        temperature = 5000 | units.milli(units.K)
        pi       = 3.14 | units.none
        converter = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.kpc)
        
        set_printing_strategy("SI", nbody_converter = converter)
        self.assertEqual(str(mass), "1.98892e+30 kg")
        self.assertEqual(str(acc), "3.03659755643e-10 m * s**-2")
        self.assertEqual(str(position), "[3.08567758128e+18, 6.17135516256e+18, 9.25703274384e+18] m")
        self.assertEqual(str(energy), "10.0 kg * m**2 * s**-2")
        self.assertEqual(str(constants.G), "6.67428e-11 m**3 * kg**-1 * s**-2")
        self.assertEqual(str(temperature), "5.0 K")
        self.assertEqual(str(pi), "3.14 none")
        
        set_printing_strategy("SI", nbody_converter = converter, print_units = False)
        self.assertEqual(str(mass), "1.98892e+30")
        self.assertEqual(str(acc), "3.03659755643e-10")
        self.assertEqual(str(position), "[3.08567758128e+18, 6.17135516256e+18, 9.25703274384e+18]")
        self.assertEqual(str(energy), "10.0")
        self.assertEqual(str(constants.G), "6.67428e-11")
        self.assertEqual(str(temperature), "5.0")
        self.assertEqual(str(pi), "3.14")
        set_printing_strategy("default")
    
    def test9(self):
        print "Testing custom printing strategy"
        mass     = 2.0 | 0.5 * units.MSun
        acc      = (0.0098 | nbody_system.length) * (1 | units.Myr**-2).as_quantity_in(units.s**-2)
        position = [0.1, 0.2, 0.3] | nbody_system.length
        power   = 10 | units.W
        temperature = 5000 | units.K
        pi       = 3.14 | units.none
        converter = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.kpc)
        
        set_printing_strategy("custom", nbody_converter = converter, preferred_units = 
            [units.amu, units.AU, units.minute, units.milli(units.K), units.erg], precision = 3, 
            prefix = "(> ", separator = " <|> ", suffix = " <)")
        self.assertEqual(str(mass), "(> 1.2e+57 <|> amu <)")
        self.assertEqual(str(acc), "(> 7.31e-18 <|> AU * min**-2 <)")
        self.assertEqual(str(position), "(> [2.06e+07, 4.13e+07, 6.19e+07] <|> AU <)")
        self.assertEqual(str(power), "(> 6e+09 <|> erg / min <)")
        self.assertEqual(str(constants.G), "(> 1.19e-67 <|> AU**3 * amu**-1 * min**-2 <)")
        self.assertEqual(str(constants.kB), "(> 1.38e-19 <|> erg * mK**-1 <)")
        self.assertEqual(str(temperature), "(> 5e+06 <|> mK <)")
        self.assertEqual(str(pi), "(> 3.14 <|> none <)")
        set_printing_strategy("default")
    
    def test10(self):
        print "Testing custom printing strategy with precision keyword"
        mass     = 2.0 | 0.5 * units.MSun
        acc      = 0.23456 | 0.54321 * units.m * units.s**-2
        velocity = [-0.12345]*3 | units.km / units.s
        position = [0.1234567890123456789, 0.2, 3.0] | units.AU
        positions = [position.number]*2 | position.unit
        multi_dimensional = [positions.number]*2 | positions.unit
        pi       = 3.1415926535 | units.none
        
        set_printing_strategy("custom", precision = 3)
        self.assertEqual(str(mass), "2 0.5 * MSun")
        self.assertEqual(str(acc), "0.235 0.54321 * m * s**-2")
        self.assertEqual(str(velocity), "[-0.123 -0.123 -0.123] km / s")
        tmp = "[ 0.123  0.2    3.   ]"
        self.assertEqual(str(position), tmp + " AU")
        self.assertEqual(str(positions), "["+tmp+"\n "+tmp+"] AU")
        self.assertEqual(str(multi_dimensional), "[["+tmp+"\n  "+tmp+
            "]\n\n ["+tmp+"\n  "+tmp+"]] AU")
        self.assertEqual(str(pi), "3.14 none")
        set_printing_strategy("default")

