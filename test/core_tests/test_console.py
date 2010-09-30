from amuse.support.units import units
from amuse.test import amusetest
from amuse.support.data import console
from amuse.support.data.console import set_printing_stategy

class TestPrintingStrategy(amusetest.TestCase):

    def test1(self):
        self.assertEqual(console.current_printing_stategy.__class__, console.DefaultPrintingStrategy)
        set_printing_stategy("no_unit")
        self.assertEqual(console.current_printing_stategy.__class__, console.NoUnitsPrintingStrategy)
        set_printing_stategy("default")
        self.assertEqual(console.current_printing_stategy.__class__, console.DefaultPrintingStrategy)
        set_printing_stategy("no_units")
        self.assertEqual(console.current_printing_stategy.__class__, console.NoUnitsPrintingStrategy)
        set_printing_stategy("with_units")
        self.assertEqual(console.current_printing_stategy.__class__, console.DefaultPrintingStrategy)
        set_printing_stategy(console.NoUnitsPrintingStrategy)
        self.assertEqual(console.current_printing_stategy.__class__, console.NoUnitsPrintingStrategy)
        set_printing_stategy(console.DefaultPrintingStrategy)
        self.assertEqual(console.current_printing_stategy.__class__, console.DefaultPrintingStrategy)
    
    def test2(self):
        mass     = 1.0 | units.kg
        acc      = 9.8 | units.m / units.s**2
        position = [1, 2.0, 3] | units.m
        mass_in_g = mass.as_quantity_in(units.g)
        
        self.assertEqual(str(mass),      "1.0 kg")
        self.assertEqual(str(acc),       "9.8 m / (s**2)")
        self.assertEqual(str(position),  "[1.0, 2.0, 3.0] m")
        self.assertEqual(str(mass_in_g), "1000.0 g")
        
        set_printing_stategy("no_unit")
        self.assertEqual(console.current_printing_stategy.__class__, console.NoUnitsPrintingStrategy)
        self.assertEqual(str(mass),      "1.0")
        self.assertEqual(str(acc),       "9.8")
        self.assertEqual(str(position),  "[1.0, 2.0, 3.0]")
        self.assertEqual(str(mass_in_g), "1000.0")
        set_printing_stategy("default")
   

