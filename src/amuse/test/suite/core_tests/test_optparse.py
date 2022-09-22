from amuse.test import amusetest

import numpy
import sys


from amuse.support.exceptions import AmuseException
from amuse.units.quantities import *
from amuse.units import si
from amuse.units import units
from amuse.units import nbody_system

from amuse import datamodel
from amuse.units import optparse

class TestQuantities(amusetest.TestCase):

    def test1(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, dest = "mass", type = float)
        options, args = x.parse_args(['-m', '2.0'])
        self.assertAlmostRelativeEquals(options.mass, 2.0 | units.MSun)
        
    def test2(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, default = 1 | units.MSun, dest = "mass", type = float)
        options, args = x.parse_args(['bla'])
        self.assertAlmostRelativeEquals(options.mass, 1.0 | units.MSun)
        self.assertEqual(args[0], 'bla')
        
    def test3(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, default = '1.5', dest = "mass", type = float)
        options, args = x.parse_args(['bla'])
        self.assertAlmostRelativeEquals(options.mass, 1.5 | units.MSun)
        self.assertEqual(args[0], 'bla')
        
    def test4(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, help = "(unit: %unit)", default = '1.5', dest = "mass", type = float)
        helpstr = x.format_help()
        
        print(helpstr)
        self.assertTrue('(unit: MSun)' in helpstr)
        
        x = optparse.OptionParser()
        x.add_option('-m', unit = nbody_system.mass, help = "(unit: %unit)", default = '1.5', dest = "mass", type = float)
        helpstr = x.format_help()
        
        print(helpstr)
        self.assertTrue('(unit: mass)' in helpstr)
        
    def test5(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, default = 1.5, dest = "mass", type = float)
        options, args = x.parse_args(['bla'])
        self.assertAlmostRelativeEquals(options.mass, 1.5 | units.MSun)
        self.assertEqual(args[0], 'bla')
        
    def test6(self):
        x = optparse.OptionParser()
        x.add_option('-m', unit = units.MSun, help = "(default: %default, unit: %unit)", default = '1.5', dest = "mass", type = float)
        helpstr = x.format_help()
        
        print(helpstr)
        self.assertTrue('unit: MSun)' in helpstr)
        self.assertTrue('(default: 1' in helpstr)
        
