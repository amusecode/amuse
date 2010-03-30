from amuse.support import io
from amuse.support.io import base
from amuse.test import amusetest
from amuse.support.units import nbody_system, units
from amuse.support.data import core

import os

class TestFileFormatProcessor(base.FileFormatProcessor):
    """
    Save files in a test format
    
    long description    
    """
    provided_formats = ['test','123']
    instance = None
    
    def __init__(self, filename = None, set=None, format=None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        TestFileFormatProcessor.instance = self
        self.stored = False
        
    @base.format_option
    def add_comma(self):
        """if True will add a comma between each value"""
        return True
    
    @base.format_option
    def save_fast(self):
        """if True will save faster but less accurate"""
        return False

    def store(self):
        self.stored = True
        
    def load(self):
        return 10
        
class FrameworkTests(amusetest.TestCase):
    
    def test1(self):
        options = TestFileFormatProcessor.get_options()
        self.assertTrue('add_comma' in options)
        self.assertTrue('save_fast' in options)
        TestFileFormatProcessor.register()
        base.write_set_to_file("test.txt", None, format="123")
        self.assertTrue(TestFileFormatProcessor.instance.stored)
        self.assertTrue(TestFileFormatProcessor.instance.add_comma)
        
    def test2(self):
        TestFileFormatProcessor.register()
        base.write_set_to_file("test.txt", None, format="123", add_comma = False)
        self.assertFalse(TestFileFormatProcessor.instance.add_comma)
        
    def test3(self):
        TestFileFormatProcessor.register()
        documentation =  base.write_set_to_file.__doc__
        print documentation
        self.assertTrue("**123**,\n      Save files in a test format" in documentation)
    
    def test4(self):
        options = base.get_options_for_format('123')
        name, description, default = options[0]
        self.assertEquals(name, 'add_comma')
        self.assertEquals(description, 'if True will add a comma between each value')
        self.assertEquals(default, True)
        name, description, default = options[1]
        self.assertEquals(name, 'save_fast')
        self.assertEquals(description, 'if True will save faster but less accurate')
        self.assertEquals(default, False)
        

class FormatTests(amusetest.TestCase):
    
    def test1(self):
        x = core.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.radius = [3.0, 4.0] | nbody_system.length
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        io.write_set_to_file(x, "test.tsf","tsf")
        y = io.read_set_from_file("test.tsf","tsf")
        
        self.assertEquals(x[0].mass, y[0].mass)

        os.remove("test.tsf")
        
    
    def test2(self):
        x = core.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.radius = [3.0, 4.0] | nbody_system.length
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        io.write_set_to_file(x, "test.dyn","dyn")
        y = io.read_set_from_file("test.dyn","dyn")
        
        self.assertEquals(x[0].mass, y[0].mass)
        
        os.remove("test.dyn")
        
    def test3(self):
        x = core.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.radius = [3.0, 4.0] | units.m
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        io.write_set_to_file(x, "test_unit.dyn","dyn", nbody_to_si_converter = convert)
        y = io.read_set_from_file("test_unit.dyn","dyn", nbody_to_si_converter = convert)
        
        self.assertAlmostEquals(x[0].mass, y[0].mass, 8)
        self.assertAlmostEquals(x[1].position.x, y[1].position.x,8)
        self.assertAlmostEquals(x[1].position.y, y[1].position.y,8)
        self.assertAlmostEquals(x[1].position.z, y[1].position.z,8)
        self.assertAlmostEquals(x[0].velocity.y, y[0].velocity.y,8)
        
        
        os.remove("test_unit.dyn")
        
    
    def test4(self):
        
        x = core.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.radius = [3.0, 4.0] | units.m
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        io.write_set_to_file(x, "test_unit.tsf","tsf", nbody_to_si_converter = convert)
        y = io.read_set_from_file("test_unit.tsf","tsf", nbody_to_si_converter = convert)
        
        self.assertAlmostEquals(x[0].mass, y[0].mass, 8)
        self.assertAlmostEquals(x[1].position.x, y[1].position.x,8)
        self.assertAlmostEquals(x[1].position.y, y[1].position.y,8)
        self.assertAlmostEquals(x[1].position.z, y[1].position.z,8)
        self.assertAlmostEquals(x[0].velocity.y, y[0].velocity.y,8)
        
        
        os.remove("test_unit.tsf")
        
    
    def test5(self):
        options = base.get_options_for_format('tsf')
        name, description, default = options[0]
        self.assertEquals(name, 'nbody_to_si_converter')
        self.assertEquals(description, 'tsf datafiles store nbody data, provide a converter to store si data (None means no converter)')
        self.assertEquals(default, None)
        options = base.get_options_for_format('dyn')
        name, description, default = options[0]
        self.assertEquals(name, 'nbody_to_si_converter')

        
        
