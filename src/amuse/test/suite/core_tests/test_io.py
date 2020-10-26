from amuse.test import amusetest
from amuse.support.exceptions import AmuseException

import os
from amuse import io
from amuse.io import base
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
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
        base.write_set_to_file(None, "test.txt", format="123")
        self.assertEqual(TestFileFormatProcessor.instance.set, None)
        self.assertEqual(TestFileFormatProcessor.instance.filename, "test.txt")
        self.assertEqual(TestFileFormatProcessor.instance.format, "123")
        self.assertTrue(TestFileFormatProcessor.instance.stored)
        self.assertTrue(TestFileFormatProcessor.instance.add_comma)
        
    def test2(self):
        TestFileFormatProcessor.register()
        base.write_set_to_file(None, "test.txt", format="123", add_comma = False)
        self.assertFalse(TestFileFormatProcessor.instance.add_comma)
        
    def test3(self):
        TestFileFormatProcessor.register()
        documentation =  base.write_set_to_file.__doc__
        print(documentation)
        self.assertTrue("**123**,\n      Save files in a test format" in documentation)
    
    def test4(self):
        options = base.get_options_for_format('123')
        options.sort(key = lambda x: x[0])
        for x in options:
            print(x)
        name, description, default = options[0]
        self.assertEqual(name, 'add_comma')
        self.assertEqual(description, 'if True will add a comma between each value')
        self.assertEqual(default, True)
        name, description, default = options[2]
        self.assertEqual(name, 'save_fast')
        self.assertEqual(description, 'if True will save faster but less accurate')
        self.assertEqual(default, False)
    
    def test5(self):
        self.assertRaises(AmuseException, io.read_set_from_file, "non_existent","test", 
            expected_message = "IO exception: Error: file 'non_existent' does not exist.")
        
        processor = base.FileFormatProcessor(format="test")
        self.assertRaises(base.CannotSaveException, processor.store, 
            expected_message = "You tried to save a file with fileformat 'test', but"
                " this format is not supported for writing files")
        
        self.assertRaises(base.CannotLoadException, processor.load, 
            expected_message = "You tried to load a file with fileformat 'test', but"
                " this format is not supported for reading files")


class FormatTests(amusetest.TestCase):
    
    def test1(self):
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.radius = [3.0, 4.0] | nbody_system.length
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        io.write_set_to_file(x, "test.tsf","tsf")
        y = io.read_set_from_file("test.tsf","tsf")
        
        self.assertAlmostEqual(x.mass, y.mass, 8)
#        self.assertAlmostEquals(x.radius, y.radius, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)

        os.remove("test.tsf")
        
    
    def test2(self):
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.radius = [3.0, 4.0] | nbody_system.length
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        io.write_set_to_file(x, "test.dyn","dyn")
        y = io.read_set_from_file("test.dyn","dyn")
        
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        self.assertRaises(AttributeError, lambda: y.radius, 
            expected_message = "You tried to access attribute 'radius' but this "
                "attribute is not defined for this set.")
        
        os.remove("test.dyn")
        
    def test3(self):
        x = datamodel.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.radius = [3.0, 4.0] | units.m
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        io.write_set_to_file(x, "test_unit.dyn","dyn", nbody_to_si_converter = convert)
        y = io.read_set_from_file("test_unit.dyn","dyn", nbody_to_si_converter = convert)
        
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        self.assertRaises(AttributeError, lambda: y.radius, 
            expected_message = "You tried to access attribute 'radius' but this "
                "attribute is not defined for this set.")
        
        os.remove("test_unit.dyn")
        
    
    def test4(self):
        
        x = datamodel.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.radius = [3.0, 4.0] | units.m
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        io.write_set_to_file(x, "test_unit.tsf","tsf", nbody_to_si_converter = convert)
        y = io.read_set_from_file("test_unit.tsf","tsf", nbody_to_si_converter = convert)
        
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        self.assertRaises(AttributeError, lambda: y.radius, 
            expected_message = "You tried to access attribute 'radius' but this "
                "attribute is not defined for this set.")
        
        os.remove("test_unit.tsf")
    
    def test5(self):
        print("Testing HDF5 io")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        x.radius = [3.0, 4.0] | units.m
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        io.write_set_to_file(x, "test_unit.hdf5","hdf5")
        y = io.read_set_from_file("test_unit.hdf5","hdf5")
        
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.radius, y.radius, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        
        os.remove("test_unit.hdf5")
        
    def test6(self):
        print("Testing HDF5 io, with options")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="1.0")
        x.mass = [10.0, 20.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="1.0")
        x.mass = [100.0, 200.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="1.0", append_to_file=True)
        y = io.read_set_from_file("test_unit.hdf5","hdf5", copy_history=False, close_file=False)
        y = y.previous_state() # weirdness for version="1.0"
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual([10.0, 20.0] | units.kg, y.previous_state().mass, 8)
        self.assertAlmostEqual([1.0, 2.0] | units.kg, y.previous_state().previous_state().mass, 8)
        self.assertEqual(y.previous_state().previous_state().previous_state(), None)
        
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=False, overwrite_file=True,version="1.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5")
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertEqual(y.previous_state().previous_state(), None)
        
        os.remove("test_unit.hdf5")

    def test6b(self):
        print("Testing HDF5 io, with options")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="1.0")
        x.mass = [10.0, 20.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="1.0")
        x.mass = [100.0, 200.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="1.0", append_to_file=True)
        y = io.read_set_from_file("test_unit.hdf5","hdf5", copy_history=True, )
        y = y.previous_state() # weirdness for version="1.0"
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual([10.0, 20.0] | units.kg, y.previous_state().mass, 8)
        self.assertAlmostEqual([1.0, 2.0] | units.kg, y.previous_state().previous_state().mass, 8)
        self.assertEqual(y.previous_state().previous_state().previous_state(), None)
        
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=False, overwrite_file=True,version="1.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5")
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertEqual(y.previous_state().previous_state(), None)
        
        os.remove("test_unit.hdf5")

    def test6c(self):
        print("Testing HDF5 io, with options, version=2.0")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="2.0")
        x.mass = [10.0, 20.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="2.0")
        x.mass = [100.0, 200.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="2.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5", close_file=False, copy_history=False)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual([10.0, 20.0] | units.kg, y.previous_state().mass, 8)
        self.assertAlmostEqual([1.0, 2.0] | units.kg, y.previous_state().previous_state().mass, 8)
        self.assertEqual(y.previous_state().previous_state().previous_state(), None)
        
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=False, overwrite_file=True, version="2.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5", copy_history=False)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertEqual(y.previous_state(), None)
        
        os.remove("test_unit.hdf5")

    def test6d(self):
        print("Testing HDF5 io, with options, version=2.0")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", version="2.0")
        x.mass = [10.0, 20.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="2.0")
        x.mass = [100.0, 200.0] | units.kg
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=True, version="2.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5", close_file=True, copy_history=True)
        y = y.previous_state() # weirdness for version="1.0"
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual([10.0, 20.0] | units.kg, y.previous_state().mass, 8)
        self.assertAlmostEqual([1.0, 2.0] | units.kg, y.previous_state().previous_state().mass, 8)
        self.assertEqual(y.previous_state().previous_state().previous_state(), None)
        
        io.write_set_to_file(x, "test_unit.hdf5","hdf5", append_to_file=False, overwrite_file=True, version="2.0")
        y = io.read_set_from_file("test_unit.hdf5","hdf5", copy_history=False)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertEqual(y.previous_state(), None)
        
        os.remove("test_unit.hdf5")


    def test7(self):
        print("Testing HDF5 io with a ParticlesSuperset")
        if os.path.exists("test_unit.hdf5"):
            os.remove("test_unit.hdf5")
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(2)
        superset = datamodel.ParticlesSuperset([set1, set2])
        superset.mass = [1.0, 2.0, 3.0, 4.0] | units.kg
        superset.radius = [3.0, 4.0, 5.0, 6.0] | units.m
        superset.position = [[1,2,3], [3,5,6], [3,2,1], [-3,-5,-6]] | units.m
        superset.velocity = [[1,2,3], [3,5,6], [3,2,1], [-3,-5,-6]] | units.m / units.s
        io.write_set_to_file(superset, "test_unit.hdf5","hdf5")
        y = io.read_set_from_file("test_unit.hdf5","hdf5")
        
        self.assertAlmostEqual(superset.mass, y.mass, 8)
        self.assertAlmostEqual(superset.radius, y.radius, 8)
        self.assertAlmostEqual(superset.position, y.position,8)
        self.assertAlmostEqual(superset.velocity, y.velocity,8)
        
        os.remove("test_unit.hdf5")
    
    def test8(self):
        options = base.get_options_for_format('tsf')
        options.sort(key = lambda x: x[0])
        name, description, default = options[1]
        self.assertEqual(name, 'nbody_to_si_converter')
        self.assertEqual(description, 'NEMO datafiles store nbody data, provide a '
            'converter to store si data (None means no converter)')
        self.assertEqual(default, None)
        options = base.get_options_for_format('dyn')
        options.sort(key = lambda x: x[0])
        name, description, default = options[1]
        self.assertEqual(name, 'dynamics_mass_units')
        
        options = base.get_options_for_format('hdf5')
        options.sort(key = lambda x: x[0])
        name, description, default = options[1]
        self.assertEqual(name, 'append_to_file')
        self.assertTrue(description.find('If set to True, new data is appended to HDF5 files.') >= 0)
        self.assertTrue(description.find('If set to False, the existing file is removed and overwritten.') >= 0)
        self.assertEqual(default, False)
    
    def test9(self):
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.kg
        self.assertRaises(AmuseException, io.write_set_to_file, x, "test_unit.bogus", "bogus", 
            expected_message = "You tried to load or save a file with fileformat 'bogus'"
                ", but this format is not in the supported formats list")
    
    def test10(self):
        print("Testing saving/loading timestamp in Starlab")
        x = datamodel.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        current_time = 1.0 | units.Myr
        io.write_set_to_file(x.savepoint(current_time), "time_test_unit.dyn","dyn", nbody_to_si_converter = convert)
        y = io.read_set_from_file("time_test_unit.dyn","dyn", nbody_to_si_converter = convert)
        self.assertAlmostEqual(current_time, y.previous_state().get_timestamp(), 8, in_units=units.Myr)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        current_time = 1.0 | nbody_system.time
        io.write_set_to_file(x.savepoint(current_time), "time_test_unit.dyn","dyn")
        y = io.read_set_from_file("time_test_unit.dyn","dyn")

        self.assertAlmostEqual(current_time, y.previous_state().get_timestamp(), 8, in_units=nbody_system.time)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        
        os.remove("time_test_unit.dyn")
    
    def test11(self):
        print("Testing saving/loading timestamp in NEMO")
        x = datamodel.Particles(2)
        convert = nbody_system.nbody_to_si(1 | units.kg, 2 | units.m)
        x.mass = [1.0, 2.0] | units.kg
        x.position = [[1,2,3], [3,5,6]] | units.m
        x.velocity = [[1,2,3], [3,5,6]] | units.m / units.s
        current_time = 1.0 | units.Myr
        io.write_set_to_file(x.savepoint(current_time), "time_test_unit.tsf","tsf", nbody_to_si_converter = convert)
        y = io.read_set_from_file("time_test_unit.tsf","tsf", nbody_to_si_converter = convert)
        
        self.assertAlmostEqual(current_time, y.previous_state().get_timestamp(), 8, in_units=units.Myr)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | nbody_system.mass
        x.position = [[1,2,3], [3,5,6]] | nbody_system.length
        x.velocity = [[1,2,3], [3,5,6]] | nbody_system.speed
        current_time = 1.0 | nbody_system.time
        io.write_set_to_file(x.savepoint(current_time), "time_test_unit.tsf","tsf")
        y = io.read_set_from_file("time_test_unit.tsf","tsf")

        self.assertAlmostEqual(current_time, y.previous_state().get_timestamp(), 8, in_units=nbody_system.time)
        self.assertAlmostEqual(x.mass, y.mass, 8)
        self.assertAlmostEqual(x.position, y.position,8)
        self.assertAlmostEqual(x.velocity, y.velocity,8)
        
        os.remove("time_test_unit.tsf")
        
    def test12(self):
        all_formats = sorted(base.registered_fileformat_processors.keys())
        for x in all_formats:
            options = base.get_options_for_format('txt')
            self.assertTrue(len(options) >= 0)

