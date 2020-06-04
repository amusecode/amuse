from amuse.test import amusetest

import pickle

from amuse.support.exceptions import AmuseException

from amuse.units import core
from amuse.units import si
from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units.quantities import zero
from amuse.units.units import *
from amuse.units.constants import *

from amuse.datamodel import Particles, parameters

import subprocess
import pickle
import sys
import os

class TestPicklingOfUnitsAndQuantities(amusetest.TestCase):

    def test1(self):
        km = 1000 * m
        self.assertEqual(1000, km.value_in(m))
        pickled_km = pickle.dumps(km)
        print(pickled_km)
        unpickled_km = pickle.loads(pickled_km)
        self.assertEqual(1000, unpickled_km.value_in(m))


    def test2(self):
        km = 1000 * m
        quantity = 12.0 | km
        self.assertEqual(12000, quantity.value_in(m))
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(12000, unpickled_quantity.value_in(m))
        self.assertEqual(quantity, unpickled_quantity)
    
    

    def test3(self):
        print(si.system)
        pickled_si_sytem = pickle.dumps(si.system)
        unpickled_si_sytem = pickle.loads(pickled_si_sytem)
        self.assertTrue(unpickled_si_sytem is si.system)
    
    

    def test4(self):
        quantity = 12.0 | nbody_system.energy
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
    
    

    def test5(self):
        quantity = 12.0 | parsec
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
        self.assertEqual(12.0, unpickled_quantity.value_in(parsec))
    
    

    def test6(self):
        quantity = [12.0, 15.0] | parsec
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
    
    

    def test7(self):
        quantity = zero
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertTrue(quantity is unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
    
    def test8(self):
        quantity = 1 | nbody_system.time
        pickled_quantity = pickle.dumps(quantity)
        print(pickled_quantity)
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))

    def test9(self):
        quantity = 1.3 | nbody_system.time
        path=os.path.abspath(os.path.join(self.get_path_to_results(), "test9.pickle"))

        with open(path, "wb") as stream: 
            pickle.dump(quantity, stream)
        
        pythonpath = os.pathsep.join(sys.path)
        env = os.environ.copy()
        env['PYTHONPATH'] = pythonpath
        code = "import pickle;stream = open('{0}', 'rb'); print str(pickle.load(stream));stream.close()".format(path)
        if sys.hexversion > 0x03000000:
            code = "import pickle;stream = open('{0}', 'rb'); print(str(pickle.load(stream)));stream.close()".format(path)
       
        process = subprocess.Popen([
                sys.executable,
                "-c",
                code
            ]
            , stdout=subprocess.PIPE
            , stderr=subprocess.PIPE
            ,env = env
        )
        unpickled_quantity_string, error_string = process.communicate()
        self.assertEqual(process.returncode, 0)        
        self.assertEqual(str(quantity),unpickled_quantity_string.strip().decode('utf-8'))
        
        
    def test10(self):
        quantity = 1  | parsec
        path=os.path.abspath(os.path.join(self.get_path_to_results(), "test10.pickle"))
        with open(path, "wb") as stream: 
            pickle.dump(quantity, stream)
               
        pythonpath = os.pathsep.join(sys.path)
        env = os.environ.copy()
        env['PYTHONPATH'] = pythonpath
        code = "import pickle;stream = open('{0}', 'rb'); print str(pickle.load(stream));stream.close()".format(path)
        if sys.hexversion > 0x03000000:
            code = "import pickle;stream = open('{0}', 'rb'); print(str(pickle.load(stream)));stream.close()".format(path)
       
        process = subprocess.Popen([
                sys.executable,
                "-c",
               code
            ]
            , stdout=subprocess.PIPE
            , stderr=subprocess.PIPE
            , env = env
        )
        unpickled_quantity_string, error_string = process.communicate()
        self.assertEqual(process.returncode, 0)        
        self.assertEqual(str(quantity),unpickled_quantity_string.strip().decode('utf-8'))
        
    
    def test11(self):
        value = 1 | stellar_type
        print(value)
        self.assertEqual(1 | stellar_type, value)
        pickled = pickle.dumps(value)
        unpickled_value = pickle.loads(pickled)
        self.assertEqual(1 | stellar_type, unpickled_value)

        
class TestPicklingOfParticleSets(amusetest.TestCase):

    def test1(self):
        particles = Particles(4)
        particles.mass = [1,2,3,4] | km
        pickled_particles = pickle.dumps(particles)
        unpickled_particles = pickle.loads(pickled_particles)
        self.assertAlmostRelativeEquals(unpickled_particles.mass, [1,2,3,4] | km)
    
    def test2(self):
        particles = Particles(4)
        particles.mass = [1, 2, 3, 6] | kg
        particles.position = [[0, 0, 0], [3, 0, 0], [0, 4, 0], [3, 4, 0]] | m
        self.assertEqual(particles.center_of_mass(), [2, 3, 0] | m)
        pickled_particles = pickle.dumps(particles)
        unpickled_particles = pickle.loads(pickled_particles)
        self.assertAlmostRelativeEquals(unpickled_particles.mass, [1, 2, 3, 6] | kg)
        self.assertEqual(unpickled_particles.center_of_mass(), [2, 3, 0] | m)
    

    def test3(self):
        particles = Particles(4)
        particles.mass = [1, 2, 3, 6] | kg
        particles.position = [[0, 0, 0], [3, 0, 0], [0, 4, 0], [3, 4, 0]] | m
        self.assertEqual(particles.center_of_mass(), [2, 3, 0] | m)
        pickled_particles = pickle.dumps(particles)
        unpickled_particles = pickle.loads(pickled_particles)
        pickled_particles = pickle.dumps(particles)     #dump it twice! 
        unpickled_particles = pickle.loads(pickled_particles)
        self.assertAlmostRelativeEquals(unpickled_particles.mass, [1, 2, 3, 6] | kg)
        self.assertEqual(unpickled_particles.center_of_mass(), [2, 3, 0] | m)

    def test4(self):
        particles = Particles(4)
        particles.mass = [1, 2, 3, 6] | kg
        particles.position = [[0, 0, 0], [3, 0, 0], [0, 4, 0], [3, 4, 0]] | m
        pickled_particles = pickle.dumps(list(particles))
        print(len(pickled_particles))
        unpickled_particles = pickle.loads(pickled_particles)
        self.assertEqual(len(unpickled_particles) , 4)
        unpickled_particles = Particles(particles=unpickled_particles)
        self.assertAlmostRelativeEquals(unpickled_particles.mass, [1, 2, 3, 6] | kg)
        self.assertEqual(unpickled_particles.center_of_mass(), [2, 3, 0] | m)




class BaseTestModule(object):
    def before_get_parameter(self):
        return
        
    def before_set_parameter(self):
        return
        
class TestPicklingOfParameters(amusetest.TestCase):
    
    def test1(self):
        definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | m
        )
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        set.test_name = 10| m
        
        self.assertEqual(o.x, 10|m)
        self.assertEqual(set.test_name, 10|m)
        
        memento = set.copy()
        self.assertEqual(memento.test_name, 10|m)

        pickled_memento=pickle.dumps(memento)
        unpickled_memento=pickle.loads(pickled_memento)
        
        self.assertEqual(memento.test_name, unpickled_memento.test_name)
        
        
