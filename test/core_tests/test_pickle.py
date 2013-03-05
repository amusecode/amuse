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

from amuse.datamodel import Particles

import subprocess
import pickle
import sys
import os

class TestPicklingOfUnitsAndQuantities(amusetest.TestCase):

    def test1(self):
        km = 1000 * m
        self.assertEqual(1000, km.value_in(m))
        pickled_km = pickle.dumps(km)
        print pickled_km
        unpickled_km = pickle.loads(pickled_km)
        self.assertEqual(1000, unpickled_km.value_in(m))


    def test2(self):
        km = 1000 * m
        quantity = 12.0 | km
        self.assertEqual(12000, quantity.value_in(m))
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(12000, unpickled_quantity.value_in(m))
        self.assertEqual(quantity, unpickled_quantity)
    
    

    def test3(self):
        print si.system
        pickled_si_sytem = pickle.dumps(si.system)
        unpickled_si_sytem = pickle.loads(pickled_si_sytem)
        self.assertTrue(unpickled_si_sytem is si.system)
    
    

    def test4(self):
        quantity = 12.0 | nbody_system.energy
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
    
    

    def test5(self):
        quantity = 12.0 | parsec
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
        self.assertEqual(12.0, unpickled_quantity.value_in(parsec))
    
    

    def test6(self):
        quantity = [12.0, 15.0] | parsec
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
    
    

    def test7(self):
        quantity = zero
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertTrue(quantity is unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))
    
    def test8(self):
        quantity = 1 | nbody_system.time
        pickled_quantity = pickle.dumps(quantity)
        print pickled_quantity
        unpickled_quantity = pickle.loads(pickled_quantity)
        self.assertEqual(quantity, unpickled_quantity)
        self.assertEqual(str(quantity), str(unpickled_quantity))

    def test9(self):
        quantity = 1.3 | nbody_system.time
        with open("test9.pickle", "w") as stream: 
            pickle.dump(quantity, stream)
        
        pythonpath = os.pathsep.join(sys.path)
        env = os.environ.copy()
        env['PYTHONPATH'] = pythonpath
        process = subprocess.Popen([
                sys.executable,
                "-c",
                "import pickle;stream = open('test9.pickle', 'r'); print str(pickle.load(stream));stream.close()"
            ]
            , stdout=subprocess.PIPE
            , stderr=subprocess.PIPE
            ,env = env
        )
        unpickled_quantity_string, error_string = process.communicate()
        self.assertEqual(process.returncode, 0)        
        self.assertEqual(str(quantity),unpickled_quantity_string.strip())
        
        
    def test10(self):
        quantity = 1  | parsec
        with open("test10.pickle", "w") as stream: 
            pickle.dump(quantity, stream)
               
        pythonpath = os.pathsep.join(sys.path)
        env = os.environ.copy()
        env['PYTHONPATH'] = pythonpath
        process = subprocess.Popen([
                sys.executable,
                "-c",
                "import pickle; stream = open('test10.pickle', 'r'); print str(pickle.load(stream));stream.close()"
            ]
            , stdout=subprocess.PIPE
            , stderr=subprocess.PIPE
            , env = env
        )
        unpickled_quantity_string, error_string = process.communicate()
        self.assertEqual(process.returncode, 0)        
        self.assertEqual(str(quantity),unpickled_quantity_string.strip())
        
class TestPicklingOfParticleSets(amusetest.TestCase):

    def test1(self):
        particles = Particles(4)
        particles.mass = [1,2,3,4] | km
        pickled_particles = pickle.dumps(particles)
        unpickled_particles = pickle.loads(pickled_particles)
        self.assertAlmostRelativeEquals(unpickled_particles.mass, [1,2,3,4] | km)
