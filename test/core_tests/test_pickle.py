from amuse.test import amusetest

import pickle

from amuse.support.exceptions import AmuseException
from amuse.support.units import core
from amuse.support.units.units import *
from amuse.support.units import si
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_system
from amuse.support.units.constants import *
from amuse.support.units.values import zero

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
    
    
