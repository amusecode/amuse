from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system

from amuse.support.data.incode_storage import *


import numpy


class TestGrids(amusetest.TestCase):
    
    def test1(self):
        class Code(object):
            def get_range(self):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k):
                return units.m(i), units.m(j), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
        )
        
        self.assertEquals(storage.storage_shape(), (10, 4, 4))
        self.assertEquals(storage._get_attribute_names(), set(["i","j","k"]))
        
        values = storage._get_values((0,1,1), ("i",))
        self.assertEquals(len(values), 1)
        print values
        self.assertEquals(values[0], 1 | units.m)
        
        values = storage._get_values((0,1,1), ("k","j","i",))
        self.assertEquals(values[0], 4 | units.m)
        self.assertEquals(values[1], 3 | units.m)
        self.assertEquals(values[2], 1 | units.m)
        
    
    def test2(self):
        class Code(object):
            def get_range(self):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k):
                return units.m(i), units.m(j), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
        )
        values = storage._get_values(numpy.s_[0:2], ("i",))
        print values
        self.assertEquals(len(values), 1)
        self.assertEquals(len(values[0]), 2)
        self.assertEquals(values[0].number.shape, (2,4,4))
        self.assertEquals(values[0][0][0][0], 1 | units.m)
        self.assertEquals(values[0][1][0][0], 2 | units.m)
        
    
    def test3(self):
        
        shape = (11,5,5)
            
        class Code(object):
            
            def __init__(self):
                self.storage = numpy.arange(shape[0]*shape[1]*shape[2]).reshape(shape)
                
            def get_range(self):
                return (0,shape[0]-1,0,shape[1]-1,0,shape[2]-1)
                
            def get_a(self,i_s,j_s,k_s):
                return units.m.new_quantity(numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)]))
                
            def set_a(self, i_s, j_s, k_s, values):
                print i_s, j_s, k_s
                print "VALUES:", values
                index = 0
                for i,j,k in zip(i_s, j_s, k_s):
                    self.storage[i][j][k] = values[index].value_in(units.m)
                    index += 1
                    print index
                    
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        values = storage._get_values(None, ("a",))
        print values[0].value_in(units.m)
        self.assertTrue(numpy.all(values[0].value_in(units.m) == code.storage))
        #self.assertTrue(False)
        values = storage._get_values((0,0,0), ("a",))
        self.assertEquals(values[0], 0 | units.m)
        storage._set_values((0,0,0), ("a",), [11.0 | units.m,])
        values = storage._get_values((0,0,0), ("a",))
        self.assertEquals(values[0], 11.0 | units.m)
        values = storage._get_values((0,0), ("a",))
        storage._set_values((0,0), ("a",), [[11.0, 12.0, 13.0, 14.0, 15.0]| units.m,])
        print code.storage[0][0]
        self.assertTrue(numpy.all(code.storage[0][0] == [11.0, 12.0, 13.0, 14.0, 15.0]))
        
