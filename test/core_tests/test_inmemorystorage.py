import numpy
import time

from amuse.test import amusetest
from amuse.units import units
from amuse.datamodel.memory_storage import InMemoryGridAttributeStorage
from amuse.datamodel.memory_storage import InMemoryVectorQuantityAttribute
        
class TestInMemoryGridAttributeStorage(amusetest.TestCase):
    
    def test1(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x.set_values_in_store(
            (i,j,k), 
            ['a','b'], 
            [2.0 | units.kg, 1.0 | units.m]
        )
        
        b, a = x.get_values_in_store(None, ['b','a'])
        print b.shape, a.shape
        self.assertEquals(b[0][1][0], 1.0 | units.m)
        self.assertEquals(b[0][0][0], 0.0 | units.m)
        self.assertEquals(a[0][1][0], 2.0 | units.kg)
        self.assertEquals(a[1][3][2], 2.0 | units.kg)
        self.assertEquals(a[1][2][2], 0.0 | units.kg)
                
        (b,) = x.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0 | units.kg)
        self.assertEquals(b[0][0][2], 0.0 | units.kg)
        self.assertEquals(b[1][2][2], 2.0 | units.kg)
        self.assertEquals(b[2][0][0], 2.0 | units.kg)
        self.assertEquals(b[3][2][2], 2.0 | units.kg)
        
        self.assertEquals(b.sum(), 8.0 | units.kg)
        
        self.assertEquals(sorted(x.get_defined_attribute_names()), ["a", "b"])
    
    def test2(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x.set_values_in_store(
            (i,j,k), 
            ['a','b'], 
            [2.0 , 1.0]
        )
        
        b, a = x.get_values_in_store(None, ['b','a'])
        print a
        self.assertEquals(b[0][1][0], 1.0 )
        self.assertEquals(b[0][0][0], 0.0 )
        self.assertEquals(a[0][1][0], 2.0 )
        self.assertEquals(a[1][3][2], 2.0  )
        self.assertEquals(a[1][2][2], 0.0  )
        
        (b,) = x.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0  )
        self.assertEquals(b[0][0][2], 0.0  )
        self.assertEquals(b[1][2][2], 2.0  )
        self.assertEquals(b[2][0][0], 2.0  )
        self.assertEquals(b[3][2][2], 2.0  )
        
        self.assertEquals(b.sum(), 8.0  )
        
        self.assertEquals(sorted(x.get_defined_attribute_names()), ["a", "b"])
        
        y = x.copy()
        
        (b,) = y.get_values_in_store((numpy.s_[0:4], numpy.s_[1:4], numpy.s_[:]), ['a'])
        
        self.assertEquals(b[0][0][0], 2.0  )
        self.assertEquals(b[0][0][2], 0.0  )
        self.assertEquals(b[1][2][2], 2.0  )
        
        
    def test3(self):
        x = InMemoryGridAttributeStorage(5,4,3)
        i = (0,1,2,3,4)
        j = (1,3,1,3,1)
        k = (0,2,0,2,0)
        x.set_values_in_store(
            (i,j,k), 
            ['a','b'], 
            [2.0 | units.kg, 1.0 | units.m]
        )
        
        b, a = x.get_values_in_store((0,1,0), ['b','a'])
        print b, a
        self.assertEquals(b, 1.0 | units.m)
        self.assertEquals(a, 2.0 | units.kg)
        b, a = x.get_values_in_store((0,0,0), ['b','a'])
        print b, a
        self.assertEquals(b, 0.0 | units.m)
        self.assertEquals(a, 0.0 | units.kg)
                
        
class TestInMemoryVectorQuantityAttribute(amusetest.TestCase):
    
    def test1(self):
        quantity  = units.m.new_quantity(numpy.array([1.0,2.0,3.0]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        
        self.assertEquals(attribute.get_length(), 3)
        self.assertEquals(attribute.get_shape(), (3,) )
        self.assertEquals(attribute.get_values([1,2]), [2.0,3.0] | units.m)
    
        attribute.increase_to_length(5)
        self.assertEquals(attribute.get_values(None), [1.0,2.0,3.0,0.0,0.0] | units.m)
        
    def test2(self):
        quantity  = units.m.new_quantity(numpy.array([1.0,2.0,3.0]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        attribute.set_values([1,2], [4.0,5.0] | units.m)
    
        attribute.increase_to_length(5)
        self.assertEquals(attribute.get_values(None), [1.0,4.0,5.0,0.0,0.0] | units.m)

    def test3(self):
        quantity  = units.m.new_quantity(numpy.array([[1.0,2.0,3.0], [4.0,5.0,6.0]]))
        attribute = InMemoryVectorQuantityAttribute('test', quantity.shape, quantity.unit)
        attribute.set_values(None, quantity)
        self.assertEquals( attribute.get_values([1]),  [4.0,5.0,6.0] | units.m)
        self.assertEquals( attribute.get_shape(), (2,3))
        attribute.increase_to_length(4)
        self.assertEquals( attribute.get_shape(), (4,3))
        self.assertEquals(attribute.get_values(None), [[1.0,2.0,3.0], [4.0,5.0,6.0], [0.0,0.0,0.0], [0.0,0.0,0.0]] | units.m)

