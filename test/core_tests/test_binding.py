import unittest

from amuse.support.units import units
from amuse.support.data import binding

import numpy

class TestBase(unittest.TestCase):
    
    class Test(object):
        
        test_property = binding.CodeProperty("get_property", units.m)
        test_property2 = binding.CodeProperty("property_value", units.m)
        
        def get_property(self):
            return (self.property_value, self.errorcode)
        
    
    def test1(self):
        o = self.Test()
        o.errorcode = 0
        o.property_value = 2.0
        
        self.assertEquals(o.test_property, 2.0 | units.m)
        try:
            o.errorcode = -1
            o.test_property
            self.fail("method returned an errorcode, so the property should raise an exception")
        except Exception, ex:
            self.assertEquals(str(ex), "calling 'get_property' to get the value for property 'test_property' resulted in an error (errorcode -1)")
            
        try:
            o.test_property = 4.0 | units.m
            self.fail("values of properties cannot be set")
        except Exception, ex:
            self.assertEquals(str(ex), "property 'test_property' of a 'Test' object is read-only, you cannot change it's value")
    
    def test2(self):
        o = self.Test()
        o.property_value = 3.0
        
        self.assertEquals(o.test_property2, 3.0 | units.m)
        
        
