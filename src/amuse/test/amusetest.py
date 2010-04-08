import unittest
from amuse.support.data import values
import numpy

class TestCase(unittest.TestCase):
    
    def failUnlessAlmostEqual(self, first, second, places=7, msg=None, in_units=None):                             
        """Fail if the two objects are unequal as determined by their                                
           difference rounded to the given number of decimal places                                  
           (default 7) and comparing to zero.                                                        
                                                                                                     
           Note that decimal places (from zero) are usually not the same                             
           as significant digits (measured from the most signficant digit).                          
                                                                                                     
           The AmuseTestcase differs from the defaulat unittest.Testcase in that                     
           it implements an abs function that can handle quantities. 
           
           If the base unit for the test (in_units) is not supplied, the
           difference of the two objects is evaluated in units of [second.unit].
        """                                                                                          
        if in_units:
            first = first.as_quantity_in(in_units)
            second = second.as_quantity_in(in_units)
        delta = second-first  
        
        if isinstance(delta, values.Quantity):
            absdif = abs(delta.value_in(delta.unit))
        else:             
            absdif = abs(delta) 
               
        tmp = numpy.array(numpy.round(absdif, places) != 0).flatten()

        if len(tmp) == 1:
            if tmp[0]:
                raise self.failureException,(msg or '%r != %r within %r places' % (first, second, places))                     
        elif any(tmp):
            err_list = ["@%i, %r != %r within %r places" % (i,first[i], second[i], places) for (i,b) in enumerate(tmp) if b]
            err = '\n'.join(err_list)
            raise self.failureException,(msg or err)                     
                                                                                                     
    assertAlmostEqual =  failUnlessAlmostEqual                                                    
    assertAlmostEquals =  failUnlessAlmostEqual
    
    
    def failUnlessAlmostRelativeEqual(self, first, second, places=7, msg=None):
        
        if abs(first) > abs(second):
            relativeError = abs((second-first) / first)
        else:
            relativeError = abs((first-second) / second)

        maxRelativeError = (relativeError/relativeError) * (0.1 ** places)
        
        if relativeError >= maxRelativeError:
            raise self.failureException,(msg or "{0!r} != {1!r} within {2!r} places".format(first, second, places))        

           
    assertAlmostRelativeEqual = failUnlessAlmostRelativeEqual
    assertAlmostRelativeEquals = failUnlessAlmostRelativeEqual
         
         
        
