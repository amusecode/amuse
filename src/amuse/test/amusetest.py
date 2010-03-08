import unittest
from amuse.support.data import values
import numpy

class TestCase(unittest.TestCase):
    
    def failUnlessAlmostEqual(self, first, second, places=7, msg=None):                             
        """Fail if the two objects are unequal as determined by their                                
           difference rounded to the given number of decimal places                                  
           (default 7) and comparing to zero.                                                        
                                                                                                     
           Note that decimal places (from zero) are usually not the same                             
           as significant digits (measured from the most signficant digit).                          
                                                                                                     
           The AmuseTestcase differs from the defaulat unittest.Testcase in that                     
           it implements an abs function that can handle quantities.                                 
        """                                                                                          
        delta = second-first  
        
        if isinstance(delta, values.Quantity): 
            absdif = abs(delta.value_in(delta.unit))
            first_nu = first.value_in(first.unit)
            second_nu = second.value_in(first.unit)
        else:             
            absdif = abs(delta) 
               
        tmp = numpy.round(absdif, places) != 0

        if numpy.isscalar(tmp):
            if tmp:
                raise self.failureException,(msg or '%r != %r within %r places' % (first, second, places))                     
        elif any(tmp):
            err_list = ["@%i, %f != %f within %r places" % (i,first_nu[i], second_nu[i], places) for (i,b) in enumerate(tmp) if b]
            err = '\n'.join(err_list)
            raise self.failureException,(msg or err)                     
                                                                                                     
    assertAlmostEqual =  failUnlessAlmostEqual                                                    
    assertAlmostEquals =  failUnlessAlmostEqual
