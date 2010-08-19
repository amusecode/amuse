import unittest

import numpy
import os
import inspect

from amuse.support import exceptions
from amuse.support.data.values import Quantity
from amuse.support.units.si import no_unit

class SkipTest(exceptions.AmuseException):
    pass
    
class TestCase(unittest.TestCase):
    PRECISION = int(round(numpy.log10(2.0/(numpy.finfo(numpy.double).eps))))-1
    
    def _convert_to(self, in_units, first, second):
        if in_units:
            return (first.as_quantity_in(in_units), second.as_quantity_in(in_units))
        else:
            return (first, second)
    
    def _check_comparable(self, first, second):
        if isinstance(first, Quantity) is not isinstance(second, Quantity):
            raise TypeError("Cannot compare quantity: {0} with non-quantity: {1}.".format(*(first,second)
                if isinstance(first, Quantity) else (second,first)))
    
    def _check_different(self, first, second):
        if isinstance(first, Quantity):
            different = (first.value_in(second.unit) != second.value_in(second.unit))
        else:
            different = (first != second)
        return numpy.array(different).flatten()
    
    def failUnlessEqual(self, first, second, msg=None, in_units=None):
        self._check_comparable(first, second)
        first, second = self._convert_to(in_units, first, second)
        unequal       = self._check_different(first, second)
        if len(unequal) == 1:
            if unequal[0]:
                raise self.failureException,(msg or '%r != %r' % (first, second))
        elif any(unequal):
            err_list = ["@%i, %r != %r" % (i,first[i], second[i]) for (i,b) in enumerate(unequal) if b]
            err = '\n'.join(err_list)
            raise self.failureException,(msg or err)
    assertEqual = failUnlessEqual
    assertEquals = failUnlessEqual
    
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
        self._check_comparable(first, second)
        first, second = self._convert_to(in_units, first, second)
        delta = second-first  
        
        if isinstance(delta, Quantity):
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
    
    
    def failUnlessAlmostRelativeEqual(self, first, second, places=None, msg=None):
        if places is None:
            places = self.PRECISION
            
        if places <= 0:
            places = self.PRECISION + places
            
        if abs(first) > abs(second):
            relativeError = abs((second-first) / first)
        else:
            relativeError = abs((first-second) / second)
        
        try:   
            maxRelativeError = (relativeError/relativeError) * (0.1 ** places)
        except ZeroDivisionError:
            return
            
        if relativeError >= maxRelativeError:
            raise self.failureException,(msg or "{0!r} != {1!r} within {2!r} places".format(first, second, places))        

           
    assertAlmostRelativeEqual = failUnlessAlmostRelativeEqual
    assertAlmostRelativeEquals = failUnlessAlmostRelativeEqual
    
    def assertIsOfOrder(self, first, second, msg=None):
        ratio = first/second
        if isinstance(ratio, Quantity):
            if ratio.unit.base:
                raise self.failureException,(msg or "Units of {0!r} and {1!r} do not match.".format(first, second))
            ratio = ratio.value_in(no_unit)
        tmp = numpy.array(numpy.round(numpy.log10(ratio)) != 0).flatten()
        if len(tmp) == 1:
            if tmp[0]:
                raise self.failureException,(msg or '%r is not of order %r' % (first, second))                     
        elif any(tmp):
            err_list = ["@%i, %r is not of order %r" % (i,first[i], second[i]) for (i,b) in enumerate(tmp) if b]
            err = '\n'.join(err_list)
            raise self.failureException,(msg or err)                     
    
    def run(self, result=None):
        if result is None:
            result = self.defaultTestResult()
        result.startTest(self)
        testMethod = getattr(self, self._testMethodName)
        try:
            try:
                self.setUp()
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())
                return

            ok = False
            try:
                testMethod()
                ok = True
            except self.failureException:
                result.addFailure(self, self._exc_info())
            except SkipTest, ex:
                ok = True
                pass
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())

            try:
                self.tearDown()
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, self._exc_info())
                ok = False
                
            if ok:
                result.addSuccess(self)
        finally:
            result.stopTest(self)
            
    
    def skip(self, reason):
        try:
            from nose.plugins import skip
            raise skip.SkipTest(reason)
        except ImportError:
            raise SkipTest(reason)
            
    def get_path_to_results(self):
        return get_path_to_results()
    
    def get_amuse_root_dir(self):
        return get_amuse_root_dir()

            
    


class TestWithMPI(TestCase):
    def setUp(self):
        from amuse.support.legacy.core import is_mpd_running
        self.assertTrue(is_mpd_running(), "MPICH2 mpd deamon process not running, cannot run this test as it requires MPI")
            
    def tearDown(self):
        from amuse.support.legacy.core import stop_interfaces
        stop_interfaces()
    
    def new_instance(self, factory, *arguments, **kwarguments):
        try:
            return factory(*arguments, **kwarguments)
        except Exception as message:
            if os.path.exists(os.path.join(os.path.dirname(inspect.getfile(factory)),'src')):
                raise            
            self.skip("Tried to instantiate a new object of the code with type '{0}', but this code is not available".format(factory))
         
            
    def new_instance_of_an_optional_code(self, factory, *arguments, **kwarguments):
        try:
            return factory(*arguments, **kwarguments)
        except Exception as message:
            self.skip("Tried to instantiate a new object of the optional code with type '{0}', but this code is not available".format(factory))
         
         
def get_path_to_results():
    name_of_testresults_directory = 'test_results'
    if os.path.exists(os.path.abspath(name_of_testresults_directory)):
        return os.path.abspath(name_of_testresults_directory)
    
    amuse_root_dir = get_amuse_root_dir()
    test_results_dir = os.path.join(amuse_root_dir, 'test_results')
    if os.path.exists(test_results_dir):
        return test_results_dir
    else:
        return os.getcwd()

def get_amuse_root_dir():
    if 'AMUSE_ROOT_DIR' in os.environ:
        return os.environ['AMUSE_ROOT_DIR']
    result = os.path.abspath(__file__)
    while not os.path.exists(os.path.join(result,'build.py')):
        result = os.path.dirname(result)
    return result
      
