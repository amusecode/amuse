import unittest

import numpy
import os
import inspect

from amuse.support.data import values


class SkipTest(Exception):
    pass
    
class TestCase(unittest.TestCase):
    PRECISION = int(round(numpy.log10(2.0/(numpy.finfo(numpy.double).eps))))-1
    
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
        if hasattr(first, 'unit') is not hasattr(second, 'unit'):
            raise TypeError("Cannot compare quantity: {0} with non-quantity: {1}.".format(*(first,second)
                if hasattr(first,'unit') else (second,first)))
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
    
    
    def failUnlessAlmostRelativeEqual(self, first, second, places=PRECISION, msg=None):
        
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

            
    


class TestWithMPI(TestCase):
    def setUp(self):
        from amuse.legacy.support.core import is_mpd_running
        self.assertTrue(is_mpd_running(), "MPICH2 mpd deamon process not running, cannot run this test as it requires MPI")
            
    def tearDown(self):
        from amuse.legacy.support.core import stop_interfaces
        stop_interfaces()
    
    def new_instance(self, factory, *arguments, **kwarguments):
        try:
            return factory(*arguments)
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
    dir = os.path.abspath(__file__)
    while not os.path.exists(os.path.join(dir,'build.py')):
        dir = os.path.dirname(dir)
    amuse_root_dir = dir
    test_results_dir = os.path.join(amuse_root_dir, 'test_results')
    if os.path.exists(test_results_dir):
        return test_results_dir
    else:
        return './'
