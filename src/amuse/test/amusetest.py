import unittest

import numpy
import os
import sys
import inspect

from amuse.support import exceptions
from amuse.support import literature
from amuse.support import options
from amuse.units.si import no_unit
from amuse.units.quantities import Quantity
from amuse.units.quantities import to_quantity
from amuse.units.quantities import is_quantity

class SkipTest(exceptions.AmuseException):
    pass
    
class TestCase(unittest.TestCase):
    PRECISION = int(round(numpy.log10(2.0/(numpy.finfo(numpy.double).eps))))-1
    
    def setUp(self):
        literature.TrackLiteratureReferences.suppress_output()
        
    def _check_comparable(self, first, second):
        if is_quantity(first) is not is_quantity(second):
            # One exception: quantity with none_unit CAN be compared with non-quantity:
            if not to_quantity(first).unit == to_quantity(second).unit:
                raise TypeError("Cannot compare quantity: {0} with non-quantity: {1}.".format(*(first,second)
                    if isinstance(first, Quantity) else (second,first)))
    
    def _convert_to_numeric(self, first, second, in_units):
        if in_units:
            return (first.value_in(in_units), second.value_in(in_units))
        elif is_quantity(first) or is_quantity(second):
            return (to_quantity(first).value_in(to_quantity(second).unit), 
                to_quantity(second).value_in(to_quantity(second).unit))
        else:
            return (first, second)
    
    def _convert_to_vectors(self, first, second):
        x = first if is_quantity(first) else numpy.array(first)
        y = second if is_quantity(second) else numpy.array(second)
        try:
            # Using numpy broadcasting to convert the arguments to arrays with equal length:
            return (x+(y*0)).flatten(), (y+(x*0)).flatten()
        except:
            #raise
            raise TypeError("Arguments do not have compatible shapes for broadcasting")
        
    
    def _raise_exceptions_if_any(self, failures, first, second, err_fmt_string, msg, *args):
        if len(failures) == 1:
            if failures[0]:
                raise self.failureException(msg or err_fmt_string.format(first, second, *args))
        elif any(failures):
            first, second = self._convert_to_vectors(first, second)
            err_list =  [("@{0}, ".format(i)+err_fmt_string.format(first[i], second[i], *args))
                            for (i,b) in enumerate(failures) if b]
            err = '\n'.join(err_list)
            raise self.failureException(msg or err)
    
    
    def failUnlessEqual(self, first, second, msg = None, in_units = None):
        self._check_comparable(first, second)
        first_num, second_num = self._convert_to_numeric(first, second, in_units)
        failures = numpy.array(first_num != second_num).flatten()
        self._raise_exceptions_if_any(failures, first, second, '{0} != {1}', msg)
        
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
        first_num, second_num = self._convert_to_numeric(first, second, in_units)
        failures = numpy.array(numpy.round(second_num - first_num, places) != 0).flatten()
        self._raise_exceptions_if_any(failures, first, second, '{0} != {1} within {2} places', msg, places)
        
    assertAlmostEqual =  failUnlessAlmostEqual                                                    
    assertAlmostEquals =  failUnlessAlmostEqual
    
    
    def failUnlessAlmostRelativeEqual(self, first, second, places=None, msg=None):
        self._check_comparable(first, second)
        first_num, second_num = self._convert_to_numeric(first, second, None)
        
        if places is None:
            places = self.PRECISION
            
        if places <= 0:
            places = self.PRECISION + places
        
        maxRelativeError = 0.1 ** places
        
        diff = numpy.abs(second_num - first_num)
        
        is_one_zero = second_num * first_num == 0.0
        is_failure_for_one_zero =  diff >= (maxRelativeError)
        
        first_for_div = first_num + (is_one_zero * 1.0)
        second_for_div = second_num + (is_one_zero * 1.0)
        
        relative_error = numpy.maximum(numpy.abs(diff /(first_for_div)), numpy.abs(diff /(second_for_div)))        
        is_failure_for_both_nonzero = relative_error >= maxRelativeError
        failures = numpy.where(is_one_zero, is_failure_for_one_zero, is_failure_for_both_nonzero)
        
        failures = failures.flatten()
        
        self._raise_exceptions_if_any(failures, first, second, "{0!r} != {1!r} within {2!r} places", msg, places)
        
    assertAlmostRelativeEqual = failUnlessAlmostRelativeEqual
    assertAlmostRelativeEquals = failUnlessAlmostRelativeEqual
    
    def assertIsOfOrder(self, first, second, msg=None):
        ratio = first*1.0/second
        if isinstance(ratio, Quantity):
            raise self.failureException,(msg or "Units of {0!r} and {1!r} do not match.".format(first, second))
        failures = numpy.array(numpy.round(numpy.log10(ratio)) != 0).flatten()
        self._raise_exceptions_if_any(failures, first, second, '{0} is not of order {1}', msg)
    
    def assertRaises(self, exception, callable, *list_args, **keyword_args):
        if 'expected_message' in keyword_args:
            exception_message = keyword_args.pop('expected_message')
            try:
                callable(*list_args, **keyword_args)
            except exception as ex:
                self.assertEqual(str(ex), exception_message)
            else:
                raise self.failureException("Exception '{0}' was not raised.".format(exception))
        else:
            unittest.TestCase.assertRaises(self, exception, callable, *list_args, **keyword_args)
    
    def run(self, result=None):
        if hasattr(unittest, 'SkipTest'):
            return unittest.TestCase.run(self, result)
            
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
                result.addError(self, sys.exc_info())
                return

            ok = False
            try:
                testMethod()
                ok = True
            except self.failureException:
                result.addFailure(self, sys.exc_info())
            except SkipTest, ex:
                ok = True
                pass
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, sys.exc_info())

            try:
                self.tearDown()
            except KeyboardInterrupt:
                raise
            except:
                result.addError(self, sys.exc_info())
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
            if hasattr(unittest, 'SkipTest'):
                raise unittest.SkipTest(reason)
            else:
                raise SkipTest(reason)
            
    def get_path_to_results(self):
        return get_path_to_results()
    
    def get_amuse_root_dir(self):
        return get_amuse_root_dir()

            
    


class TestWithMPI(TestCase):
    def setUp(self):
        TestCase.setUp(self)
        from amuse.rfi.core import is_mpd_running
        self.assertTrue(is_mpd_running(), "MPICH2 mpd deamon process not running, cannot run this test as it requires MPI")
            
    def tearDown(self):
        from amuse.rfi.core import stop_interfaces
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

    
    def can_compile_modules(self):
        return TestDefaults().can_run_tests_to_compile_modules
        
    def check_can_compile_modules(self):
        if not self.can_compile_modules():
            self.skip('will not run tests that compile codes')
            
    
class TestDefaults(options.OptionalAttributes):
        
    @options.option(sections=['test'])
    def path_to_results(self):
        name_of_testresults_directory = self.name_of_testresults_directory
        if os.path.exists(os.path.abspath(name_of_testresults_directory)):
            if os.access(os.path.abspath(name_of_testresults_directory),  os.W_OK):
                return os.path.abspath(name_of_testresults_directory)
            else:
                result = os.getcwd()
                if not os.access(os.getcwd(),  os.W_OK):
                    raise Exception("the current directory must be writable for amuse tests")
                return os.getcwd()
        
        amuse_root_dir = self.amuse_root_dir
        test_results_dir = os.path.join(amuse_root_dir, self.name_of_testresults_directory)
        if os.path.exists(test_results_dir):
            try:
                f = open('test.txt','w')
                f.close()                    
                return test_results_dir
            except IOError as ex:
                if not os.access(os.getcwd(),  os.W_OK):
                    raise Exception("the current directory must be writable for amuse tests")
                return os.getcwd()
        else:
            if not os.access(os.getcwd(),  os.W_OK):
                raise Exception("the current directory must be writable for amuse tests")
            return os.getcwd()

    @options.option(sections=['test'])
    def name_of_testresults_directory(self):
        return 'test_results'
    
    @options.option(sections=['data'])
    def amuse_root_dir(self):
        if 'AMUSE_DIR' in os.environ:
            return os.environ['AMUSE_DIR']    
        previous = None
        result = os.path.abspath(__file__)
        while not os.path.exists(os.path.join(result,'build.py')):
            result = os.path.dirname(result)
            if result == previous:
                return os.path.dirname(os.path.dirname(__file__))
            previous = result
        return result
    
    @options.option(type='boolean',sections=['test'])
    def can_run_tests_to_compile_modules(self):
        return True



def get_path_to_results():
    return TestDefaults().path_to_results

def get_amuse_root_dir():
    return TestDefaults().amuse_root_dir
    
