from __future__ import print_function

import time
import traceback
import sys
import linecache
import inspect
import os.path
try:  # Python 2
    import Queue as queue
    from StringIO import StringIO
    func_code_attr = 'func_code'
except ImportError:  # Python 3
    import queue
    from io import StringIO
    func_code_attr = '__code__'
import subprocess
import threading
import tempfile
import shutil

from mpi4py import MPI

from nose.core import TestProgram
from nose.plugins.capture import Capture
from nose.plugins.skip import Skip, SkipTest
from nose.plugins.doctests import Doctest

from multiprocessing import Process, Queue

from . import project


def is_mpd_running():
        
    name_of_the_vendor, version = MPI.get_vendor()
    if name_of_the_vendor == 'MPICH2':
        must_check_mpd = True
        if 'AMUSE_MPD_CHECK' in os.environ:
            must_check_mpd = os.environ['AMUSE_MPD_CHECK'] == '1'
        if 'PMI_PORT' in os.environ:
            must_check_mpd = False
        if 'HYDRA_CONTROL_FD' in os.environ:
            must_check_mpd = False
        
        if not must_check_mpd:
            return True
        try:
            process = subprocess.Popen(['mpdtrace'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (output_string, error_string) = process.communicate()
            return not (process.returncode == 255)
        except OSError as ex:
            return True
    else:
        return True
        
def ensure_mpd_is_running():
    if not is_mpd_running():
        name_of_the_vendor, version = MPI.get_vendor()
        if name_of_the_vendor == 'MPICH2':
            try:
                process = subprocess.Popen(['nohup','mpd'], close_fds=True)
            except OSError as ex:
                pass



def number_str(number, singular, plural = None):
    if plural == None:
        plural = singular + 's'
    return str(number) + ' ' + (singular if number == 1 else plural)

def find_method_in_class(name_of_the_method, code, class_to_search):
    if name_of_the_method in class_to_search.__dict__:
        member =  class_to_search.__dict__[name_of_the_method]
        if inspect.isfunction(member):
            if getattr(member, func_code_attr) == code:
                return member
        if inspect.ismethoddescriptor(member):
            pass
    return None

def extract_tb(tb, limit = None):
    list = []
    n = 0
    while tb is not None and (limit is None or n < limit):
        f = tb.tb_frame
        lineno = tb.tb_lineno
        co = f.f_code
        filename = co.co_filename
        name = co.co_name
        linecache.checkcache(filename)
        line = ""
        if '__file__' in f.f_globals:
            for global_name, x in f.f_globals.items():
                if global_name.startswith('_'):
                    continue

                if inspect.isfunction(x):
                    if global_name == name and get(x, func_code_attr) == co:
                        args, varargs, varkw, defaults = inspect.getargspec(x)
                        name += inspect.formatargspec(args, varargs, varkw, defaults)
                elif inspect.isclass(x):
                    method = find_method_in_class(name,co,x)
                    if not method is None:
                        args, varargs, varkw, defaults = inspect.getargspec(method)
                        name += inspect.formatargspec(args, varargs, varkw, defaults)
                        name = x.__name__ + '.' + name
                            
        if line:
            line = line.strip()
        else: 
            line = None
        list.append((filename, lineno, name, line))
        tb = tb.tb_next
        n = n+1
    return list


class TestCaseReport(object):
    type = 'unit-report'
    
    def __init__(self, test):
        self.id = test.id()
        self.address = test.address()
        
        self.start_time = 0.0
        self.end_time = 0.0
        self.total_time = 0.0
        self.number_of_runs = 0.0
        
        self.number_of_suite_runs = 0.0
        self.is_paused = False
        
        
        if hasattr(test.test, "_dt_test"):
            self.lineno = test.test._dt_test.lineno
        elif hasattr(test.test, "_testMethodName"):
            method = getattr(test.test, getattr(test.test, "_testMethodName"))
            self.lineno = method.getattr(func_code_attr).co_firstlineno
        else:
            self.lineno = test.test.descriptor.compat_co_firstlineno
            
        self.failed = False
        self.errored = False
        self.skipped = False
        self.found = False
    
    def __str__(self):
        return str(self.id)
    
    def start(self):
        self.start_time = time.time()
        self.end_time = 0.0
        self.found = True
    
    def end(self):
        self.end_time = time.time()
        self.number_of_runs += 1.0
        self.total_time += (self.end_time - self.start_time)
        self.failed = False
        self.errored = False
        self.skipped = False
        
    def mean_time(self):
        if self.number_of_runs == 0:
            return 0.0
        return self.total_time / self.number_of_runs

    def add_syntax_error_file_to_traceback(self, error_type, error_value):
        if not issubclass(error_type, SyntaxError):
            return
            
        try:
            msg, (filename, lineno, offset, badline) = error_value.args
        except Exception:
            pass
        
        self.traceback.insert(0, (filename, lineno, "", ""))
        
    def end_with_error(self, error_tuple):
        self.end_time = time.time()
        error_type, error_value, error_traceback = error_tuple        
        self.error_string = traceback.format_exception_only(error_type, error_value)
        self.traceback = list(reversed(extract_tb(error_traceback)))
        
        self.add_syntax_error_file_to_traceback(error_type, error_value)
        
        self.reset_timing()
        
        self.failed = False
        self.errored = True
        self.skipped = False
        
    def end_with_failure(self, error_tuple):
        self.end_time = time.time()
        error_type, error_value, error_traceback = error_tuple        
        self.error_string = traceback.format_exception_only(error_type, error_value)
        self.traceback = list(reversed(traceback.extract_tb(error_traceback)))
        self.add_syntax_error_file_to_traceback(error_type, error_value)
        
        self.reset_timing()
        
        self.failed = True
        self.errored = False
        self.skipped = False
        
    def end_with_skip(self):
        self.end_time = time.time()
        
        self.skipped = True
        self.failed = False
        self.errored = False
        
    def reset_timing(self):
        self.number_of_runs = 0
        self.total_time = 0.0
            
    def to_dict(self):
        result = self.__dict__.copy()
        result['mean_time'] = self.mean_time()
        return result
        
        
class MakeAReportOfATestRun(object):
    score = 2000
    
    def __init__(self, previous_report = None, reports_queue = None):
        self.errors = 0
        self.failures = 0
        self.tests = 0
        self.start_time = 0
        self.end_time = 0
        self.skipped = 0
        self.report_id = -1
        self.last_test_run = None
        self.crashed = False
        
        if reports_queue is None:
            self._queue_report = self.ignore_report
        else:
            self._queue_report = self.store_report
            self._reports_queue = reports_queue
        
        if previous_report is None:
            self.address_to_report = {}
        else:
            self.address_to_report = previous_report.address_to_report
           
            
        self.name = 'report on a test'
        self.enabled = True

    def __getstate__(self):
        result = {}
        for key, value in self.__dict__.iteritems():
            if not key.startswith('_'):
                result[key] = value
        return result
        
        
    def ignore_report(self, report):
        pass
        
    def store_report(self, report):
        if self.total_number_of_tests() % 1 == 0 or report.failed or report.errored:
            self._reports_queue.put((report.type,report.to_dict(), self.to_information_dict() ))
        
    
    def prepareTest(self, suite):
        pass
        #all_tests = list(suite)
        #number_of_tests = 0
        #self._reports_queue.put(('start-report',{"number_of_tests": number_of_tests}))
        
    def addSuccess(self,test):
        report = self.get_report(test)
        report.end()
        self._queue_report(report)
        
    def addError(self,test,error_tuple):
        report = self.get_report(test)
        error_class, ignore_1, ignore_2 = error_tuple
        if issubclass(error_class, SkipTest):
            self.skipped += 1
            self.tests -= 1
            report.end_with_skip()
        else: 
            report.end_with_error(error_tuple)
            self.errors += 1
        
        self._queue_report(report)
    
    def addFailure(self, test, error_tuple):
        self.failures += 1
        report = self.get_report(test)
        report.end_with_failure(error_tuple)
        
        self._queue_report(report)
           
    def get_report(self, test):
        address = test.address()
        if not address in self.address_to_report:
            self.address_to_report[address] = TestCaseReport(test)
        return self.address_to_report[address] 
        
    def options(self, parser, env):
        pass
        
    def configure(self, parser, env):
        pass
        
        
    def beforeTest(self,test):
        self.tests += 1
        x = self.get_report(test)
        x.start()
        x.number_of_suite_runs += 1
    
    def begin(self):
        self.start_time = time.time()
        self.end_time = 0
        for x in self.address_to_report.values():
            x.errored = False
            x.failed = False
            x.skipped = False
            x.found = False
        
    def finalize(self, x):    
        self.end_time = time.time()
        for key, report in list(self.address_to_report.items()):
            if not report.found:
                del self.address_to_report[key]        

    def startTest(self, test):
        if self.is_test_able_to_run(test):
            report = TestCaseStartReport(test)
            self._queue_report(report)
            return
        else:
            raise SkipTest

    def is_test_able_to_run(self, test):
        report = self.get_report(test)
        if report.is_paused:
            return False
            
        time_taken = report.mean_time()
        if time_taken < 0.1:
            return True
        if time_taken < 0.25:
            return (report.number_of_suite_runs % 5) == 0
        if time_taken < 0.5:
            return (report.number_of_suite_runs % 10) == 0
        if time_taken < 1.0:
            return (report.number_of_suite_runs % 15) == 0
        if time_taken < 3.0:
            return (report.number_of_suite_runs % 30) == 0
            
        return (report.number_of_suite_runs % 40) == 0
    
    def has_errors(self):
        return self.errors > 0
    
    def has_failures(self):
        return self.failures > 0
        
    def has_skipped_tests(self):
        return self.skipped > 0
        
    def total_number_of_tests(self):
        return self.tests + self.skipped  + self.failures + self.errors
        
    def title_string(self):
        title = '';
        if self.has_failures() or self.has_errors():
            title += "FAIL "
        else:
            title += "OK "

        
        if self.has_errors():
            title += 'E' + str(self.errors) + ' '
        
        if self.has_failures():
            title += 'F' + str(self.failures) + ' '
        
        
        title += '+' + str(self.tests) + ' '
        
        if self.has_skipped_tests():
            title += 'S' + str(self.skipped) + ' '
        
        title += 'T' + str(self.total_number_of_tests()) + ' '
        
        if not self.is_running():
            delta_time = self.end_time - self.start_time
            title += ("%10.3f" % delta_time) + 's';
        else:
            cur_time = time.time()
            delta_time =  cur_time - self.start_time
            title += "(running {0:10.3f})".format(delta_time)
            
        title += ' ';
        title += time.strftime("%H:%M:%S", time.gmtime(self.start_time))
        return title;

    def is_running(self):
        return self.end_time <= 1.0
        
    def to_information_dict(self):
        result = {}
        if not self.is_running():
            result['report_id'] = self.report_id
        else:
            result['report_id'] = -1
        result['title'] = self.title_string()
        result['success'] = not self.has_errors() and not self.has_failures()
        return result
    
    def to_dict(self):
        result = self.to_information_dict()
        
        for x in [ 
            'errors', 'failures', 'tests', 
            'start_time', 'end_time', 'skipped',
            ]:
            result[x] = getattr(self, x)
        
        testcases = list(self.address_to_report.values())
        for x in testcases:
            if x.address[0] is None:
                x.address = list(x.address)
                x.address[0] = ''
                
        testcases.sort(key=lambda x: os.path.basename(x.address[0]))
        result['testcases'] = [x.to_dict() for x in testcases]

        return result



class SelectOneTestAndStoreOutput(object):
    name = 'select one test'
    enabled = True
    score = 10
    
    def __init__(self, address):
        self.address = address
        self.capture_stdout = None
        self.stdout = None
        
    def options(self, parser, env):
        pass
        
    def configure(self, parser, env):
        pass
    
    def afterTest(self, test):
        if not self.stdout is None:
            sys.stdout = self.stdout
            self.stdout = None
            
    def startTest(self, test):
        if test.address() == self.address:
            self.stdout = sys.stdout
            self.capture_stdout = StringIO()
            sys.stdout = self.capture_stdout 
            return
        else:
            raise SkipTest
        
    @property
    def buffer(self):
        if self.capture_stdout:
            return self.capture_stdout.getvalue()
        else:
            return "none"

class RunTests(object):
    DIRECTORIES = None
    WORKING_DIRECTORY = None
    
    def __init__(self):
        self.test_is_running = False
        self.report_queue = queue.Queue()
        self.report_info = {}
        self.life_reports = []
        self.get_life_reports_semaphore = threading.Semaphore()
        
    def _perform_the_testrun(self, directories, results_queue, previous_report = None):
        try:
            ensure_mpd_is_running()
            
            null_device = open('/dev/null')
            os.stdin = null_device
            report = MakeAReportOfATestRun(previous_report, results_queue)
            doctest = Doctest()
            doctest.enabled = True
            plugins = [doctest, report, Skip(), Capture()] 
            argv = ['nose', '-v']
            old_working_directory = os.getcwd()
            if not self.WORKING_DIRECTORY is None:
                argv.extend(['-w', self.WORKING_DIRECTORY])
                os.chdir(self.WORKING_DIRECTORY)
                
            argv.extend( directories)
            argv.extend(['--with-doctest', '--doctest-extension=txt'])
            
            result = TestProgram(exit = False, argv=argv, plugins=plugins);
            
            os.chdir(old_working_directory)
            results_queue.put(('test-report', report,) )
        except :
            results_queue.put(('test-error', 'Exception happened: ' + str(sys.exc_info()[0]) + " - " + str(sys.exc_info()[1]), ))
        finally:
            results_queue.put(None)
            MPI.Finalize()


    def _perform_one_test(self, directories, results_queue, address):
        try:
            print("start test run")
            null_device = open('/dev/null')
            os.stdin = null_device
            select = SelectOneTestAndStoreOutput(address)
            plugins = [select, Skip()]  
            argv = ['nose']
            argv.extend( directories)
            result = TestProgram(
                exit = False, 
                argv=argv, 
                plugins=plugins)
            success = result.success
            if success:
                result = 'Success'
            else:
                result = 'Failure'
            
            if select.buffer:
                result += ' - '
                if len(select.buffer) > 1000:
                    result += select.buffer[:min(1000, len(result) - 1)] 
                    result += ' ' + str(result - 1000) + ' more ...'
                else:
                    result += select.buffer
            
            results_queue.put(('test-output',result))
        except:
            results_queue.put(('test-error', 'Exception happened: ' + str(sys.exc_info()[0]) + " - " + str(sys.exc_info()[1]), ))
        finally:
            results_queue.put(None)
            MPI.Finalize()
            print("calling finalize done")

    def run_test_with_address(self, address):
        typestring, result = self.run_in_another_process(self._perform_one_test, address)
        return result
         
    def run_tests(self, previous_report):
        typestring, result = self.run_in_another_process(self._perform_the_testrun, previous_report)
        if typestring == 'test-error':
            print(result)
            raise Exception(result)
        else:
            return result  
            
    def can_start_a_test(self):
        return not self.test_is_running
        
    def run_in_another_process(self, method, argument):
        if self.test_is_running:
            raise Exception("Test is already running")
        
        self.life_reports = []
        self.clear_reports_queue()
        self.test_is_running = True
        self.last_test_started = None
        try:
            result_queue = Queue()
            cwd = os.getcwd()
            if self.DIRECTORIES is None:
                paths = [os.path.abspath(x) for x in project.DIRECTORIES]
            else:
                paths = [os.path.abspath(x) for x in self.DIRECTORIES]
        
            process = Process(
                target=method, 
                args=(
                    paths, 
                    result_queue, 
                    argument
                ))
            
            previous_output = sys.stdout
            previous_error = sys.stderr
            
            temp_stdout = tempfile.NamedTemporaryFile('w')
            temp_stderr = temp_stdout

            print(temp_stdout.name)

            sys.stdout = temp_stdout
            sys.stderr = temp_stderr
            
            process.start()
            
            sys.stdout = previous_output
            sys.stderr = previous_error
            last_message = None
            while True:
                message = result_queue.get(True, 60)
                if message is None:
                    break;
                if message[0] == 'unit-report':
                    print(message[1])
                    self.report_info = message[2]
                    if(message[1]['failed'] or message[1]['errored']):
                        self.report_queue.put(message[1])
                   
                if message[0] == 'start-report':
                    self.last_test_started = message[1]
                    print(self.last_test_started)
                    continue
                    
                last_message = message
            result = last_message
            self.last_test_started = None
        except queue.Empty:
            print("No message recieved from process for 60 seconds")
            report = MakeAReportOfATestRun()
            report.last_test_run = self.last_test_started
            report.start_time = report.end_time = time.time()
            report.crashed = True
            result = 'failed', report
            temp_stdout.flush()
            
            shutil.copyfile(temp_stdout.name, "test_output.txt")
            
            process.terminate()
        finally:
            temp_stdout.close()
            self.test_is_running = False
            self.clear_reports_queue()
            
        process.join(2)
        del result_queue
        return result
     
    def clear_reports_queue(self):
        
        must_clear = not self.report_queue.empty()
        while must_clear:
            try:
                self.report_queue.get_nowait()
                must_clear = not self.report_queue.empty()
            except queue.Empty:
                must_clear = False
                
    def get_reports(self):
        self.get_life_reports_semaphore.acquire()
        try:
            must_clear = not self.report_queue.empty()
            while must_clear:
                try:
                    self.life_reports.append(self.report_queue.get_nowait())
                    must_clear = not self.report_queue.empty()
                except queue.Empty:
                    must_clear = False
            return self.life_reports
        finally:
            self.get_life_reports_semaphore.release()
        
            




class TestCaseStartReport(object):
    type = 'start-report'
    
    def __init__(self, test):
        self.id = test.id()
        self.address = test.address()
        self.start_time = time.time()
            
    def to_dict(self):
        result = self.__dict__.copy()
        return result

