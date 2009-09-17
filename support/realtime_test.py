import time
import urlparse
import threading
import traceback
import json
import nose
import sys
import linecache
import inspect
import os.path
import BaseHTTPServer
import SocketServer
import Queue as queue

from mpi4py import MPI
from nose.plugins.capture import Capture
from nose.plugins.skip import Skip, SkipTest
from nose.core import TestProgram
from multiprocessing import Process, Queue
from optparse import OptionParser
from subprocess import call
from StringIO import StringIO

test_is_running = False

def number_str(number, singular, plural = None):
        if plural == None:
            plural = singular + 's'
        return str(number) + ' ' + (singular if number == 1 else plural)

def find_method_in_class(name_of_the_method, code, class_to_search):
    if name_of_the_method in class_to_search.__dict__:
        member =  class_to_search.__dict__[name_of_the_method]
        if inspect.isfunction(member):
            if member.func_code == code:
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
        #line = linecache.getline(filename, lineno, f.f_globals)
        line = ""
        if '__file__' in f.f_globals:
            for global_name, x in f.f_globals.iteritems():
                if global_name.startswith('_'):
                   continue
                   
                if inspect.isfunction(x):
                    if global_name == name and x.func_code == co:
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

class MonitoredFile(object):
    def __init__(self, path, container):
        self.path = path
        self.container = container
        self.timestamp = self.get_last_modification_time()
      
    def is_file(self):
        return True
  
    def check(self, monitor):
        if not os.path.exists(self.path):
            self.container.remove(self)
            monitor.deleted(self)
            return
        
        measured_timestamp = self.get_last_modification_time()
        if self.timestamp < measured_timestamp:
            self.timestamp = measured_timestamp
            monitor.updated(self)
            return
        
        monitor.unchanged(self)
    
    def get_last_modification_time(self):
        statinfo = os.stat(self.path)
        return statinfo.st_mtime

    
        
class MonitoredDirectory(object):
    def __init__(self, path, container = None):
        self.path = path
        self.elements = []
        self.container = container
        self.path_to_element = {}
        self.setup_from_filesystem()
         
    def is_file(self):
        return False
        
    def setup_from_filesystem(self):
        names = os.listdir(self.path)
        for name in names:
            path = os.path.join(self.path, name)
            element = self.new_element(path)
            self.elements.append(element)
            self.path_to_element[path] = element
    
    def new_element(self, path):
        if os.path.isdir(path):
            return MonitoredDirectory(path, self)
        else:
            return MonitoredFile(path, self)
                
    def remove(self, element):
        self.elements.remove(element)
        del self.path_to_element[element.path]
        
    def check(self, monitor):
        if not os.path.exists(self.path):
            if not self.container is None:
                self.container.remove(self)
            monitor.deleted(self)
            return
            
        for x in self.elements:
            x.check(monitor)
            
        names = os.listdir(self.path)
        for name in names:
            path = os.path.join(self.path, name)
            if not path in self.path_to_element:
                element = self.new_element(path)
                monitor.created(element)
                self.elements.append(element)
                self.path_to_element[path] = element
            
class MonitorDirectories(object):
    def __init__(self, paths):
        self.elements = map(lambda x : MonitoredDirectory(x), paths)
        self.changed = False
        
    def check(self):
        self.changed = False
        for x in self.elements:
            x.check(self)
    
    def deleted(self, monitored_element):
        if monitored_element.path.endswith('.pyc'):
            return
        self.changed = True
        
    def created(self, monitored_element):
        if monitored_element.path.endswith('.pyc'):
            return
            
        self.changed = True
        
    def unchanged(self, monitored_element):
        pass
    
    def updated(self, monitored_element):
        if monitored_element.path.endswith('.pyc'):
            return
        print monitored_element.path
        self.changed = True



class TestCaseReport(object):
    def __init__(self, test):
        self.id = test.id()
        self.address = test.address()
        
        self.start_time = 0.0
        self.end_time = 0.0
        self.total_time = 0.0
        self.number_of_runs = 0.0
        
        self.number_of_suite_runs = 0.0
        
        if hasattr(test.test, "_testMethodName"):
            method = getattr(test.test, getattr(test.test, "_testMethodName"))
            self.lineno = method.func_code.co_firstlineno
        else:
            self.lineno = test.test.descriptor.compat_co_firstlineno
            
        self.failed = False
        self.errored = False
        self.skipped = False
        self.found = False
        
    def start(self):
        self.start_time = time.time()
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

    def end_with_error(self, error_tuple):
        self.end_time = time.time()
        error_type, error_value, error_traceback = error_tuple        
        self.error_string = traceback.format_exception_only(error_type, error_value)
        self.traceback = list(reversed(extract_tb(error_traceback)))
        
        self.number_of_runs = 0
        self.total_time = 0.0
        
        self.failed = False
        self.errored = True
        self.skipped = False
        
    def end_with_failure(self, error_tuple):
        self.end_time = time.time()
        error_type, error_value, error_traceback = error_tuple        
        self.error_string = traceback.format_exception_only(error_type, error_value)
        self.traceback = list(reversed(traceback.extract_tb(error_traceback)))
        
        self.number_of_runs = 0
        self.total_time = 0.0
        
        self.failed = True
        self.errored = False
        self.skipped = False
        
    def end_with_skip(self):
        self.end_time = time.time()
        
        self.skipped = True
        self.failed = False
        self.errored = False
        
    def to_dict(self):
        return self.__dict__.copy()
        
        
class MakeAReportOfATestRun(object):
    score = 2000
    
    def __init__(self, previous_report = None):
        self.errors = 0
        self.failures = 0
        self.tests = 0
        self.start_time = 0
        self.end_time = 0
        self.skipped = 0
        self.problem_text = ""
        
        if previous_report is None:
            self.address_to_report = {}
        else:
            self.address_to_report = previous_report.address_to_report
           
            
        self.name = 'report on a test'
        self.enabled = True

    def addSuccess(self,test):
        report = self.get_report(test)
        report.end()
        
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
    
    def addFailure(self, test, error_tuple):
        self.failures += 1
        report = self.get_report(test)
        report.end_with_failure(error_tuple)
           
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
        for x in self.address_to_report.values():
            x.errored = False
            x.failed = False
            x.skipped = False
            x.found = False
        
    def finalize(self, x):    
        self.end_time = time.time()
        for key, report in list(self.address_to_report.iteritems()):
            if not report.found:
                del self.address_to_report[key]
        pass
        

    def startTest(self, test):
        if self.is_test_able_to_run(test):
           return
        else:
           raise SkipTest

    def is_test_able_to_run(self, test):
        report = self.get_report(test)
        time_taken = report.mean_time() 
        if time_taken < 0.1:
            return True
        if time_taken < 1.0:
            return (report.number_of_suite_runs % 5) == 0
        return (report.number_of_suite_runs % 10) == 0
    
    def to_dict(self):
        result = {}
        for x in [ 'errors', 'failures', 'tests' , 'start_time', 'end_time', 'skipped', 'problem_text']:
            result[x] = getattr(self, x)
        
        testcases = list(self.address_to_report.values())
        testcases.sort(key=lambda x: os.path.basename(x.address[0]))
        result['testcases'] = map(lambda x: x.to_dict(),testcases )
        
        return result  



class Select(object):
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

def _perform_the_testrun(directory, results_queue, previous_report = None):
    try:
        print "start test run"
        null_device = open('/dev/null')
        os.stdin = null_device
        report = MakeAReportOfATestRun(previous_report)
        plugins = [report , Skip(), Capture()]  
        result = TestProgram(exit = False, argv=['nose', directory], plugins=plugins);
        results_queue.put(report)
    except:
        results_queue.put('exception happened')
    finally:
        MPI.Finalize()


def _perform_one_test(directory, results_queue, address):
    try:
        print "start test run"
        null_device = open('/dev/null')
        os.stdin = null_device
        select = Select(address)
        plugins = [select, Skip()]  
        result = TestProgram(exit = False, argv=['nose', directory], plugins=plugins);
        print "end test run"
        result = select.buffer
        result = result[:min(1000, len(result) - 1)]
        results_queue.put(result)
    except:
        results_queue.put('exception happened')
    finally:
        MPI.Finalize()

def _run_test_with_address(address):
    server.run_all_tests.paused = False
    global test_is_running
    if test_is_running:
        return 'test is already running'
    test_is_running = True
    try:
        result_queue = Queue()
        process = Process(
            target=_perform_one_test, 
            args=(
                os.getcwd(), 
                result_queue,
                address
            ))
        process.start()
        
        print "star joined"
        process.join()
        print "process joined"
        result = result_queue.get()
        
        print "test finished:", result
        return result 
    finally:
        test_is_running = False
    

class RunAllTestsWhenAChangeHappens(object):
    
    def __init__(self, server):
        self.must_run = False
        self.server = server
        self.paused = False
        
    def start(self):
        self.must_run = True;
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True;
        self.thread.start()
    
    def stop(self):
        self.must_run = False;
        
    def run(self):
        global test_is_running
        paths = [os.path.join(os.getcwd(), 'src'), os.path.join(os.getcwd(), 'test')]
        monitor = MonitorDirectories(paths)
        monitor.check()
        monitor.changed = True
        while self.must_run:
            if monitor.changed:
                if test_is_running or self.paused:
                    monitor.check()
                    time.sleep(0.5)
                    continue
                    
                test_is_running = True
                try:
                    result_queue = Queue()
                    process = Process(
                        target=_perform_the_testrun, 
                        args=(
                            os.getcwd(), 
                            result_queue, 
                            self.server.last_report
                        ))
                    process.start()
                    process.join()
                    report = result_queue.get() 
                    print report
                    self.server.last_report = report
                    del result_queue
                    self.server.tests_finished.set()
                    
                    print "end test run"
                finally:
                    test_is_running = False
                    
                monitor.check()
            else:
                time.sleep(0.5)
                monitor.check()
                
def open_file(path, lineno = 1):
    call(['geany', path, '+'+str(lineno)])
    
    
class HandleRequest(BaseHTTPServer.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    def do_GET(self):
        parsed_path = urlparse.urlparse(self.path)
        
        if parsed_path.path == '/start':

            thread = threading.Thread(target=self.server.start_run_tests)
            thread.daemon = True;
            thread.start()
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/stop':
            self.stop_server()
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/pause':
            self.pause()
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/get_last_report':
            if self.server.last_report is None:
                string = "null"
            else:
                string = json.dumps(self.server.last_report.to_dict())
            content_type = 'text/javascript'
        elif parsed_path.path == '/open_file':
            parameters = urlparse.parse_qs(parsed_path.query)
            path = parameters['path'][0]
            lineno = int(parameters['lineno'][0])
            open_file(path, lineno)
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/run_test':
            parameters = urlparse.parse_qs(parsed_path.query)
            a0 = parameters['a0'][0]
            a1 = parameters['a1'][0]
            a2 = parameters['a2'][0]
            address = (a0, a1, a2)
            string = json.dumps(_run_test_with_address(address))
            content_type = 'text/javascript'
        elif parsed_path.path == '/events':
            self.do_long_poll()
            return
        else:
            string, content_type = self.index_file()
        
        self.send_response(200)
        self.send_header("Content-type", content_type)
        self.send_header("Content-Length", str(len(string)))
        #self.send_header("Connection", "keep-alive")
        
        self.end_headers()
        self.wfile.write(string)
    
    def do_long_poll(self):
        print self.request_version
        self.send_response(200)
        self.send_header("Content-Type", "text/javascript")   
        #self.send_header("connection", "keep-alive")
        self.send_header("Transfer-Encoding", "chunked")
        self.send_header("Cache-Control", "no-cache, no-store")
        self.send_header("Pragma", "no-cache")
        self.end_headers()
        while True:
            self.server.tests_finished.wait(10.0)
            if self.server.tests_finished.is_set():
                self.send_chunk('true')
                self.server.tests_finished.clear()
            else:
                self.send_chunk('false')
        self.wfile.write('0\r\n\r\n')
        self.wfile.flush()
                    
    def send_chunk(self, string):
        hex_length = hex(len(string))[2:]
        #print hex(len(string)), hex_length
        self.wfile.write('%s \r\n' % hex_length)
        self.wfile.flush()

        self.wfile.write(string)
        self.wfile.write('\r\n')
        self.wfile.flush()
       
    def stop_server(self):
        thread = threading.Thread(target=self.server.stop)
        thread.daemon = True;
        thread.start()
        
    def pause(self):
        server.run_all_tests.paused = True
        
    def index_file(self):
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "realtime_test.html")
        with open(filename, "r") as file:
            contents = file.read()
            return contents, 'text/html'
        

class ContinuosTestWebServer(SocketServer.ThreadingMixIn, BaseHTTPServer.HTTPServer):
    
    def __init__(self, port):
        BaseHTTPServer.HTTPServer.__init__(self, ('',port), HandleRequest)
        self.last_report = None
        self.run_all_tests = RunAllTestsWhenAChangeHappens(self)
        self.run_all_tests.start()
        self.daemon_threads = True
        self.tests_finished = threading.Event()
        
    def start(self):
        self.serve_forever()
        
    def start_run_tests(self):
        server.run_all_tests.paused = False
        global test_is_running
        if test_is_running:
            return
        
        test_is_running = True
        try:
            result_queue = Queue()
            process = Process(
                target=_perform_the_testrun, 
                args=(os.getcwd(), result_queue))
            process.start()
            process.join()
            report = result_queue.get() 
            self.last_report = report
            del result_queue
            self.tests_finished.set()
        finally:
            test_is_running = False
    def stop(self):
        self.run_all_tests.stop()
        self.shutdown()
        
if __name__ == '__main__':
    parser = OptionParser() 
    
    parser.add_option("-p", "--port", 
      dest="serverport",
      help="start serving on PORT", 
      metavar="PORT", 
      default=9070,
      type="int")
      
    (options, args) = parser.parse_args()
    
    print "starting server on port: ", options.serverport
    
    server = ContinuosTestWebServer(options.serverport)
    server.start()
    
