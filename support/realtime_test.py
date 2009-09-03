import sys
import sqlite3
import urlparse
import urllib
import cgi
import os.path
import BaseHTTPServer, SocketServer

from mpi4py import MPI
import json
import nose
import nose.plugins
from nose.plugins.capture import Capture
import time
from nose.core import TestProgram
from nose.plugins.skip import Skip, SkipTest
import threading

from multiprocessing import Process, Queue
import Queue as queue

from optparse import OptionParser
from subprocess import call


def number_str(number, singular, plural = None):
        if plural == None:
            plural = singular + 's'
        return str(number) + ' ' + (singular if number == 1 else plural)


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
        self.changed = True
        
    def created(self, monitored_element):
        self.changed = True
        
    def unchanged(self, monitored_element):
        pass
    
    def updated(self, monitored_element):
        self.changed = True


class ReportOnATestRun(object):
    score = 2000
    
    def __init__(self, original = None):
        if original is None:
            self.errors = 0
            self.failures = 0
            self.tests = 0
            self.start_time = 0
            self.end_time = 0
            self.skipped = 0
            self.problem_text = ""
        else:
            self.errors = original.errors
            self.failures = original.failures
            self.tests = original.tests
            self.start_time = original.start_time
            self.end_time = original.end_time
            self.skipped = original.skipped
            self.problem_text = original.problem_text
            
        self.name = 'report on a test'
        self.enabled = True

    def addError(self,test,err):
        error_class, u1, u2 = err
        if issubclass(error_class, SkipTest):
            self.skipped += 1
            self.tests -= 1
        else: 
            self.errors += 1
            self.problem_text += '\nerror:'
            #self.problem_text += str(error_class)
            self.problem_text += str(u1)
            self.problem_text += '\n   '
            self.problem_text += str(test)
            #self.problem_text += test.shortDescription()
            pass
            
    def to_dict(self):
        result = {}
        for x in ['errors', 'failures', 'tests' , 'start_time', 'end_time', 'skipped', 'problem_text']:
            result[x] = getattr(self, x)
        return result
        
    def options(self, parser, env):
        pass
        
    def configure(self, parser, env):
        pass
        
    def addFailure(self, test, err):
        error_class, u1, u2 = err
        self.failures += 1
        self.problem_text += '\nassertion failed:'
        self.problem_text +=  str(u1)
        self.problem_text += '\n   '
        self.problem_text += str(test)
        
    def beforeTest(self,test):
        self.tests += 1
        
    def begin(self):
        self.start_time = time.time()

        
    def finalize(self, x):    
        self.end_time = time.time()
        pass
        
    def __str__(self):
        w = []
        if self.failures > 0:
            w.append(number_str(self.failures,'failure'))
            w.append(' ')
        if self.errors > 0:
            w.append(number_str(self.errors,'error'))
            w.append(' ')
        if self.errors == 0 and self.failures == 0:
            w.append('all test passed')
        w.append(' ')
        w.append('- ')
        delta_t = self.end_time - self.start_time
        delta_t = round(delta_t, 3)
        w.append(str(delta_t))
        w.append(' s')
        w.append('\n')
        w.append(number_str(self.tests,'test'))
        if self.skipped > 0:
            w.append(', ')
            w.append(number_str(self.skipped,'test'))
            w.append(' skipped')
        if self.problem_text:
            w.append('\n\n')
            w.append(problem_text)
        return ''.join(w)
        
    def title(self):
        w = []
        w.append(number_str(self.tests,'test'))
        w.append(' run, ')
        if self.failures > 0:
            w.append(number_str(self.failures,'failure')) 
            w.append(' ')
        if self.errors > 0:
            w.append(number_str(self.errors,'error'))
            w.append(' ')
        if self.errors == 0 and self.failures == 0:
            w.append('all tests passed')
        w.append(' - ')
        delta_t = self.end_time - self.start_time
        delta_t = round(delta_t, 3)
        w.append(str(delta_t))
        w.append(' seconds')
        return ''.join(w)

class TimingsOfOneTest(object):
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
        
    def start(self):
        self.start_time = time.time()
    
    def end(self):
        self.end_time = time.time()
        self.number_of_runs += 1.0
        self.total_time += (self.end_time - self.start_time)
        
    def mean_time(self):
        if self.number_of_runs == 0:
            return 0.0
        return self.total_time / self.number_of_runs

class TimeATest(object):
    
            
    def __init__(self, id_to_timings =  {}):
        self.id_to_timings =  id_to_timings
        self.enabled = True
        
    def beforeTest(self, test):
        timings = self.get_timings(test)
        timings.start()
        timings.number_of_suite_runs += 1
    
    def is_test_able_to_run(self, test):
        timings = self.get_timings(test)
        time_taken = timings.mean_time() 
        if time_taken < 0.1:
            return True
        if time_taken < 1.0:
            return (timings.number_of_suite_runs % 5) == 0
        return (timings.number_of_suite_runs % 10) == 0
    
    def addSuccess(self,test):
        timings = self.get_timings(test)
        timings.end()
        timings.failed = False
        timings.errored = False
        
    def addFailure(self, test, err):
        timings = self.get_timings(test)
        timings.end_time = time.time()
        timings.number_of_runs = 0
        timings.total_time = 0.0
        timings.errored = False
        timings.failed = True
        
    def addError(self, test, err):
        timings = self.get_timings(test)
        timings.end_time = time.time()
        timings.number_of_runs = 0
        timings.total_time = 0.0
        timings.failed = False
        timings.errored = True
        
    def startTest(self, test):
        if True or self.is_test_able_to_run(test):
           return
        else:
           raise SkipTest
           
    def get_timings(self, test):
            id = test.address()
            if not id in self.id_to_timings:
                self.id_to_timings[id] = TimingsOfOneTest(test)
            return self.id_to_timings[id] 



def perform_the_testrun(directory, results_queue, id_to_timings):
    try:
        print "start test run"
        null_device = open('/dev/null')
        os.stdin = null_device
        report = ReportOnATestRun()
        time = TimeATest(id_to_timings)
        plugins = [report , Skip(), Capture() , time]  
        result = TestProgram(exit = False, argv=['nose', directory], plugins=plugins);
        results_queue.put((report, time.id_to_timings))
    finally:
        MPI.Finalize()

class RunAllTestsWhenAChangeHappens(object):
    
    def __init__(self, server):
        self.must_run = False
        self.server = server
        self.id_to_timings = {}
        
    def start(self):
        self.must_run = True;
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True;
        self.thread.start()
    
    def stop(self):
        self.must_run = False;
        
    def run(self):
        paths = [os.path.join(os.getcwd(), 'src'), os.path.join(os.getcwd(), 'test')]
        monitor = MonitorDirectories(paths)
        monitor.check()
        while self.must_run:
            monitor.check()
            if monitor.changed:
                result_queue = Queue()
                p = Process(target=perform_the_testrun, args=(os.getcwd(), result_queue, self.id_to_timings))
                p.start()
                p.join()
                report, self.id_to_timings = result_queue.get() 
                self.server.last_report = report
                self.server.last_timings =  self.id_to_timings
                del result_queue
            else:
                time.sleep(0.5)
                
        
class HandleRequest(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_GET(self):
        parsed_path = urlparse.urlparse(self.path)
        
        if parsed_path.path == '/start':
            self.start_run_tests()
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/stop':
            self.stop_server()
            string = 'null'
            content_type = 'text/javascript'
        elif parsed_path.path == '/get_last_report':
            if self.server.last_report is None:
                string = "null"
            else:
                string = json.dumps(self.server.last_report.to_dict())
            content_type = 'text/javascript'
        else:
            string, content_type = self.index_file()
        
        self.send_response(200)
        self.send_header("Content-type", content_type)
        self.send_header("Content-Length", str(len(string)))
        #self.send_header("Connection", "keep-alive")
        
        self.end_headers()
        self.wfile.write(string)
        
       
    def stop_server(self):
        thread = threading.Thread(target=self.server.stop)
        thread.daemon = True;
        thread.start()
        
    def index_file(self):
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "continuous_test.html")
        with open(filename, "r") as file:
            contents = file.read()
            return contents, 'text/html'
        

class ContinuosTestWebServer(SocketServer.ThreadingMixIn, BaseHTTPServer.HTTPServer):
    
    def __init__(self, port):
        BaseHTTPServer.HTTPServer.__init__(self, ('',port), HandleRequest)
        self.last_report = None
        self.last_timings = None
        self.run_all_tests = RunAllTestsWhenAChangeHappens(self)
        self.run_all_tests.start()
        self.daemon_threads = True
        
    def start(self):
        self.serve_forever()
        
    def start_run_tests(self):
        p = Process(target=perform_the_testrun, args=(os.getcwd(),self.queue))
        p.start()
    
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
    
