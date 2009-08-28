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

def number_str(number, singular, plural = None):
        if plural == None:
            plural = singular + 's'
        return str(number) + ' ' + (singular if number == 1 else plural)

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



def _run_the_tests(directory, results_queue, id_to_timings):
    #from subprocess import Popen
    
    #self.Popen(["svn", "update"])
    try:
        print "start test run"
        null_device = open('/dev/null')
        os.stdin = null_device
        report = ReportOnATestRun()
        time = TimeATest(id_to_timings)
        plugins = [report , Skip() , time] #Capture() , 
        result = TestProgram(exit = False, argv=['nose', directory], plugins=plugins);
        results_queue.put((report, time.id_to_timings))
    finally:
        MPI.Finalize()

class RunAllTestsOnInterval(object):
    
    def __init__(self, server):
        self.interval_in_seconds = 90 * 60;
        self.must_run = False;
        self.queue = Queue()
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
        while self.must_run:
            p = Process(target=_run_the_tests, args=(os.getcwd(),self.queue,  self.id_to_timings))
            p.start()
            p.join()
            report, self.id_to_timings = self.queue.get() 
            self.server.last_report = report
            self.server.last_timings =  self.id_to_timings
            
            WriteTestReportOnTestingBlog(report, self.id_to_timings).start()
            
            time.sleep(self.interval_in_seconds)

class WriteTestReportOnTestingBlog(object):
    from subprocess import Popen
    
    def __init__(self, report, timings):
        self.report = report
        self.timings = timings
        self.base_directory = os.path.split(__file__)[0]
        self.remote_directory = 'blogs/testing/entries'
        self.local_directory = os.path.join(self.base_directory, "entries")
        if not os.path.exists(self.local_directory):
            os.makedirs(self.local_directory)
        
    def start(self):
        thread = threading.Thread(target=self.run)
        thread.daemon = True;
        thread.start()
        
    def run(self):
        time_struct = time.gmtime(self.report.start_time)
        filename = time.strftime("%Y%m%d_%H_%M.txt", time_struct)
        path = os.path.join(self.local_directory, filename)
        with open(path,"w") as file:
            
            file.write(self.report.title())
            file.write('\n\n')
            
            if self.report.failures > 0:
                file.write('<p>Failed tests:</p>')
                file.write('<ul>')
                for x in self.timings.values():
                    if x.failed:
                        print x
                        filename = x.address[0][len(os.getcwd()):]
                        file.write('<li>')
                        file.write('<a href="/trac/amuse/browser/trunk/'+filename+'#L'+str(x.lineno)+'">')
                        file.write(str(x.id))
                        file.write('</a>')
                        file.write(' - ')
                        delta_t = x.end_time - x.start_time
                        delta_t = round(delta_t, 3)
                        file.write(str(delta_t))
                        file.write(' seconds')
                        file.write('\n')
                        file.write('</li>')
                file.write('</ul>')
                
            if self.report.errors > 0:
                file.write('<p>Errored tests:</p>')
                file.write('<ul>')
                for x in self.timings.values():
                    if x.errored:
                        print x
                        filename = x.address[0][len(os.getcwd()):]
                        file.write('<li>')
                        file.write('<a href="/trac/amuse/browser/trunk/'+filename+'#L'+str(x.lineno)+'">')
                        file.write(str(x.id))
                        file.write('</a>')
                        file.write(' - ')
                        delta_t = x.end_time - x.start_time
                        delta_t = round(delta_t, 3)
                        file.write(str(delta_t))
                        file.write(' seconds')
                        file.write('\n')
                        file.write('</li>')
                file.write('</ul>')
            
            file.write('<p>Tests run:</p>')
            file.write('<ul>')
            for x in self.timings.values():
                filename = x.address[0][len(os.getcwd()):]
                file.write('<li>')
                file.write('<a href="/trac/amuse/browser/trunk/'+filename+'#L'+str(x.lineno)+'">')
                file.write(str(x.id))
                file.write('</a>')
                file.write(' - ')
                delta_t = x.end_time - x.start_time
                delta_t = round(delta_t, 3)
                file.write(str(delta_t))
                file.write(' seconds')
                file.write('\n')
                file.write('</li>')
            file.write('</ul>')
        self.Popen(["scp", path, "doctor:"+self.remote_directory])
        
    
        
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
        self.run_all_tests = RunAllTestsOnInterval(self)
        self.run_all_tests.start()
        self.daemon_threads = True
        
    def start(self):
        self.serve_forever()
        
    def start_run_tests(self):
        p = Process(target=_run_the_tests, args=(os.getcwd(),self.queue))
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
    
