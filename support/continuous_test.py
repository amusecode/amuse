import sys
import os
import signal
import time
import datetime
import nose
import nose.plugins
import traceback
import linecache
import inspect

from nose.plugins.capture import Capture
from nose.core import TestProgram
from nose.plugins.skip import Skip, SkipTest
from nose.plugins.doctests import Doctest
from optparse import OptionParser
from subprocess import call, Popen, PIPE


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
    

class TestCaseReport(object):
    def __init__(self, test):
        self.id = test.id()
        self.address = test.address()
        
        self.start_time = 0.0
        self.end_time = 0.0
        self.total_time = 0.0
        self.number_of_runs = 0.0
        
        self.number_of_suite_runs = 0.0
        
        if hasattr(test.test, "_dt_test"):
            self.lineno = test.test._dt_test.lineno
        elif hasattr(test.test, "_testMethodName"):
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
        for x in [ 
            'errors', 'failures', 'tests' , 
            'start_time', 'end_time', 'skipped',
            'problem_text']:
            result[x] = getattr(self, x)
        
        testcases = list(self.address_to_report.values())
        testcases.sort(key=lambda x: os.path.basename(x.address[0]))
        result['testcases'] = map(lambda x: x.to_dict(),testcases )
        
        return result  



def _run_the_tests(directory):
    
    print "updating the code"
    call(["svn", "update"])
    call(["make", "clean"])
    call(["make", "all"])
    
    print "start test run"
    null_device = open('/dev/null')
    os.stdin = null_device
    report = MakeAReportOfATestRun()
    plugins = [report , Skip() ,Capture(), Doctest()] 
    result = TestProgram(exit = False, argv=['nose', directory, '--with-doctest'], plugins=plugins);
    return report
    
class MakeSVNStatusReport(object):
    from xml.dom import minidom 
    
    def start(self):
        process = Popen(['svn','info' , '--xml'], stdout = PIPE, stderr = PIPE)
        stdoutstring, stderrstring = process.communicate()
        doc = self.minidom.parseString(stdoutstring)
        commit_node = list(doc.getElementsByTagName('commit'))[0]
        revision = commit_node.getAttribute('revision')
        self.result = 'SVN revision: {0}'.format(revision)
        print self.result
        return self
        
class MakePlatformReport(object):
    import platform
    import mpi4py.MPI
    
    def start(self):
        w = []
        w.append('<table>')
        items = ['uname', 'system', 'node', 'release', 'version', 'machine', 'processor', 'architecture']
        if self.platform.system() == 'Linux':
            items.append('linux_distribution')
        if self.platform.system() == 'Darwin':
            items.append('mac_ver')
                    
        for x in items:
            self.addRow(w, x, str(getattr(self.platform, x)()))
            
        name_of_the_vendor, version = self.mpi4py.MPI.get_vendor()
        self.addRow(w, 'MPI vendor', str(name_of_the_vendor))
        self.addRow(w, 'MPI version', str(version))
        
        w.append('</table>')
        self.result = ''.join(w)
        return self
    
    def addRow(self, w,  string1, string2):
        w.append('<tr>') 
        w.append('<td>')
        w.append(string1)
        w.append('</td>')
        w.append('<td>')
        w.append(string2)
        w.append('</td>')
        w.append('</tr>') 
        w.append('\n')


    

        
class WriteTestReportOnTestingBlog(object):
    
    def __init__(self, report):
        self.report = report
        self.base_directory = os.path.split(__file__)[0]
        self.remote_directory = 'blogs/testing/entries'
        self.local_directory = os.path.join(self.base_directory, "entries")
        if not os.path.exists(self.local_directory):
            os.makedirs(self.local_directory)
    
    def xtitlestring_for_the_report(self):
        delta_time = self.report.end_time - self.report.start_time;
        title = '';
        title += str(datetime.datetime.fromtimestamp(self.report.start_time));
        title += ' ';
        if report.errors > 0:
            title += 'E' + str(report.errors) + ' ';
        
        if report.failures > 0:
            title += 'F' +  str(report.failures) + ' ';
        
        title += '+' + str(report.tests) + ' ';
        if report.skipped > 0:
            title += 'S' +  str(report.skipped) + ' ';
        
        title += str(delta_time) + 's';
        print title
        return title
        
    def titlestring_for_the_report(self):
        w = []                                                                                         
        w.append(number_str(self.report.tests,'test'))                                                        
        w.append(' run, ')                                                                             
        if self.report.failures > 0:                                                                          
            w.append(number_str(self.report.failures,'failure'))                                              
            w.append(' ')                                                                              
        if self.report.errors > 0:                                                                            
            w.append(number_str(self.report.errors,'error'))                                                  
            w.append(' ')                                                                              
        if self.report.errors == 0 and self.report.failures == 0:                                                    
            w.append('all tests passed')                                                               
        w.append(' - ')                                                                                
        delta_t = self.report.end_time - self.report.start_time                                                      
        delta_t = round(delta_t, 3)                                                                    
        w.append(str(delta_t))                                                                         
        w.append(' seconds')                                                                           
        return ''.join(w)  
        
    def start(self):
        time_struct = time.gmtime(self.report.start_time)
        filename = time.strftime("%Y%m%d_%H_%M.txt", time_struct)
        path = os.path.join(self.local_directory, filename)
        self.parts = []
        
        self.parts.append(self.titlestring_for_the_report())
        self.parts.append('\n\n')
        
        self.parts.append('<p>')
        self.parts.append(MakeSVNStatusReport().start().result)
        self.parts.append('</p>')
        self.parts.append('\n\n')
        self.parts.append(MakePlatformReport().start().result)
        self.parts.append('\n\n')
        
        testcases = list(self.report.address_to_report.values())
        testcases.sort(key=lambda x: os.path.basename(x.address[0]))
        
        if self.report.failures > 0:
            self.parts.append('<p>Failed tests:</p>')
            self.parts.append('<ul>')
            for x in testcases:
                if x.failed:
                    self.write_item_on_test(x)
            self.parts.append('</ul>')
            
        if self.report.errors > 0:
            self.parts.append('<p>Errored tests:</p>')
            self.parts.append('<ul>')
            for x in testcases:
                if x.errored:
                    self.write_item_on_test(x)
            self.parts.append('</ul>')
        
        self.parts.append('<p>Tests run:</p>')
        self.parts.append('<ul>')
        for x in testcases:
            self.write_item_on_test(x)
        self.parts.append('</ul>')
        
        with open(path,"w") as file:
            file.write(''.join(self.parts))
            
        call(["scp", path, "castle.strw.leidenuniv.nl:"+self.remote_directory])
        
    def write_item_on_test(self, x):
        self.parts.append('<li>')
        self.write_link_to_file_in_svn(x.address[0], x.lineno, x.id)
        self.parts.append(' - ')
        delta_t = x.end_time - x.start_time
        delta_t = round(delta_t, 3)
        self.parts.append(str(delta_t))
        self.parts.append(' seconds')
        self.parts.append('\n')
        if x.failed or x.errored:
            self.write_error_and_traceback(x)
        self.parts.append('</li>')
    
    def write_error_and_traceback(self, x):
        self.parts.append('<ul>')
        self.parts.append('<li>')
        for y in x.error_string:
            self.parts.append(y)
        self.parts.append('</li>')
        self.parts.append('<li>') 
        self.parts.append('<ul>')
        for filename, lineno, name, line in x.traceback:
            self.parts.append('<li>') 
            self.write_link_to_file_in_svn(filename, lineno)
            self.parts.append(' - ')
            self.parts.append(name)
            self.parts.append(' - ')
            self.parts.append(str(line))
            self.parts.append('</li>')
        self.parts.append('</ul>')
        self.parts.append('</li>')
        self.parts.append('</ul>')
        
    def write_link_to_file_in_svn(self, filename, lineno, id = None):
        if id is None:
            id = filename + ':' + str(lineno)
            
        filename = filename[len(os.getcwd()):]
        if filename.endswith('.pyc'):
            filename = filename[:-3] + 'py'
        
        self.parts.append('<a href="/trac/amuse/browser/trunk'+filename+'#L'+str(lineno)+'">')
        self.parts.append(str(id))
        self.parts.append('</a>')
        

def handler(signum, frame):
    if signum == signal.SIGALRM:
        sys.exit("tests took too much time")
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--dir", 
      dest="directory",
      help="run tests in DIRECTORY", 
      metavar="DIRECTORY", 
      default=os.getcwd(),
      type="string")
      
    (options, args) = parser.parse_args()
    
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(1200) #building and testing must be done in 20 minutes

    report = _run_the_tests(options.directory) 
    WriteTestReportOnTestingBlog(report).start()
    
    
