from __future__ import print_function

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


from . import background_test
from . import project
            
background_test.RunTests.instance = background_test.RunTests()


def number_str(number, singular, plural = None):
    if plural == None:
        plural = singular + 's'
    return str(number) + ' ' + (singular if number == 1 else plural)


def get_svn_credentials():
    if 'RVPW' in os.environ:
        svnpassword = os.environ['RVPW']
    else:
        svnpassword = ''
    if 'RVUSR' in os.environ:
        svnusername = os.environ['RVUSR']
    else:
        svnusername = 'reviewboard'
    return svnpassword, svnusername
        
def _run_the_tests(directory, do_update = False):
    svnpassword, svnusername = get_svn_credentials()
    arguments = ["svn", "--no-auth-cache"]
    if svnpassword:
        arguments.extend(('--password', svnpassword))
    if svnusername:
        arguments.extend(('--username', svnusername))
    arguments.append("update")
    call(arguments)
    call(["make", "clean"])
    call(["make", "all"])
    
    print("start test run")
    report = background_test.RunTests.instance.run_tests(None)
    return report
    
class MakeSVNStatusReport(object):
    from xml.dom import minidom 
    
    def start(self):
        
        svnpassword, svnusername = get_svn_credentials()
        arguments = ["svn", "--no-auth-cache"]
        if svnpassword:
            arguments.extend(('--password', svnpassword))
        if svnusername:
            arguments.extend(('--username', svnusername))
        arguments.extend(('info', '--xml'))
        
        process = Popen(arguments, stdout = PIPE, stderr = PIPE)
        stdoutstring, stderrstring = process.communicate()
        doc = self.minidom.parseString(stdoutstring)
        commit_node = list(doc.getElementsByTagName('commit'))[0]
        revision = commit_node.getAttribute('revision')
        self.result = 'SVN revision: {0}'.format(revision)
        print(self.result)
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
        print(title)
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

        testcases.sort(key=lambda x: '' if x.address[0] is None else os.path.basename(x.address[0]))
        
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
            
        if self.report.skipped > 0:
            self.parts.append('<p>Skipped tests:</p>')
            self.parts.append('<ul>')
            for x in testcases:
                if x.skipped:
                    self.write_item_on_test(x)
            self.parts.append('</ul>')
        
        self.parts.append('<p>Tests run:</p>')
        self.parts.append('<ul>')
        for x in testcases:
            self.write_item_on_test(x)
        self.parts.append('</ul>')
        
        with open(path,"w") as file:
            file.write(''.join(self.parts))
            
        call(["scp", "-i", "forreportkey.id", path, "castle.strw.leidenuniv.nl:"+self.remote_directory])
        
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
        if filename is None:
            self.parts.append(str(id))
            return
            
        if id is None:
            id = filename + ':' + str(lineno)
            
        filename = filename[len(os.getcwd()):]
        if filename.endswith('.pyc'):
            filename = filename[:-3] + 'py'
        
        self.parts.append('<a href="/trac/browser/trunk'+filename+'#L'+str(lineno)+'">')
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
      
    parser.add_option("-u", "--update", 
      dest="do_update",
      help="update the code before running the tests", 
      default=True,
      action="store_true")
    parser.add_option("", "--no-update", 
      dest="do_update",
      help="do not update the code before running the tests", 
      action="store_false")
      
    (options, args) = parser.parse_args()
    
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(120 * 60) #building and testing must be done in 2 hours
    os.setpgrp()
    try:
        report = _run_the_tests(options.directory) 
        WriteTestReportOnTestingBlog(report).start()
    finally:
        #os.killpg(0, signal.SIGKILL)
        pass
