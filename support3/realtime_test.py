
import sys
if sys.hexversion < 0x03000000:
    from unittest.case import SkipTest
    raise SkipTest('incorrect python version')
    
import time
import urllib.parse
import threading
import json
import os.path
import webbrowser


from optparse import OptionParser

from . import webserver
from . import monitor
from . import background_test
from . import project

background_test.RunTests.instance = background_test.RunTests()

class RunAllTestsWhenAChangeHappens(object):
   
    
    def __init__(self, server):
        self.must_run = False
        self.server = server
        
    def start(self):
        self.must_run = True;
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True;
        self.thread.start()
    
    def stop(self):
        self.must_run = False;
        
    def run(self):
        cwd = os.getcwd()
        paths = [os.path.join(cwd, x) for x in project.DIRECTORIES] 

        monitor_directories = monitor.MonitorDirectories(paths)
        monitor_directories.check()
        monitor_directories.changed = True
        while self.must_run:
            if monitor_directories.changed:
                if not background_test.RunTests.instance.can_start_a_test():
                    monitor_directories.check()
                    time.sleep(0.5)
                    continue
                
                if not self.server.last_report is None:
                    for element in monitor_directories.updated_elements:
                        if not element.is_file():
                            continue
                        for x in list(self.server.last_report.address_to_report.values()):
                            path, module, testcase =  x.address
                            if path == element.path:
                                x.reset_timing()
                                print("will rerun: ", module, testcase)
                        
                
                print("Changed files:")
                number_of_changed_files = 0
                for element in monitor_directories.updated_elements:
                    if not element.is_file():
                        continue
                    print(element.path)
                    number_of_changed_files += 1
                    
                if number_of_changed_files > 0 or self.server.last_report is None:
                    last_report = self.server.last_report
                    
                    self.server.set_last_report(None)
                    
                    report = background_test.RunTests.instance.run_tests(last_report)
                
                    self.server.set_last_report(report)
             
                monitor_directories.check()
            else:
                time.sleep(0.5)
                monitor_directories.check()
                
                
    
class HandleRequest(webserver.HandleRequest):
    
    def do_start(self):
        self.server.restart_testrunner()
        string = 'null'
        content_type = 'text/javascript'
        return string, content_type
    
    
    def do_pause(self):
        self.server.stop_testrunner()
        return 'null', 'text/javascript' 
    
    def do_get_last_report(self):
        string = json.dumps(self.server.get_last_report_as_dict())
        content_type = 'text/javascript'
        return string, content_type

    def do_get_last_report_information(self):
        string = json.dumps(self.server.get_last_report_information())
        content_type = 'text/javascript'
        return string, content_type
    
    def do_run_test(self):
        parameters = urllib.parse.parse_qs(self.parsed_path.query)
        a0 = parameters['a0'][0]
        a1 = parameters['a1'][0]
        a2 = parameters['a2'][0]
        address = (a0, a1, a2)
        result = background_test.RunTests.instance.run_test_with_address(address)
        string = json.dumps(result)
        content_type = 'text/javascript'
        self.server.continue_testrunner()
        return string, content_type
        
    def index_file(self):
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "realtime_test.html")
        with open(filename, "r") as file:
            contents = file.read()
            return contents, 'text/html'
            

class ContinuosTestWebServer(webserver.WebServer):
    
    def __init__(self, port):
        webserver.WebServer.__init__(self,  port, HandleRequest)
        self.last_report = None
        self.run_all_tests = RunAllTestsWhenAChangeHappens(self)
        self.run_all_tests.start()
        self.report_id = 0
        
        
    def stop(self):
        self.run_all_tests.stop()
        self.shutdown()
        
    def restart_testrunner(self):
        self.run_all_tests.stop()
        self.last_report = None
        self.run_all_tests = RunAllTestsWhenAChangeHappens(self)
        self.run_all_tests.start()
        
    def stop_testrunner(self):
        self.run_all_tests.stop()        
    
    def continue_testrunner(self):
        if not self.run_all_tests.must_run:
            self.restart_testrunner()
        
    def get_last_report_as_dict(self):
        if self.last_report is None:
            return None
        else:
            result = self.last_report.to_dict()
            return result
            
    def get_last_report_information(self):
        if self.last_report is None:
            result =  self.get_live_report_info()
        else:
            result =  self.last_report.to_information_dict()
        result['reports'] = self.get_live_reports()
        return result
    
    def get_live_reports(self):
        return background_test.RunTests.instance.get_reports()
    
    def get_live_report_info(self):
        return background_test.RunTests.instance.report_info
        
    def set_last_report(self, report):
        self.last_report = report
        if not report is None:
            self.report_id += 1
            self.last_report.report_id = self.report_id
        

def start_browser(serverport):
    time.sleep(2.0)
    webbrowser.open("http://localhost:{0}/".format(serverport))
            
if __name__ == '__main__':
    parser = OptionParser() 
    
    
    parser.add_option("-p", "--port", 
      dest="serverport",
      help="start serving on PORT", 
      metavar="PORT", 
      default=9070,
      type="int")
      
    parser.add_option("-e", "--editor", 
      dest="editor",
      help="preferred EDITOR for editing the files", 
      metavar="EDITOR", 
      default="geany",
      type="string")

    parser.add_option("-b", "--browser", 
      dest="startbrowser",
      help="automatically start a browser", 
      metavar="PORT", 
      default="yes",
      type="string")
      
    (options, args) = parser.parse_args()
    
    print("starting server on port: ", options.serverport)
    print("will use editor: ", options.editor)
    webserver.EDITOR = options.editor
    
    if options.startbrowser == "yes":
        thread = threading.Thread(target = start_browser, args = (options.serverport,))
        thread.start()
    
    server = ContinuosTestWebServer(options.serverport)
    server.start()
    
