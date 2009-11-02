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
from subprocess import call, Popen, PIPE
from StringIO import StringIO

EDITOR = None


                
                
osascript_to_open_xcode = """on run argv
 set linenumber to (item 1 of argv) as integer
 set filename_string to item 2 of argv
 set file_to_open to POSIX file filename_string
 tell application "Xcode"
  activate
  set doc_to_edit to (open file_to_open)
  tell doc_to_edit 
   set its selection to item linenumber of paragraph  of it
  end tell
 end tell
end run"""

def open_file(path, lineno = 1):
    global EDITOR
    
    if sys.platform == 'darwin':        
        program = Popen(
            ['osascript', '-', str(lineno), os.path.join(os.getcwd(), path) ], 
            stdin = PIPE, stdout = PIPE, stderr = PIPE)
        out, err = program.communicate(osascript_to_open_xcode)
    else:
        possible_programs = (
            ['geany', path, '+'+str(lineno)],
            ['kate', '-u', '--line',str(lineno),path],
            ['emacs', '+'+str(lineno), path],
            ['nedit-client','-line', str(lineno), path],
        )
        
        for program in possible_programs:
            if program[0] == EDITOR:
                returncode = call(['which', program[0]])
                if returncode == 0:
                    call(program)
                    return 
        
        for program in possible_programs:
            returncode = call(['which', program[0]])
            if returncode == 0:
                call(program)
                return 
        
        call([EDITOR, path])
    
class HandleRequest(BaseHTTPServer.BaseHTTPRequestHandler):
   
    
    def do_GET(self):
        self.parsed_path = urlparse.urlparse(self.path)
        path = self.parsed_path.path[1:]
        method_name = 'do_' + path
        if hasattr(self, method_name):
           method = getattr(self,method_name)
           string, content_type =  method()
        else:
           string, content_type = self.index_file()
        
        self.send_response(200)
        self.send_header("Content-type", content_type)
        self.send_header("Content-Length", str(len(string)))        
        self.end_headers()
        self.wfile.write(string)
    
    def do_long_poll(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/javascript")  
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
        self.wfile.write('%s \r\n' % hex_length)
        self.wfile.flush()

        self.wfile.write(string)
        self.wfile.write('\r\n')
        self.wfile.flush()
       
    def index_file(self):
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "realtime_test.html")
        with open(filename, "r") as file:
            contents = file.read()
            return contents, 'text/html'
            
    def log_message(self, format, *args):
        pass
        #sys.stderr.write("%s - - [%s] %s\n" %
        #                 (self.address_string(),
        #                  self.log_date_time_string(),
        #                  format%args))
        
    def do_stop(self):
        thread = threading.Thread(target=self.server.stop)
        thread.daemon = True;
        thread.start()
        return 'null', 'text/javascript'
        
    def do_events(self):
        new_events = self.server.get_all_events_since_previous_query()
        string = json.dumps(new_events)
        content_type = 'text/javascript'
        return string, content_type
        
    def do_open_file(self):
        parameters = urlparse.parse_qs(self.parsed_path.query)
        path = parameters['path'][0]
        lineno = int(parameters['lineno'][0])
        open_file(path, lineno)
        string = 'null'
        content_type = 'text/javascript'
        return string, content_type
        

class WebServer(SocketServer.ThreadingMixIn, BaseHTTPServer.HTTPServer):
    
    def __init__(self, port, request_handler):
        BaseHTTPServer.HTTPServer.__init__(self, ('', port), request_handler)
        self.daemon_threads = True
        self.events_queue = queue.Queue()
        
    def start(self):
        self.serve_forever()
        
    def stop(self):
        self.shutdown()
        
    def get_all_events_since_previous_query(self):
        try:
            events = []
            while True:
                events.append(self.events_queue.get(False))
        except queue.Empty:
            pass
            
        return events
            
            

    
