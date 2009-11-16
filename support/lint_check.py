import subprocess
import xml.dom
import xml.dom.minidom
import re
import os
import os.path
from optparse import OptionParser
from collections import namedtuple
import monitor
import webserver
import json
import urlparse
import threading

ReportMessageLine = namedtuple('ReportMessageLine', 'filename lineno state message_id method_id message_string')
ReportOnAFile = namedtuple('ReportOnAFile', 'filename state_to_number_of_messages')

class PylintReport(object):
    
    def __init__(self, messages):
        self.messages = messages
        
    def filenames(self):
        result = set([])
        for message in self.messages:
            result.add(message.filename)
        return sorted(result)
    
    def filereports(self):
        filename_to_report = {}
        for message in self.messages:
            report = filename_to_report.setdefault(
                message.filename, 
                ReportOnAFile(
                    message.filename, 
                    self.new_state_to_number_of_messages_dict()
                )
            )
            report.state_to_number_of_messages[message.state] += 1
        return sorted(filename_to_report.values(), key=lambda x : x.filename)
        
    def get_messages_of_file(self, filename):
        result = []
        for message in self.messages:
            if message.filename == filename:
                result.append(message)
        return result
        
    def to_dict(self):
        result = {}
        result['filereports'] = [x._asdict() for x in self.filereports()]
        return result
    
    def new_state_to_number_of_messages_dict(self):
        return {'C':0, 'W':0, 'F':0, 'R':0, 'E':0}
            
class InterfaceToPyLint(object):
    NAME_OF_THE_COMMAND = 'pylint' 
    
    def is_available(self):
        process = subprocess.Popen(
            ['which', self.NAME_OF_THE_COMMAND],
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE,
        )
        process.communicate()
        return process.returncode == 0
        
    def run_onfile(self, path, extra_python_paths):
        environment = self.get_environment_with_pythonpath(extra_python_paths)
        
        process = subprocess.Popen(
            [self.NAME_OF_THE_COMMAND, '-f', 'parseable', '-i', 'y', path], 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE,
            env = environment
        )
        stdout_string, stderr_string  = process.communicate()
        return self.parse_report(stdout_string)
        
    def parse_report(self, string):
        messages = []
        lines = string.splitlines()
        regexp = re.compile("(.+?)\:(\d+)\: \[(.)(\d\d\d\d)(?:, (.+?))?\] (.*)")
        number_of_statements_re = re.compile("(\d+) statements analysed.")
        for line in lines:
            match = regexp.match(line)
            if not match is None:
                messages.append(ReportMessageLine(*match.groups()))
            match = number_of_statements_re.match(line)
            if not match is None:
                print match.groups()

                print match.group(1)
                number_of_statements = int(match.group(1))
            
        return PylintReport(messages)
        
    def get_environment_with_pythonpath(self, extra_python_paths):
        environment = os.environ.copy()
        new_pythonpath =  ':'.join(extra_python_paths)
        
        if 'PYTHONPATH' in environment:
            new_pythonpath = environment['PYTHONPATH'] + ':' + new_pythonpath
        
        environment['PYTHONPATH'] = new_pythonpath
        return environment

CSS_STRING = """
.numbers {background-color: #eee;}
.lines {background-color: #ffffe0 ;}
.state {background-color: #eee; width: 1.5em; text-align: center;}
.state a {padding-left: 0.5em; padding-right: 0.5em; text-decoration: none;}
.set-C {background-color: lightblue;}
.unset-C {}
.set-R {background-color: #666;}
.unset-R {}
.set-W {background-color: yellow;}
.unset-W {}
.set-E {background-color: red;}
.unset-E {}
.set-F {background-color: orange;}
.unset-F {}


.state  a:visited, a {
    color: black;
}
.state  a:hover {
    color: #CD5C5C;
}

"""

class MakeHTMLDomFromFile(object):
    
    def __init__(self, path, messages):
        self.path = path
        self.filecontents = self.get_filecontents()
        self.messages = messages
        self.number_of_lines_in_the_file = self.get_number_of_lines_in_the_file()
        
    def get_filecontents(self):
        with open(self.path, "r") as file:
            return file.read()
            
    def get_number_of_lines_in_the_file(self):
        lines = self.filecontents.splitlines()
        return len(lines)
        
    def start(self):
        self.document = self.new_document()
        
        self.head = self.document.createElement("head")
        script = self.document.createElement("script")
        script.setAttribute("type","text/javascript")
        script.setAttribute("src", "http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.js")
        self.head.appendChild(script)
        self.document.documentElement.appendChild(self.head)
        
        self.body = self.document.createElement("body")
        self.document.documentElement.appendChild(self.body)
        
        table = self.new_table(1,7)
        tr = table.firstChild
        
        td1 = tr.firstChild
        td2 = tr.lastChild
        td1.setAttribute("class", "numbers")
        td2.setAttribute("class", "lines")
        for i in range(1,6):
            tr.childNodes[i].setAttribute("class", "state")
        
        state_characters = ('C', 'R', 'W', 'E', 'F')
        for i, state_character in enumerate(state_characters):
            self.add_state_line(state_character, tr.childNodes[i+1])
        
        counts = [str(x+1) for x in range(self.number_of_lines_in_the_file)]
        pre = self.document.createElement("pre")
        pre.appendChild(self.document.createTextNode('\n'.join(counts)))
        td1.appendChild(pre)
        
        pre = self.document.createElement("pre")
        file_data = self.document.createCDATASection(self.filecontents)
        pre.appendChild(file_data)
        self.document.documentElement.setAttribute("xmlns",xml.dom.XHTML_NAMESPACE)
        td2.appendChild(pre)
        self.body.appendChild(table)
        
        self.add_stylesheet()
        
        self.result = self.document.toxml()
        self.document.unlink()
    
    def add_state_line(self, state_character, td):
        line_number_to_messages = {}
        for message in self.messages:
            if message[2] == state_character:
                line_number = int(message[1])
                message_strings = line_number_to_messages.setdefault(line_number, [])
                message_strings.append(message[5])
                line_number_to_messages[line_number] = message_strings
        
        pre = self.document.createElement("pre")
        for line_number in range(1, self.number_of_lines_in_the_file+1):
            span = self.document.createElement("a")
            span.appendChild(self.document.createTextNode(state_character))
            pre.appendChild(span)
            pre.appendChild(self.document.createTextNode('\n'))
            filename = os.path.join(os.getcwd(), self.path)
            span.setAttribute("href","#");
            span.setAttribute("onclick","javascript:$.getJSON('/open_file', {'path':'"+filename+"', 'lineno':"+str(line_number)+"}, function handle_result(data){}); false");
            
            if line_number in line_number_to_messages:
                span.setAttribute("class", "set-"+state_character)
                span.setAttribute("title", str(line_number_to_messages[line_number]))
            else:
                span.setAttribute("class", "unset-"+state_character)
        
        td.appendChild(pre)
                
        
    def new_table(self, number_of_rows, number_of_columns):
        table = self.document.createElement("table")
        for row in range(number_of_rows):
            tr = self.document.createElement("tr")
            table.appendChild(tr)
            for column in range(number_of_columns):
                td = self.document.createElement("td")
                tr.appendChild(td)
        return table
        
    def new_document(self):
        implementation = xml.dom.minidom.getDOMImplementation()
        return implementation.createDocument(xml.dom.XHTML_NAMESPACE, "xhtml", None)

    
    def add_stylesheet(self):
        self.style = self.document.createElement("style")
        self.style.setAttribute("type", "text/css")
        self.head.appendChild(self.style)
        stylesheet_string = CSS_STRING
        stylesheet_cdata = self.document.createCDATASection(stylesheet_string)
        self.style.appendChild(stylesheet_cdata)
        
        

class HandleRequest(webserver.HandleRequest):

    def index_file(self):
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "lint_check.html")
        with open(filename, "r") as file:
            contents = file.read()
            return contents, 'text/html'
    
    def do_start(self):
        self.server.start_lint()
        return 'null', 'text/javascript'
    
    def do_show_file(self):
        parameters = urlparse.parse_qs(self.parsed_path.query)
        filename = parameters['path'][0]
        messages = self.server.last_report.get_messages_of_file(filename)
        x = MakeHTMLDomFromFile(filename, messages)
        x.start()
        return x.result, 'application/xhtml+xml'
            
    def do_get_last_report(self):
        string = json.dumps(self.server.get_last_report_as_dict())
        content_type = 'text/javascript'
        return string, content_type
        
        
class LintWebServer(webserver.WebServer):
    
    def __init__(self, port):
        webserver.WebServer.__init__(self,  port, HandleRequest)
        self.run_lint()
        
    def start_lint(self):
        thread = threading.Thread(target=self.run_lint)
        thread.start()
        
    def run_lint(self):
        print "running lint"
        self.set_last_report(InterfaceToPyLint().run_onfile(
            'src/amuse', 
            [os.path.join(os.getcwd(), 'src/amuse/support')]))
        print "done..."
        
    def get_last_report_as_dict(self):
        return self.last_report.to_dict()
    
    def set_last_report(self, report):
        self.last_report = report
        self.events_queue.put('done')
        
        
        
def Run():
    #monitor_directories = monitor.MonitorDirectories(['src', 'support', 'test'])
    #monitor_directories.check()
    #monitor_directories.walk(lambda x : a.append(x))
    pylint = InterfaceToPyLint()
    if not pylint.is_available():
        print "Error, pylint is not available."
        print "please install pylint first, this can be done with 'easy_install pylint'"
        sys.exit(1)
        
    parser = OptionParser() 
    
    
    parser.add_option("-p", "--port", 
      dest="serverport",
      help="start serving on PORT", 
      metavar="PORT", 
      default=9071,
      type="int")
      
    parser.add_option("-e", "--editor", 
      dest="editor",
      help="preferred EDITOR for editing the files", 
      metavar="EDITOR", 
      default="geany",
      type="string")
      
    (options, args) = parser.parse_args()
    
    print "starting server on port: ", options.serverport
    webserver.EDITOR = options.editor
    
    server = LintWebServer(options.serverport)
    server.start()
    

    
        
if __name__ == '__main__':
    Run()
    #print InterfaceToPyLint().run_onfile('src/amuse/support/data', [os.path.join(os.getcwd(), 'src')]).to_dict()
