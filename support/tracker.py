import subprocess
import socket, os

import time
import urlparse
import threading
import json
import os.path
import os
import sys
import re

import webbrowser

from xml.dom import minidom

from optparse import OptionParser

import webserver
import background_test
import project
import pickle
import textwrap

import Queue

background_test.RunTests.instance = background_test.RunTests()

class late(object):
    def __init__(self, initializer):
        
        self.initializer = initializer
        self.__doc__ = self.initializer.__doc__
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        value = self.initializer(instance)
        setattr(instance,self.initializer.__name__, value)
        return value
        
class SendAnEmail(object):
    
    def __init__(self, **keyword_arguments):
        if len(keyword_arguments) > 0:
            for key, value in keyword_arguments.iteritems():
                setattr(self, key, value)
                self.start()

    def start(self):
        
        call = ['mail','-s',self.subject, '-r', self.sender_email_address]
        call.extend(self.recipients)
        
        print call
        
        process = subprocess.Popen(
            call,
            stdout=subprocess.PIPE,
            stdin =subprocess.PIPE,
            stderr=subprocess.PIPE
        )
            
        stdout, stderr = process.communicate(self.mail_contents)
        
        if process.returncode != 0:
            raise Exception("Could not send e-mail, error output was {0}".format(stderr))
        
    @late
    def sender_email_address(self):
        if "EMAIL" in os.environ:
            return os.environ["EMAIL"]
        else:
            return 'noreply@{0}'.format(socket.getfqdn())
    
    @late
    def mail_contents(self):
        return 'Mail send by the SendAnEmail class\nContents not provided\n'
    
    @late
    def subject(self):
        return 'Automatic mail subject'
        


def get_first_element_with_tag(parent, name):
    for node in parent.childNodes:
        if node.nodeType == minidom.Node.ELEMENT_NODE and \
            (name == "*" or node.tagName == name):
            return node
    return None
    
header = """\
Dear {name},

"""

errored_start = """\
This is to inform you that your commit had errors.
"""
success_start = """\
This is to inform you that your commit had no errors, well done!
"""

commit_info = """\
AUTHOR   : {author}
REVISION : {revision}
DATE     : {date}
MESSAGE  : {msg}
"""

error_info = """\
NUMBER OF ERRORS        : {number_of_errors:>5d}
"""

skip_info = """\
NUMBER OF SKIPPED TESTS : {number_of_skips:>5d}
"""

tests_info = """\
NUMBER OF TESTS         : {number_of_tests:>5d}
TIME TAKEN(seconds)     : {number_of_seconds:>6.1f}
"""
footer = """\
Regards,

The AMUSE automated testing system
"""

more_tests = """\
The number of tests has increased with : {delta_tests}. Good job!
"""
less_tests = """\
Watch out, the number of tests has decreased! There are now {0} less tests. 
"""
same_tests = """\
The number of tests has not increased, is your code tested?
"""

errored_email_subject = """Found {number_of_errors} error(s) in revision {revision}"""
success_email_subject = """For revision revision {revision}, all {number_of_tests} tests were successful!"""


def run_command(arguments):
    print "running :" + ' '.join(arguments)
    process = subprocess.Popen(arguments, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdoutstr, stderrstr = process.communicate()
    print stderrstr
    return stdoutstr

class RequestACodeReview(object):
    
    def __init__(self, **keyword_arguments):
        if len(keyword_arguments) > 0:
            for key, value in keyword_arguments.iteritems():
                setattr(self, key, value)
                self.start()

    def start(self):
        if not self.log_string:
            pass
        
        if subprocess.call(['which', 'post-review']) != 0:
            return
            
        revision = self.revision
        
        repository_url =  '--repository-url=' + self.svn_repository_url
        password = '--password=' + self.reviewboard_password
        username = '--username=' + self.reviewboard_username
        description = "--description=(In [%s]) %s" % (revision, self.log_string)
        submitas = '--submit-as=' + self.author
        revision = '--revision-range=%s:%s' % (int(revision)-1, revision)
        server = '--server=' + self.reviewboard_url
        
        arguments = [
            'post-review',
            repository_url,
            password,
            username,
            submitas,
            revision,
            server,
            self.review_identifier_argument,
            self.publish_argument,
            self.testing_done_argument
        ]
        
        if len(self.review_identifier_argument) == 0:
            arguments += [self.summary_argument, description]
        print run_command(arguments)
        
    @late
    def svn_repository_url(self):
        return 'http://www.amusecode.org/svn/'
    
    @late
    def reviewboard_url(self):
        return 'http://www.amusecode.org/review/'
        
    @late
    def reviewboard_username(self):
        return 'svn'
        
    @late
    def reviewboard_password(self):
        return os.environ['RVPW']

    @late
    def author(self):
        return 'unknown'
    
    @late
    def log_string(self):
        return ''
    
    @late
    def review_identifier_argument(self):
        m = re.search(r'update(?: )?review:([0-9]+)', self.log_string, re.M | re.I)
        if m:
            return '--review-request-id=' + m.group(1)
        else:
            return ''
            
    @late
    def publish_argument(self):
        if re.search(r'draft(?: )?review', self.log_string, re.M | re.I):
            return ''
        else:
            return '-p'
    
    @late
    def summary_argument(self):
        return '--summary=' + self.log_string[:250].splitlines().pop(0).split('. ').pop(0)
    
    @late
    def testing_done_argument(self):
        if self.test_report:
            return '--testing-done=' + self.test_report
        else:
            return ''



            
class RunAllTestsOnASvnCommit(object):
    DEFAULT = None
    
    def __init__(self):
        self.queue =  Queue.Queue()
        self.must_run = False

    @classmethod
    def default(cls):
        if cls.DEFAULT is None:
            cls.DEFAULT = cls()
        return cls.DEFAULT
    
    @late
    def working_directory(self):
        path = os.path.dirname(os.path.dirname(__file__))
        path = os.path.abspath(os.path.join(path, 'working-copy'))
        return path
    
    @late
    def directories(self):
        return [os.path.join(self.working_directory, x) for x in project.DIRECTORIES]
    
    @late
    def mapping_from_author_to_email(self):
        path = os.path.dirname(__file__)
        path = os.path.join(path, "authors.map")
        with open(path, 'r') as f:
            return pickle.load(f)
    
    def cleanup_compiled_python_files(self):
        for path in self.find_pyc_files(self.working_directory):
            os.remove(path)

    def find_pyc_files(self, rootname):
        for dirname, subdirectories, files in os.walk(rootname):
            for filename in files:
                if filename.endswith('.pyc'):
                    yield os.path.join(dirname, filename)


    
    def update_from_svn(self, revision):
        subprocess.call(['svn', '--force', 'update', '-r', revision], cwd = self.working_directory)
        
    def build_code(self):
        subprocess.call(['make','clean'], cwd = self.working_directory)
        subprocess.call(['make'], cwd = self.working_directory)
        
    def get_author_date_and_msg_for(self, revision):
        process = subprocess.Popen(['svn','log', '-r', revision, '--xml'], cwd = self.working_directory, stdout = subprocess.PIPE)
        result, ignore = process.communicate()
        print result
        if not process.returncode == 0:
            raise Exception("could not retrieve log for revision {0}" + revision)
            
        doc = minidom.parseString(result)
        results = []
        entry = list(doc.getElementsByTagName('logentry'))[0]
        author  =  get_first_element_with_tag(entry, "author").firstChild.data
        date_string = get_first_element_with_tag(entry, "date").firstChild.data
        msg_string = get_first_element_with_tag(entry, "msg").firstChild.data
            
        return author, date_string, msg_string
        
    def run_all_tests(self):
        background_test.RunTests.DIRECTORIES = self.directories
        background_test.RunTests.WORKING_DIRECTORY = self.working_directory
        return background_test.RunTests.instance.run_tests(None)
        
    def send_report_as_email_to(self, report, mail_content_string, recipient):
        uc = SendAnEmail()
        uc.mail_contents = mail_content_string
        
        if report["number_of_errors"] > 0:
            uc.subject = errored_email_subject.format(**report)
        else:
            uc.subject = success_email_subject.format(**report)
        
        if not recipient is None:
            uc.recipients = [recipient]
            uc.start()
        
        uc.recipients = [self.admin_email_address]
        uc.start()
        
    def new_mail_content_string(self, report):
        contents = []
        contents.append(header.format(**report))
        
        if report["number_of_errors"] > 0:
            contents.append(errored_start.format(**report))
        else:
            contents.append(success_start.format(**report))
            
        
        if report['delta_tests'] > 0:
            contents.append(more_tests.format(**report))
        elif report['delta_tests'] < 0:
            contents.append(less_tests.format(-1 * report['delta_tests']))
        else:
            contents.append(same_tests.format(**report))
            
        contents.append(commit_info.format(**report))
        
        if report["number_of_errors"] > 0:
            contents.append(error_info.format(**report))
            
        if report["number_of_skips"] > 0:
            contents.append(skip_info.format(**report))
            
            
        contents.append(tests_info.format(**report))
            
        
        if report['errors']:
            for location_line, error_string in report['errors']:
                contents.append(location_line)
                contents.append(textwrap.fill('\n'.join(error_string), 80, initial_indent = "  "))
            contents.append('')
                
        contents.append(footer.format(**report))
        
        return '\n'.join(contents)
        
    def check_svn_commit(self,revision):
        
        self.cleanup_compiled_python_files()
        
        
        author, date, msg = self.get_author_date_and_msg_for(revision)
        if author in self.mapping_from_author_to_email:
            name, email = self.mapping_from_author_to_email[author]
        else:
            name, email = 'Admin', None
        
        self.update_from_svn(revision)
        self.build_code()
        
        test_report = self.run_all_tests()
        
        report = {}
        report['author'] = author
        report['date'] = date
        report['msg'] = msg
        report['name'] = name
        report['number_of_errors'] = test_report.errors + test_report.failures
        report['number_of_skips'] = test_report.skipped
        report['number_of_tests'] = test_report.tests
        report['number_of_seconds'] = test_report.end_time - test_report.start_time
        report['revision']  = revision
        
        previous_revision = int(revision) - 1
        while not previous_revision in self.mapping_from_revision_to_report and previous_revision > 0:
            previous_revision -= 1
            
        if previous_revision in self.mapping_from_revision_to_report:
            previous_report = self.mapping_from_revision_to_report[previous_revision]
            previous_number_of_tests = previous_report['number_of_tests']
            delta_in_number_of_tests = test_report.tests - previous_number_of_tests
        else:
            delta_in_number_of_tests = 0
            
        report['delta_tests'] = delta_in_number_of_tests
        
        
        testcases = list(test_report.address_to_report.values())
        testcases.sort(key=lambda x: os.path.basename("" if x.address[0] is None else x.address[0]))
       
      
        errors = []
        for x in testcases:
            if x.errored or x.failed:
                if not x.address[0] is None:
                    filename = os.path.basename(x.address[0])
                else:
                    filename = 'unknown-file'
                
                location_line = "{0}:{1} {2}".format(filename, x.lineno, x.address[2])
                error_string = x.error_string
                
                errors.append((location_line,error_string,))
                
        report['errors'] = errors
        
        self.mapping_from_revision_to_report[int(revision)] = report
        self.dump_revision_reports()
        
        content = self.new_mail_content_string(report)
        self.send_report_as_email_to(report, content,  email)
        
        requestACodeReview = RequestACodeReview()
        requestACodeReview.log_string = msg
        requestACodeReview.author = author
        requestACodeReview.revision = revision
        requestACodeReview.test_report = content
        requestACodeReview.start()
        
        
    @late
    def admin_email_address(self):
        if "EMAIL" in os.environ:
            return os.environ["EMAIL"]
        else:
            return None
            
    def start(self):
        self.thread = threading.Thread(target = self.runloop)
        self.thread.daemon = True
        self.must_run = True
        self.thread.start()
    
          
    def runloop(self): 
        while self.must_run:
            revision = self.queue.get()
            if revision is None:
                self.must_run = False
            else:
                self.check_svn_commit(revision)
        
        
    def queue_check(self, revision):
        self.queue.put(revision)
        
    def stop(self):
        self.queue.put(None)
    
    def number_of_tested_revisions(self):
        return len(self.mapping_from_revision_to_report)
        
    def dump_revision_reports(self):
        path = os.path.dirname(__file__)
        path = os.path.join(path, "reports.map")
        with open(path, 'w') as f:
            pickle.dump(self.mapping_from_revision_to_report, f)
    
    @late
    def mapping_from_revision_to_report(self):
        path = os.path.dirname(__file__)
        path = os.path.join(path, "reports.map")
        if os.path.exists(path):
            try:
                with open(path, 'r') as f:
                    return pickle.load(f)
            except IOError:
                return {}
        else:
            return {}
        
        
class HandleRequest(webserver.HandleRequest):
   
    def do_check_svn_commit(self):
        parameters = urlparse.parse_qs(self.parsed_path.query)
        revision = parameters['rev'][0]
        
        
        self.server.tracker.queue_check(revision)
        
        string = json.dumps("test for revision {0} queued".format(revision))
        content_type = 'text/javascript'
        return string, content_type
        
    
    def do_ping(self):
        string = json.dumps(True)
        content_type = 'text/javascript'
        return string, content_type
        
    def index_file(self):
        if True:
            return ("nothing here", "text/html")
            
        base = os.path.split(__file__)[0]
        filename = os.path.join(base, "tracker.html")
        with open(filename, "r") as file:
            contents = file.read()
        return contents, 'text/html'
            

class ContinuosTestWebServer(webserver.WebServer):
    
    def __init__(self, port):
        webserver.WebServer.__init__(self,  port, HandleRequest)
        self.tracker = RunAllTestsOnASvnCommit()
        
        
    def stop(self):
        self.tracker.stop()
        self.shutdown()
        
            
if __name__ == '__main__':
    parser = OptionParser() 
    
    parser.add_option("-p", "--port", 
      dest="serverport",
      help="start serving on PORT", 
      metavar="PORT", 
      default=9075,
      type="int")
    
    parser.add_option("-a", "--admin", 
      dest="admin_email",
      help="e-mail address of the admin", 
      default=None,
      type="string")
      
    (options, args) = parser.parse_args()
    
    print "starting server on port: ", options.serverport
      
    server = ContinuosTestWebServer(options.serverport)
    
    if options.admin_email:
        server.tracker.admin_email_address = options.admin_email
    else:
        parser.error("Must set the admin e-mail address using -a")
        
    server.tracker.start()
    server.start()

    

        

        


    
    
