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
        
        if measured_timestamp < 0:
            self.container.remove(self)
            monitor.errored(self)
            return
            
        if self.timestamp < measured_timestamp:
            self.timestamp = measured_timestamp
            monitor.updated(self)
            return
        
        monitor.unchanged(self)
    
    def get_last_modification_time(self):
        try:
            statinfo = os.stat(self.path)
            return statinfo.st_mtime
        except:
            return -1

    def walk(self, monitor):
        monitor.found(self)
    
        
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
            if os.path.islink(path):
                continue
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
    
    def walk(self, monitor):
        monitor.found(self)
         
        for x in self.elements:
            x.walk(monitor) 
            
        
            
class MonitorDirectories(object):
    def __init__(self, paths):
        self.elements = map(lambda x : MonitoredDirectory(x), paths)
        self.changed = False
        
    def check(self):
        self.changed = False
        for x in self.elements:
            x.check(self)
    
    def deleted(self, monitored_element):
        if not self.must_monitor_file(monitored_element):
            return
        self.changed = True
        
    def created(self, monitored_element):
        if not self.must_monitor_file(monitored_element):
            return
        self.changed = True
        
    def unchanged(self, monitored_element):
        pass
        
    def errored(self, monitored_element):
        print "error while monitoring file: ", monitored_element.path
        pass
    
    def updated(self, monitored_element):
        if not self.must_monitor_file(monitored_element):
            return
        self.changed = True

    def walk(self, callback_function):
        self.callback_function = callback_function
        for x in self.elements:
            x.walk(self)
            
    def found(self, monitored_element):
        if not self.must_monitor_file(monitored_element):
            return
        self.callback_function(monitored_element)
        
    def must_monitor_file(self, monitored_element):
        return (
            monitored_element.path.endswith('.py') and
            not os.path.basename(monitored_element.path).startswith('.')
        )
