import sys, os, re, subprocess

from distutils.core import Command
from distutils import log

try:
    from nose.core import TestProgram
except ImportError:
    TestProgram = None
    

class run_tests(Command):

    description = "run add unit tests"

    user_options = [
        ('test-dir', 't', "root directory of the tests cases")
    ]

    boolean_options = []


    def initialize_options (self):
        self.test_dir = None

    def finalize_options (self):
        self.announce("test directory: {0}".format(self.test_dir), level = log.INFO)
        if not self.test_dir:
            self.test_dir = 'test'

    def run (self):
        if TestProgram is None:
            self.error("no nosetests framework found, please install python nose packge first!")
            return
            
        TestProgram(argv=[self.test_dir], exit=False) 
            
            
            
            
            
    
            
