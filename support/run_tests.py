import sys, os, re, subprocess

from distutils.core import Command
from distutils import log

from nose.core import TestProgram

# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')

class run_tests(Command):

    description = "run add unit tests"

    user_options = [
        ('test-dir', 't', "root directory of the tests cases")
        ]

    boolean_options = []


    def initialize_options (self):
        self.test_dir = None

    def finalize_options (self):
        if not self.test_dir:
            self.test_dir = 'test'

    def run (self):
        TestProgram(argv=[self.test_dir], exit=False) 
            
            
            
            
            
    
            
