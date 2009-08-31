__revision__ = "$Id:$"

import sys, os, re, subprocess
from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.util import convert_path
from distutils import log

# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')

class LegacyCommand(Command):
    user_options = [
        ('legacy-dir', 'd', "directory containing legacy codes"),
        ]

    boolean_options = ['force']


    def initialize_options (self):
        self.legacy_dir = None
        self.amuse_src_dir =  os.path.join('src','amuse')

    def finalize_options (self):
        if self.legacy_dir is None:
            self.legacy_dir = os.path.join(self.amuse_src_dir,'legacy')

    def subdirs_in_legacy_dir(self):
        names = os.listdir(self.legacy_dir)
        for name in names:
            path = os.path.join(self.legacy_dir, name)
            if os.path.isdir(path):
                yield path
                
    def makefile_paths(self):
        for x in self.subdirs_in_legacy_dir():
            for name in ('makefile', 'Makefile'):
                makefile_path = os.path.join(x, name)
                if os.path.exists(makefile_path):
                    yield x
    

class BuildLegacy(LegacyCommand):

    description = "build interfaces to legacy codes"
    
    def run (self):
        for x in self.makefile_paths():
            self.announce("building " + x)
            self.spawn(['make','-C', x, 'all'])
 
 
class CleanLegacy(LegacyCommand):

    description = "clean build products in legacy codes"

    def run (self):
        for x in self.makefile_paths():
            self.announce("cleaning " + x)
            self.spawn(['make','-C', x, 'clean'])
            
            
            
            
            
    
            
           
            
            
            
            
    
            
