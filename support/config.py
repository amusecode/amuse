import sys
import os
import re
import subprocess
import ConfigParser

from distutils.core import Command
from distutils.command.config import config as original_config
from distutils import log
from distutils import spawn


class config(original_config):
    
    description = "find AMUSE prerequisite libraries"

    user_options = original_config.user_options + [
        ('cuda-dir=', None,
         "root directory where cuda is installed, for example /opt/cuda"),
        ('cuda-sdk-dir=', None,
         "root directory where the cuda sdk is installed, for example /opt/cuda-sdk"),
        ('compile', None,
         "compile a test program to see if library files are correct"),
        ]

    def initialize_options(self):
        original_config.initialize_options(self)
        self.library_dirs = []
        self.include_dirs = []
        
    def default_search_dirs(self):
        if os.name == "posix":
            dir = self.python_root_dir()
            return set(['/usr', '/usr/local', '/opt', '/opt/local', dir])
        else:
            return set([])
            
    def default_include_dirs(self):
        return map(lambda x : os.path.join(x, 'include'), self.default_search_dirs())
        
    def default_lib_dirs(self):
        result = []
        for x in self.default_search_dirs():
            result.append( os.path.join(x, 'lib'))
            result.append( os.path.join(x, 'lib64'))
        return set(result)
        
    def python_root_dir(self):
        str = sys.executable
        index = str.find('bin')
        return str[:index]
        
    def finalize_options(self):
        original_config.finalize_options(self)
        pass
    
    def find_mpich2_directory(self):
        for dir in self.default_lib_dirs():
            name = os.path.join(dir, 'libmpich.so')
            if os.path.exists(name):
                return dir
        log.error("Could not find the MPICH-2 library file)")
        
    def add_mpich2_options(self, config):
        mpich2_dir = self.find_mpich2_directory()
        config.set('tests','mpi-lib-dir',mpich2_dir)
    def find_cuda(self):
        dir = spawn.find_executable('nvcc')
        print os.path.dirname(os.path.dirname(dir))
        print dir
        
    
    def run(self):
        config = ConfigParser.RawConfigParser()
        log.info("configuring AMUSE")
        config.add_section('tests')
        config.set('tests', 'test-dir', '15')
        self.add_mpich2_options(config)
        self.find_cuda()
        file = open('setup.cfg', 'wb')
        config.write(file)
        file.close()
    
    
        
    