import sys, os, re, subprocess

from distutils.core import Command
from distutils.command import config as _original_config
from distutils import log


class config(_original_config):
    
    description = "find AMUSE prerequisite libraries"

    user_options.extend([
        ('cuda-dir=', None,
         "root directory where cuda is installed, for example /opt/cuda"),
        ('cuda-sdk-dir=', None,
         "root directory where the cuda sdk is installed, for example /opt/cuda-sdk"),
        ('compile', None,
         "compile a test program to see if library files are correct"),
        ])

    def initialize_options(self):
        _original_config.initialize_options(self)
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
        result = set([])
        for x in self.default_search_dirs():
            result.append( os.path.join(x, 'lib'))
            result.append( os.path.join(x, 'lib64'))

    def python_root_dir(self):
        str = sys.executable
        index = str.find('bin')
        return str[:index]
    def finalize_options(self):
        _original_config.finalize_options(self)
        pass
    def run(self):
        config = ConfigParser.RawConfigParser()
    