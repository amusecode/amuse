from amuse.support.core import late, print_out

from amuse.support.options import option
from amuse.support.options import OptionalAttributes
from amuse.support import get_amuse_root_dir

import os
import inspect
import sys
        


class CreateAPythonWorker(OptionalAttributes):
    
    @option(sections=['data'])
    def amuse_root_dir(self):
        return get_amuse_root_dir()
        
    @late
    def channel_type(self):
        return 'mpi'
        
    @late
    def template_dir(self):
        return os.path.dirname(__file__)
        
    @late
    def worker_dir(self):
        return os.path.abspath(os.path.curdir)
        
    @late
    def template_string(self):
        path = self.template_dir
        path = os.path.join(path, 'python_code_script.template')
            
        with open(path, "r") as f:
            template_string = f.read()
        
        return template_string
    
    
    @late
    def worker_name(self):
        filename = os.path.basename(inspect.getfile(self.implementation_factory))
        filename = filename.split('.')[0]
        filename.replace(os.sep, '_')
        path = os.path.join(self.worker_dir, filename)
        
        return path
        
    @late
    def output_name(self):
        executable_path = self.worker_name
        return executable_path
    
    @late
    def interface_class(self):
        return self.specification_class
        
    def new_executable_script_string(self):
        return self.template_string.format(
            executable = sys.executable,
            syspath = ','.join(map(repr, sys.path)),
            factory_module = inspect.getmodule(self.implementation_factory).__name__,
            factory = self.implementation_factory.__name__,
            interface_module = inspect.getmodule(self.interface_class).__name__,
            interface = self.interface_class.__name__,
        )
    
    @property
    def result(self):
        return self.new_executable_script_string()
        
    def start(self):
        string = self.new_executable_script_string()
            
        with open(self.output_name, 'w') as f:
            f.write(string)
            
        os.chmod(self.output_name, 0o777)
        
