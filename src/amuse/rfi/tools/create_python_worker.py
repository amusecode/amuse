from amuse.support.core import late, print_out

from amuse.support.options import option
from amuse.support.options import OptionalAttributes

import os
import inspect
import sys
        


class CreateAPythonWorker(OptionalAttributes):
    
    @option(sections=['data'])
    def amuse_root_dir(self):
        if 'AMUSE_DIR' in os.environ:
            return os.environ['AMUSE_DIR']    
        previous = None
        result = os.path.abspath(__file__)
        while not os.path.exists(os.path.join(result,'build.py')):
            result = os.path.dirname(result)
            if result == previous:
                return os.path.dirname(os.path.dirname(__file__))
            previous = result
        return result
        
    @late
    def channel_type(self):
        return 'mpi'
        
    @late
    def template_dir(self):
        return os.path.dirname(__file__)
        
    @late
    def template_string(self):
        path = self.template_dir
        if self.channel_type in ('sockets', 'ibis'):
            path = os.path.join(path, 'python_socket_code_script.template')
        else:
            path = os.path.join(path, 'python_code_script.template')
            
        with open(path, "r") as f:
            template_string = f.read()
        
        return template_string
    
    
    @late
    def worker_name(self):
        filename = os.path.basename(inspect.getfile(self.implementation_factory))
        filename = filename.split('.')[0]
        filename.replace(os.sep, '_')
        path = os.path.abspath(os.path.curdir)
        path = os.path.join(path, filename)
        
        return path
        
    @late
    def output_name(self):
        executable_path = self.worker_name
        if self.channel_type in ('sockets', 'ibis'):
            executable_path += '_sockets'
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
            
        os.chmod(self.output_name, 0777)
        
