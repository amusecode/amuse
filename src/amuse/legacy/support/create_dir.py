from amuse.support.core import late, print_out

from amuse.legacy import *

import os

interface_file_template = """\
from amuse.legacy import *

class {0.name_of_the_legacy_interface_class}({0.name_of_the_superclass_for_the_legacy_interface_class}):
    
    include_headers = []
    
    def __init__(self, **keyword_arguments):
        {0.name_of_the_superclass_for_the_legacy_interface_class}.__init__(self, **keyword_arguments)
    
class {0.name_of_the_code_interface_class}({0.name_of_the_superclass_for_the_code_interface_class}):

    def __init__(self):
        {0.name_of_the_superclass_for_the_code_interface_class}.__init__(self,  {0.name_of_the_legacy_interface_class}())
    
    
    
"""

class CreateADirectoryAndPopulateItWithFilesForALegacyCode(object):
    
    @late
    def path_of_the_root_directory(self):
        return os.path.dirname(os.path.dirname(__file__))
        
    @late
    def name_of_the_legacy_code(self):
        return self.name_of_the_code_interface_class.lower()
        
    @late
    def name_of_the_python_module(self):
        return 'interface.py'
        
    @late
    def name_of_the_code_interface_class(self):
        return 'MyCode'
        
    @late
    def name_of_the_legacy_interface_class(self):
        return self.name_of_the_code_interface_class + 'Interface'
    
    @late
    def name_of_the_code_directory(self):
        return 'src'
    
    @late
    def path_of_the_legacy_code(self):
        return os.path.join(self.path_of_the_root_directory, self.name_of_the_legacy_code)
        
    @late
    def path_of_the_source_code(self):
        return os.path.join(self.path_of_the_legacy_code, self.name_of_the_code_directory)
        
    @late
    def path_of_the_init_file(self):
        return os.path.join(self.path_of_the_legacy_code, '__init__.py')
        
    @late
    def path_of_the_interface_file(self):
        return os.path.join(self.path_of_the_legacy_code, self.name_of_the_python_module)
        
    @late
    def path_of_amuse(self):
        current = os.path.dirname(os.path.dirname(__file__))
        while not os.path.exists(os.join(current, 'build.py')):
            current = os.path.dirname(current)
        return current
        
    @late
    def reference_to_amuse_path(self):
        return os.relpath(self.path_of_amuse, self.path_of_the_root_directory)
        
    @late 
    def name_of_the_superclass_for_the_legacy_interface_class(self):
        return LegacyInterface.__name__
        
    @late
    def name_of_the_superclass_for_the_code_interface_class(self):
        return CodeInterface.__name__
        
    def start(self):
        
        self.make_directories()
        self.make_python_files()
        self.make_makefile()
        
    def make_directories(self):
        os.mkdir(self.path_of_the_legacy_code)
        os.mkdir(self.path_of_the_source_code)
        
    def make_python_files(self):
        with open(self.path_of_the_init_file, "w") as f:
            f.write("# generated file")
        
        with open(self.path_of_the_interface_file, "w") as f:
            string = interface_file_template.format(self)
            f.write(string)
        
    def make_makefile(self):
        pass
        
        
