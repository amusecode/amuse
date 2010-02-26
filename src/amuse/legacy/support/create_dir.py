from amuse.support.core import late, print_out
import os

class CreateADirectoryAndPopulateItWithFilesForALegacyCode(object):
    
    @late
    def path_of_the_root_directory(self):
        return os.path.dirname(os.path.dirname(__file__))
    @late
    def name_of_the_legacy_code(self):
        return self.name_of_the_code_interface_class.lower()
        
    @late
    def name_of_the_python_module(self):
        return 'interface'
        
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
    def path_of_amuse(self):
        current = os.path.dirname(os.path.dirname(__file__))
        while not os.path.exists(os.join(current, 'build.py')):
            current = os.path.dirname(current)
        return current
        
    @late
    def reference_to_amuse_path(self):
        return os.relpath(self.path_of_amuse, self.path_of_the_root_directory)
        
    def start(self):
        
        self.make_directories()
        self.make_python_files()
        self.make_makefile()
        
    def make_directories(self):
        os.mkdir(self.path_of_the_legacy_code)
        os.mkdir(self.path_of_the_source_code)
        
    def make_python_files(self):
        pass
        
    def make_makefile(self):
        pass
        
        
