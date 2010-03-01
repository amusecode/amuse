#!/usr/bin/env python

import sys
import os.path
import os

from optparse import OptionParser


def get_amuse_directory():
    filename_of_this_script = __file__
    directory_of_this_script = os.path.dirname(filename_of_this_script)
    if os.path.isabs(directory_of_this_script):
        return directory_of_this_script
    else:
        return os.path.abspath(directory_of_this_script)

def setup_sys_path():
    amuse_directory = get_amuse_directory()
    
    sys.path.insert(0, amuse_directory)
    sys.path.insert(0, os.path.join(amuse_directory,"src"))
    
    sys.path.append(os.getcwd())

class ParseCommandLine(object):
    usage = """usage: %prog [options] name_of_module name_of_class_in_module.
    
    Will generate code from the class with name <name_of_class_in_module>. The
    class must be defined in the module <name_of_module>. The module name
    can be a python file or the python module name.
    
    """
    
    def __init__(self):
        self.parser = OptionParser(self.usage)
        self.parser.add_option(
            "-t",
            "--type",
            choices=["c","h", "H", "f90"],
            default="c",
            dest="type",
            help="TYPE of the code to generate. Can be one of c, h, H or f90. <c> will generate c code. <h/H> will generate c/c++ header. <f90> will generate fortran 90 code. (Defaults to c)")
        self.parser.add_option(
            "-m",
            "--mode",
            choices=["mpi","interface"],
            default="mpi",
            dest="mode",
            help="MODE of the code to generate. Can be MPI or INTERFACE. Generate the MPI handling code or the INTERFACE code. Only needed when generating code. (Defaults to mpi)")
        self.parser.add_option(
            "-o",
            "--output",
            default="-",
            dest="output",
            help="Name of the OUTPUT file. Use - for standard out. ")
        
        
        self.options = None
        self.arguments = None
        
    def parse_options(self):
        (self.options, self.arguments) = self.parser.parse_args()
        
    def parse_arguments(self):
        if len(self.arguments) != 2:
            self.parser.error("incorrect number of arguments")
        try:
            self.options.name_of_module_or_python_file = self.arguments[0]
            self.options.name_of_class = self.arguments[1]
        except Exception as exception:
            self.show_error_and_exit(exception)
            
    
    def start(self):
        self.parse_options()
        self.parse_arguments()
        
    def show_error_and_exit(self, exception):
        self.parser.error(exception)
        
    
    
    
    
def module_name(string):
    if string.endswith('.py'):
        amuse_src_directory = os.path.join(get_amuse_directory(), 'src')
        if not os.path.isabs(string):
            string = os.path.join(os.getcwd(), string)
        if not os.path.exists(string):
            raise Exception("Cannot find file with name {0}".format(string))
        if not string.startswith(amuse_src_directory):
            raise Exception("File {0} must be placed under directory {1}.".format(string, amuse_src_directory))
        
        string = string[len(amuse_src_directory)+1:]
        string = string[:-len('.py')]
        string = string.replace(os.sep, '.')
    return string
    
def make_cplusplus_header():
    result = create_c.MakeACHeaderStringOfAClassWithLegacyFunctions()
    result.make_extern_c = False
    return result
    
if __name__ == "__main__":
    setup_sys_path()
    
    from amuse.legacy.support import create_c
    from amuse.legacy.support import create_fortran
    
    uc = ParseCommandLine()
    uc.start()
    
    settings = uc.options
    try:
        if settings.name_of_module_or_python_file.endswith('.py'):
            module = {}
            execfile(settings.name_of_module_or_python_file, module)
            class_with_legacy_functions = module[settings.name_of_class]
        else:
            module = __import__(settings.name_of_module,fromlist=[settings.name_of_class])
            class_with_legacy_functions = getattr(module, settings.name_of_class)
    except ImportError as exception:
        uc.show_error_and_exit(exception)
        
        
        
    
    usecases = { 
        ('c','mpi'): create_c.MakeACStringOfAClassWithLegacyFunctions,
        ('h','mpi'): create_c.MakeACHeaderStringOfAClassWithLegacyFunctions,
        ('H','mpi'): make_cplusplus_header,
        ('f90','mpi'): create_fortran.MakeAFortranStringOfAClassWithLegacyFunctions,        
    }
    
    builder = usecases[(settings.type, settings.mode)]()
    builder.class_with_legacy_functions = class_with_legacy_functions
    if settings.output == '-':
        print builder.result
    else:
        try:
            with open(settings.output, "w") as f:
                f.write(builder.result)
        except Exception as exception:
            uc.show_error_and_exit(exception)

    
    
    
        
    
