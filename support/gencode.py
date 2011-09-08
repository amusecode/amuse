import sys
import os.path
import os

from optparse import OptionParser


def get_amuse_directory():
    filename_of_this_script = __file__
    directory_of_this_script = os.path.dirname(os.path.dirname(filename_of_this_script))
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
    
    or: %prog --mode=dir name_of_the_code
    
    This script will generate code from the class with name <name_of_class_in_module>. The
    class must be defined in the module <name_of_module>. The module name
    can be a python file or the python module name.
    
    If mode is dir the script will create a directory with all files
    needed to start creating a code interface.
    
    This script handles all code generation for the AMUSE framework. It can
    be used to create C++ or Fortran code to handle the MPI messages, 
    create a header file or create stub code as a start for defining
    the interface between the code and AMUSE.
    
    Examples
    --------
    To generate code for interfacing with MPI do:
        %prog --type=c --mode=mpi test.py TestInterface
    or (for fortran):
        %prog --type=f90 --mode=mpi test.py TestInterface
        
    To generate a header file do (for C):
        %prog --type=h test.py TestInterface
    or (for C++):
        %prog --type=H test.py TestInterface
        
    To generate a stub file do:
        %prog --type=c --mode=stub test.py TestInterface
    or (for fortran):
        %prog --type=f90 --mode=stub test.py TestInterface
        
    To generate create a directory and put files in it do:
        %prog --type=c --mode=dir MyCode
    or (for fortran):
        %prog --type=f90 --mode=dir MyCode
    
    To see a description of all arguments do:
        %prog --help
    
    """
    
    def __init__(self):
        self.parser = OptionParser(self.usage)
        self.parser.prog = 'build.py' #hack to set the name, for reporting errors and help
        
        self.parser.add_option(
            "-t",
            "--type",
            choices=["c","h", "H", "f90"],
            default="c",
            dest="type",
            help="TYPE of the code to generate. Can be one of c, h, H, f90. <c> will generate c code. <h/H> will generate c/c++ header. <f90> will generate fortran 90 code. (Defaults to c)")
        self.parser.add_option(
            "-m",
            "--mode",
            choices=["mpi","stub", "dir", "sockets"],
            default="mpi",
            dest="mode",
            help="MODE of the code to generate. Can be <mpi>, <stub>, <dir> or <socket>. Generate the MPI handling code or STUB code for the link between mpi and the code (if needed). <dir> will create a directory ann populate it with the files needed to build a code. (Defaults to mpi)")
        self.parser.add_option(
            "-o",
            "--output",
            default="-",
            dest="output",
            help="Name of the OUTPUT file. Use - for standard out. ")
        self.parser.add_option(
            "-i",
            "--ignore",
            default="",
            dest="ignore",
            help="Name of the classes to ignore, functions defined on these classes will not generate code. Comma separated list")
        
        
        self.options = None
        self.arguments = None
        
    def parse_options(self):
        (self.options, self.arguments) = self.parser.parse_args()
        if self.options.ignore:
            self.options.ignore_classes = list(self.parse_ignore_classe())
        else:
            self.options.ignore_classes = []
        
    def parse_arguments(self):
        if self.options.mode == 'dir':
            if len(self.arguments) != 1:
                self.show_error_and_exit("incorrect number of arguments, need name of the code")
                
            self.options.name_of_the_code = self.arguments[0]
        else:
            if len(self.arguments) != 2:
                self.show_error_and_exit("incorrect number of arguments")
            try:
                self.options.name_of_module_or_python_file = self.arguments[0]
                self.options.name_of_class = self.arguments[1]
            except Exception as exception:
                self.show_error_and_exit(exception)
    
    def parse_ignore_classe(self):
        names = self.options.ignore.split(',')
        for name in names:
            index_of_module_classname_split = name.rfind('.')
            modulename = name[:index_of_module_classname_split]
            classname = name[index_of_module_classname_split+1:]
            
            __import__(modulename)
            class_to_ignore = getattr(sys.modules[modulename], classname)
            yield class_to_ignore
        
    
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
    result = create_c.GenerateACHeaderStringFromASpecificationClass()
    result.make_extern_c = False
    return result
   
def make_file(settings):
    
    try:
        if settings.name_of_module_or_python_file.endswith('.py'):
            module = {}
            execfile(settings.name_of_module_or_python_file, module)
            specification_class = module[settings.name_of_class]
        else:
            module = __import__(settings.name_of_module,fromlist=[settings.name_of_class])
            specification_class = getattr(module, settings.name_of_class)
    except ImportError as exception:
        uc.show_error_and_exit(exception)
    
        
    usecases = { 
        ('c','mpi'): create_c.GenerateACSourcecodeStringFromASpecificationClass,
        ('h','mpi'): create_c.GenerateACHeaderStringFromASpecificationClass,
        ('H','mpi'): make_cplusplus_header,
        ('f90','mpi'): create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass,      
        ('c','stub'): create_c.GenerateACStubStringFromASpecificationClass,    
        ('f90','stub'): create_fortran.GenerateAFortranStubStringFromASpecificationClass,
        ('c','sockets'): create_c_sockets.GenerateACSourcecodeStringFromASpecificationClass,    
        ('f90','sockets'): create_fortran_sockets.GenerateAFortranSourcecodeStringFromASpecificationClass,   
    }
    
    try:
        builder = usecases[(settings.type, settings.mode)]()
        builder.specification_class = specification_class
        builder.ignore_functions_from_specification_class = settings.ignore_classes
    except:
        uc.show_error_and_exit("'{0}' and '{1}' is not a valid combination of type and mode, cannot generate the code".format(settings.type, settings.mode))
    
    if settings.output == '-':
        print builder.result
    else:
        try:
            with open(settings.output, "w") as f:
                f.write(builder.result)
        except Exception as exception:
            uc.show_error_and_exit(exception)
            
            

def make_directory(settings):

    usecases = {
        ('c','dir'): create_dir.CreateADirectoryAndPopulateItWithFilesForACCode,    
        ('f90','dir'): create_dir.CreateADirectoryAndPopulateItWithFilesForAFortranCode, 
    }
    
    try:
        builder = usecases[(settings.type, settings.mode)]()
        builder.name_of_the_code_interface_class = settings.name_of_the_code
        builder.path_of_the_root_directory = os.getcwd()
    except:
        uc.show_error_and_exit("'{0}' and '{1}' is not a valid combination of type and mode, cannot generate the code".format(settings.type, settings.mode))
        
    builder.start()
    
    
if __name__ == '__main__':
    
    setup_sys_path()
    
    from amuse.rfi.tools import create_c
    from amuse.rfi.tools import create_fortran
    from amuse.rfi.tools import create_dir
    from amuse.rfi.tools import create_c_sockets
    from amuse.rfi.tools import create_fortran_sockets
    
    uc = ParseCommandLine()
    uc.start()
    
    settings = uc.options
    if settings.mode == 'dir':
        make_directory(settings)
    else:
        make_file(settings)
    
