import sys
import os.path
import os

from optparse import OptionParser

# Should probably use an absolute import here (support.config), but
# we're not guaranteed this script will always be in a support
# subdirectory with an __init__.py file.
#try:  # running as a module
#    from . import config
#except (ImportError, ValueError):  # running as a stand-alone script
#    import config


if sys.hexversion > 0x03000000:
    def get_amuse_directory():
        filename_of_this_script = __file__
        directory_of_this_script = os.path.dirname(os.path.dirname(filename_of_this_script))
        directory_of_this_script = os.path.join(directory_of_this_script, 'build', 'lib')
        if os.path.isabs(directory_of_this_script):
            return directory_of_this_script
        else:
            return os.path.abspath(directory_of_this_script)
else:
    
    def get_amuse_directory():
        filename_of_this_script = __file__
        directory_of_this_script = os.path.dirname(os.path.dirname(filename_of_this_script))
        if os.path.isabs(directory_of_this_script):
            return directory_of_this_script
        else:
            return os.path.abspath(directory_of_this_script)

def get_amuse_directory_root():
    filename_of_this_script = __file__
    directory_of_this_script = os.path.dirname(os.path.dirname(filename_of_this_script))
    if os.path.isabs(directory_of_this_script):
        return directory_of_this_script
    else:
        return os.path.abspath(directory_of_this_script)

def setup_sys_path():
    amuse_directory = os.environ["AMUSE_DIR"]
    sys.path.insert(0, amuse_directory)
    try:
      src_root_directory = os.environ["MUSE_PACKAGE_DIR"]
      sys.path.insert(0, src_root_directory)
    except:
      src_root_directory = amuse_directory
    sys.path.insert(0, os.path.join(src_root_directory,"src"))
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
            choices=["c","h", "H", "f90", "py", "java", "cython"],
            default="c",
            dest="type",
            help="TYPE of the code to generate. Can be one of c, h, H, f90, py, java or cython. <c> will generate c code. <h/H> will generate c/c++ header. <f90> will generate fortran 90 code. <py> will generate a python worker wrapper <java> will generate java interface or class, depending on mode. (Defaults to c)")
        
        self.parser.add_option(
            "-m",
            "--mode",
            choices=["mpi","stub", "dir", "sockets", "interface", "class", "script"],
            default="mpi",
            dest="mode",
            help="MODE of the code to generate. Can be <mpi>, <stub>, <dir>,<sockets>, <interface>, <class> or <script>. Generate the MPI handling code or STUB code for the link between mpi and the code (if needed). <dir> will create a directory ann populate it with the files needed to build a code. (Defaults to mpi)")
        
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
    
        self.parser.add_option(
            "-u",
            "--underscore",
            default="",
            dest="underscore",
            help="Name of the classes to underscore the functions of, for XL fortran compilers")
        
        self.parser.add_option(
            "-n",
            "--needs-mpi",
            default="true",
            dest="needs_mpi",
            help="If this boolean flag is set, the worker will initialize mpi, even in the sockets channel is used. Defaults to true")
        
        self.parser.add_option(
            "-x",
            "--executable",
            action="store_true", 
            default=False,
            dest="make_executable",
            help="Set the executable bit when generating the output file")
        
        self.parser.add_option(
            "--cython-import",
            default="",
            dest="cython_import",
            help="Name of the module to import for the cython worker (name of the .so file)")
        
        self.parser.add_option(
            "--prefix",
            default="",
            dest="function_name_prefix",
            help="Prefix for generated function names, relevant for cython")

        self.parser.add_option(
            "--get_amuse_dir",
            action="store_true",
            default=False,
            dest="get_amuse_dir",
            help="Only output amuse directory")
        
        self.options = None
        self.arguments = None
        
    def parse_options(self):
        (self.options, self.arguments) = self.parser.parse_args()
        if self.options.ignore:
            self.options.ignore_classes = list(self.parse_ignore_classes())
        else:
            self.options.ignore_classes = []
            
        if self.options.underscore:
            self.options.underscore_classes = list(self.parse_underscore_classes())
        else:
            self.options.underscore_classes = []
        
        self.options.name_of_implementation_class = None 
        self.options.name_of_module_or_python_file = None 
        self.options.name_of_class = None 
        self.options.name_of_the_code = None 
        
        
    def parse_arguments(self):
        if self.options.get_amuse_dir:
            return
        if self.options.mode == 'dir':
            if len(self.arguments) != 1:
                self.show_error_and_exit("incorrect number of arguments, need name of the code")
                
            self.options.name_of_the_code = self.arguments[0]
        else:
            if not len(self.arguments) in (2,3) :
                self.show_error_and_exit("incorrect number of arguments")
            try:
                self.options.name_of_module_or_python_file = self.arguments[0]
                if len(self.arguments) > 1:
                    self.options.name_of_class = self.arguments[1]
                if len(self.arguments) > 2:
                    self.options.name_of_implementation_class = self.arguments[2]
            except Exception as exception:
                self.show_error_and_exit(exception)
    

    def parse_ignore_classes(self):
        names = self.options.ignore.split(',')
        for name in names:
            index_of_module_classname_split = name.rfind('.')
            modulename = name[:index_of_module_classname_split]
            classname = name[index_of_module_classname_split+1:]
            
            __import__(modulename)
            class_to_ignore = getattr(sys.modules[modulename], classname)
            yield class_to_ignore
        
    def parse_underscore_classes(self):
        names = self.options.underscore.split(',')
        for name in names:
            index_of_module_classname_split = name.rfind('.')
            modulename = name[:index_of_module_classname_split]
            classname = name[index_of_module_classname_split+1:]
            
            __import__(modulename)
            class_to_underscore = getattr(sys.modules[modulename], classname)
            yield class_to_underscore
        
    
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

def make_a_python_worker(channel_type):
    result = create_python_worker.CreateAPythonWorker()
    result.channel_type = channel_type
    return result

def make_a_mpi_python_worker():
    return make_a_python_worker('mpi')
    
def make_a_socket_python_worker():
    return make_a_python_worker('sockets')
    
def make_file(settings):
    implementation_class = None
    try:
        if settings.name_of_module_or_python_file.endswith('.py'):
            module = {}
            # Replace with runpy in the future?
            with open(settings.name_of_module_or_python_file) as fh:
                text = fh.read()
                code = compile(text, settings.name_of_module_or_python_file, 'exec')
                exec(code, module)
            #execfile(settings.name_of_module_or_python_file, module)
            specification_class = module[settings.name_of_class]
            if not settings.name_of_implementation_class is None:
                implementation_class = module[settings.name_of_implementation_class]
        else:
            module = __import__(settings.name_of_module_or_python_file,fromlist=[settings.name_of_class])
            specification_class = getattr(module, settings.name_of_class)
            if not settings.name_of_implementation_class is None:
                implementation_class = getattr(module, settings.name_of_implementation_class)
    except ImportError as exception:
        uc.show_error_and_exit(exception)
    
        
    usecases = { 
        ('c','mpi'): create_c.GenerateACSourcecodeStringFromASpecificationClass,
        ('h','mpi'): create_c.GenerateACHeaderStringFromASpecificationClass,
        ('H','mpi'): make_cplusplus_header,
        ('f90','mpi'): create_fortran.GenerateAFortranSourcecodeStringFromASpecificationClass,      
        ('c','stub'): create_c.GenerateACStubStringFromASpecificationClass,    
        ('f90','stub'): create_fortran.GenerateAFortranStubStringFromASpecificationClass,
        ('java','interface'): create_java.GenerateAJavaInterfaceStringFromASpecificationClass,
        ('java','class'): create_java.GenerateAJavaSourcecodeStringFromASpecificationClass,
        ('java','script'): create_java.GenerateAJavaWorkerScript,
        ('py','sockets'): make_a_socket_python_worker,
        ('py','mpi'): make_a_mpi_python_worker,    
        ('cython','script'): create_cython.GenerateACythonStartScriptStringFromASpecificationClass,
        ('cython','mpi'): create_cython.GenerateACythonSourcecodeStringFromASpecificationClass,
        ('cython', 'interface'): create_cython.GenerateAFortranInterfaceSourcecodeStringFromASpecificationClass
    }
    
    try:
        builder = usecases[(settings.type, settings.mode)]()
        builder.specification_class = specification_class
        
        if not implementation_class is None:
            builder.implementation_factory = implementation_class
            
        builder.ignore_functions_from_specification_classes = settings.ignore_classes
        builder.underscore_functions_from_specification_classes = settings.underscore_classes
        builder.needs_mpi = settings.needs_mpi.lower() == 'true'
        builder.is_mpi_enabled = config.mpi.is_enabled
        builder.name_of_outputfile = settings.output
        builder.cython_import = settings.cython_import
        builder.function_name_prefix = settings.function_name_prefix
    except:
        uc.show_error_and_exit("'{0}' and '{1}' is not a valid combination of type and mode, cannot generate the code".format(settings.type, settings.mode))
    
    if settings.output == '-':
        sys.stdout.write(str(builder.result) + '\n')
    else:
        try:
            
            with open(settings.output, "w") as f:
                f.write(builder.result)
                
            if settings.make_executable:
                os.chmod(settings.output, 0o755)
                
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
    
    # setup_sys_path()

    from amuse import config

    from amuse.rfi.tools import create_c
    from amuse.rfi.tools import create_fortran
    from amuse.rfi.tools import create_java
    from amuse.rfi.tools import create_dir
    from amuse.rfi.tools import create_python_worker
    from amuse.rfi.tools import create_cython
        
    from amuse.support import get_amuse_root_dir    
        
    uc = ParseCommandLine()
    uc.start()
    
    settings = uc.options
    if settings.get_amuse_dir:
        print(get_amuse_root_dir())
        exit(0)
    elif settings.mode == 'dir':
        make_directory(settings)
    else:
        make_file(settings)
    
