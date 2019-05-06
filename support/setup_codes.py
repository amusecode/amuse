from __future__ import print_function

__revision__ = "$Id:$"

import sys, os, re, subprocess

from . import supportrc

try:
    import numpy
except ImportError:
    print( "numpy etc needed during build; operation may fail" )

try:
    import ConfigParser as configparser
    from StringIO import StringIO
except ImportError:
    import configparser
    from io import StringIO

import os.path
import datetime
import stat

from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.dir_util import create_tree
from distutils.util import convert_path
from distutils import log
from distutils import spawn
from distutils import file_util
from distutils.errors import DistutilsError
if sys.hexversion > 0x03000000:
    from distutils.util import run_2to3
    from distutils.command.build_py import build_py_2to3
from distutils.command.build import build
from distutils.command.clean import clean
from distutils.command.install import install

from subprocess import call, Popen, PIPE, STDOUT

if supportrc["framework_install"]:
    from .generate_main import generate_main
    from .build_latex import build_latex
    from .run_tests import run_tests
  
try:
    from numpy.distutils import fcompiler
except ImportError:
    fcompiler = None

# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')
    
from glob import glob

def pyfiles_in_build_dir(builddir):
    module_files = glob(os.path.join(builddir, "*.py"))
    result = []
    for x in module_files:
        result.append(os.path.abspath(x))
    return result

def run_2to3_on_build_dirs(paths, target, src):
    for dir in paths:
        buildir = os.path.join(target,  os.path.relpath(dir, src))
        run_2to3(pyfiles_in_build_dir(buildir))

class InstallLibraries(Command):
    user_options = [
        ('build-temp=', 't',
         "directory for temporary files (build by-products)"),
        ('inplace', 'i',
         "ignore build-lib and put compiled extensions into the source " +
         "directory alongside your pure Python modules", 1),
        ('no-inplace', 'k',
         "put compiled extensions into the  build temp "),
        ('lib-dir=', 'l', "directory containing libraries to build"),
        ('install-data=', None, "installation directory for data files"),
        ('root=', None, "install everything relative to this alternate root directory"),
    ]

    negative_opt = {'no-inplace':'inplace'}
    
    boolean_options = ['inplace']

    def initialize_options (self):
        self.codes_dir = None
        self.lib_dir = None
        self.inplace = False
        self.build_lib = None
        self.build_temp = None
        self.install_data = None
        self.root=None
                
    def finalize_options (self):
        self.set_undefined_options(
            'install',
           ('build_lib', 'build_lib'),
           ('root', 'root'),
           ('install_data', 'install_data'),           
        )

        self.set_undefined_options(
            'build',
           ('build_temp', 'build_temp'),
        )

        if self.lib_dir is None:
            if self.inplace:
                self.lib_dir = os.path.join('lib')
            else:
                self.lib_dir = os.path.join(self.build_temp, 'lib')
        else:
            if self.inplace:
                pass
            else:
                self.lib_dir=os.path.join(self.build_temp, 'lib')

    def run(self):
        data_dir = os.path.join(self.install_data,'share','amuse') # for the moment add to amuse..
        data_dir = os.path.abspath(data_dir)

        # copy only:
        # '*.h', '*.a', '*.mod', '*.inc', '*.so', '*.dylib'
        files=[os.path.join(dp, f) for dp, dn, fn in os.walk(self.lib_dir) for f in fn]
        ext=['.h', '.a', '.mod', '.inc', '.so', '.dylib']
        files=[f  for f in files if (os.path.splitext(f)[1] in ext)]
        files=[os.path.relpath(f,self.lib_dir) for f in files]
        create_tree(os.path.join(data_dir,'lib'), files)

        for f in files:
            src=os.path.join(self.lib_dir,f)
            target=os.path.join(data_dir,'lib', f)
            self.copy_file(src,target)

class GenerateInstallIni(Command):
    user_options =   (
        ('build-dir=', 'd', "directory to install to"),
        ('install-data=', None, "installation directory for data files"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ('root=', None, "install everything relative to this alternate root directory"),
    )
    
    boolean_options = ['force']
    
    def initialize_options(self):
        self.build_dir = None
        self.install_data = None
        self.force = False
        self.root = None
        
    def finalize_options(self):
        self.set_undefined_options(
            'install',
            ('build_lib', 'build_dir'),
            ('install_data', 'install_data'),
            ('root', 'root'),
            ('force', 'force'),
        )
        #~ print(self.install_data)
        #~ raise
        
    def run(self):
        outfilename = os.path.join(self.build_dir, supportrc["package_name"], 'amuserc')
        
        
        # this does not work for pip installs
        #~ data_dir = os.path.join(self.install_data,'share','amuse')
        #~ if not self.root is None:
            #~ data_dir = os.path.relpath(data_dir,self.root)
            #~ data_dir =  os.path.join('/',data_dir)
        #~ else:
            #~ data_dir = os.path.abspath(data_dir)

        installinilines = []
        installinilines.append('[channel]')
        installinilines.append('must_check_if_worker_is_up_to_date=0')
        installinilines.append('use_python_interpreter=1')
        #installinilines.append('worker_code_directory={0}'.format(os.path.join(data_dir, 'bin')))
        if sys.platform == 'win32':
            installinilines.append('worker_code_suffix=".exe"')
        installinilines.append('[data]')
        #~ installinilines.append('input_data_root_directory={0}'.format(os.path.join(data_dir, 'data')))
        installinilines.append('output_data_root_directory=_amuse_output_data')
        #~ installinilines.append('amuse_root_dir={0}'.format(data_dir))
        
        if 'BUILD_BINARY' in os.environ:
            installinilines.append('[test]')
            installinilines.append('can_run_tests_to_compile_modules=0')

        self.mkpath(os.path.join(self.build_dir, supportrc["package_name"]))
        file_util.write_file(outfilename, installinilines)
        
        


class CodeCommand(Command):
    user_options = [
        ('build-lib=', 'b',
         "directory for compiled extension modules"),
        ('build-temp=', 't',
         "directory for temporary files (build by-products)"),
        ('inplace', 'i',
         "ignore build-lib and put compiled extensions into the source " +
         "directory alongside your pure Python modules", 1),
        ('no-inplace', 'k',
         "put compiled extensions into the  build temp "),
        ('define=', 'D',
         "C preprocessor macros to define"),
        ('undef=', 'U',
         "C preprocessor macros to undefine"),
        ('debug', 'g',
         "compile/link with debugging information"),
        ('force', 'f',
         "forcibly build everything (ignore file timestamps)"),
        ('variant', 'V',
         "build variants of the codes (gpu versions etc)"),
        ('codes-dir=', 'd', "directory containing codes"),
        ('lib-dir=', 'l', "directory containing libraries to build"),
    ]

    negative_opt = {'no-inplace':'inplace'}
    
    boolean_options = ['force', 'inplace', 'debug', 'variant']

    def initialize_options (self):
        self.codes_dir = None
        self.lib_dir = None
        self.lib_src_dir = None
        self.amuse_src_dir =  os.path.join('src',supportrc["package_name"])
        self.environment = {}
        self.environment_notset = {}
        self.found_cuda = False
        self.found_sapporo = False
        self.variant = True
        self.inplace = False
        
        self.build_lib = None
        self.build_temp = None
        self.debug = None
        self.force = None


        
    def finalize_options (self):
        self.set_undefined_options(
            'build',
           ('build_lib', 'build_lib'),
           ('build_temp', 'build_temp'),
           ('debug', 'debug'),
           ('force', 'force'),
        )

        self.config=None

        if supportrc["framework_install"]:
            try:
                from . import config
                self.config=config
            except ImportError:
                # continue
                pass
        else:
            from amuse import config
            self.config=config

        if self.codes_dir is None:
            if self.inplace:
                self.codes_dir = os.path.join(self.amuse_src_dir,'community')
                self.codes_src_dir = self.codes_dir
            else:
                #~ self.codes_dir = os.path.join(self.build_temp, 'src', 'amuse', 'community')
                self.codes_dir = os.path.join(self.build_temp, 'codes')
                self.codes_src_dir = os.path.join(self.amuse_src_dir,'community')
        else:
            if self.inplace:
                self.codes_src_dir = self.codes_dir
            else:
                self.codes_src_dir = self.codes_dir
                #~ self.codes_dir=os.path.join(self.build_temp, 'src', 'amuse', 'community')
                self.codes_dir=os.path.join(self.build_temp, 'codes')
            
        if self.lib_dir is None:
            if self.inplace:
                self.lib_dir = os.path.join('lib')
                self.lib_src_dir = self.lib_dir
            else:
                self.lib_dir = os.path.join(self.build_temp, 'lib')
                self.lib_src_dir = os.path.join('lib')
        else:
            if self.inplace:
                self.lib_src_dir = self.codes_dir
            else:
                self.lib_src_dir = self.codes_dir
                self.lib_dir=os.path.join(self.build_temp, 'lib')

        if self.config:
            self.environment['PYTHON'] = self.config.interpreters.python
        else:
            self.environment['PYTHON'] = sys.executable
            
        self.set_cuda_variables()
        self.set_mpi_variables()
        self.set_compiler_variables()
        
        
        self.set_fortran_variables()
        
        if 'FORTRAN' in self.environment:
            self.environment['F90'] = self.environment['FORTRAN']
            self.environment['FC'] = self.environment['FORTRAN']
        self.set_java_variables()
        self.set_openmp_flags()
        self.set_libdir_variables()
        self.set_libs_variables()
        self.save_cfgfile_if_not_exists()
        
        if 'MSYSCON' in os.environ:
            pass
        else:
            if not supportrc["framework_install"]:
                try:
                    from amuse.support import get_amuse_root_dir
                except ImportError:
                    raise Exception("AMUSE framework needs to be installed and environment set up.")
                self.environment['AMUSE_DIR'] = get_amuse_root_dir()
            else:
                if self.inplace:
                   self.environment['AMUSE_DIR'] = os.path.abspath(os.getcwd())
                else:
                   self.environment['AMUSE_DIR'] = os.path.abspath(self.build_temp)

            if self.inplace:
               self.environment['MUSE_PACKAGE_DIR'] = os.path.abspath(os.getcwd())
            else:
               self.environment['MUSE_PACKAGE_DIR'] = os.path.abspath(self.build_temp)

    
    def set_fortran_variables(self):
        if 'FORTRAN' in self.environment:
            return
            
        if 'FORTRAN' in os.environ:
            self.environment['FORTRAN'] = os.environ['FORTRAN']
            return
            
        if self.config:
            self.environment['FORTRAN'] = self.config.compilers.fc
            return
        
        if 'FC' in os.environ:
            self.environment['FORTRAN'] = os.environ['FC']
            return
            
        if 'FORT' in os.environ:
            self.environment['FORTRAN'] = os.environ['FORT']
            return
            
        if 'F90' in os.environ:
            self.environment['FORTRAN'] = os.environ['F90']
            return
            
        mpif90 = os.environ['MPIF90'] if 'MPIF90' in os.environ else 'mpif90'
        
        try:
            process = Popen([mpif90,'-show'], stdout = PIPE, stderr = PIPE)
            stdoutstring, stderrstring = process.communicate()
            if process.returncode == 0:
                parts = stdoutstring.split()
                self.environment['FORTRAN']  = parts[0]
                return
            
            process = Popen([mpif90,'--showme '], stdout = PIPE, stderr = PIPE)
            stdoutstring, stderrstring = process.communicate()
            if process.returncode == 0:
                parts = stdoutstring.split()
                self.environment['FORTRAN']  = parts[0]
                return  
        except:
            pass
            
        
        if fcompiler:    
            compiler = fcompiler.new_fcompiler(requiref90=True)
            fortran_executable = compiler.executables['compiler_f90'][0]
            self.environment['FORTRAN'] = fortran_executable
    
    
    
    def is_mpi_enabled(self):
        if self.config and hasattr(self.config.mpi, 'is_enabled'):
            return self.config.mpi.is_enabled
        else:
            return True
    
    def set_cuda_variables(self):
        all_found = True
        if self.config and self.config.cuda.is_enabled:
            self.found_cuda = True
            self.environment['CUDA_LIBDIRS'] = '-L'+self.config.cuda.toolkit_path+'/lib' + ' -L'+self.config.cuda.toolkit_path+'/lib64'
            self.environment['CUDA_TK'] = self.config.cuda.toolkit_path
            self.environment['CUDA_SDK'] = self.config.cuda.sdk_path
            if hasattr(self.config.cuda, 'cuda_libs'):
                self.environment['CUDA_LIBS'] = self.config.cuda.cuda_libs
            else:
                raise DistutilsError("configuration is not up to date for cuda, please reconfigure amuse by running 'configure --enable-cuda'")
               
            return
            
        if self.config and not self.config.cuda.is_enabled:
            self.found_cuda = True
            self.environment['CUDA_LIBDIRS'] = '-L/NOCUDACONFIGURED/lib' + ' -LNOCUDACONFIGURED/lib64'
            self.environment['CUDA_LIBS'] = '-lnocuda'
            self.environment['CUDART_LIBS'] = '-lnocudart'
            self.environment['CUDA_TK'] = '/NOCUDACONFIGURED'
            self.environment['CUDA_SDK'] = '/NOCUDACONFIGURED'
            return 

        for x in ['CUDA_TK', 'CUDA_SDK']:
            if not x in self.environment:
                all_found = False
                break
        
        if all_found:
            cuda_dir = self.environment['CUDA_TK']
            self.environment['CUDA_LIBDIRS'] = '-L'+cuda_dir+'/lib' +  ' -L'+cuda_dir+'/lib64'
            self.environment['CUDA_LIBS'] = '-lcudart'
            return
            
        dir = spawn.find_executable('nvcc')
        if dir is None:
            self.found_cuda = False
            self.environment_notset['CUDA_SDK'] = '<directory>'
            self.environment_notset['CUDA_TK'] = '<directory>'
            return
        cuda_dir = os.path.dirname(os.path.dirname(dir))
        self.environment['CUDA_LIBDIRS'] = '-L'+cuda_dir+'/lib' + ' -L'+cuda_dir+'/lib64'
        self.environment['CUDA_LIBS'] = '-lcudart'
        self.environment['CUDA_TK'] = cuda_dir
        if not 'CUDA_SDK' in self.environment:
            self.environment_notset['CUDA_SDK'] = '<directory>'
        
        self.found_cuda = True

    def set_mpi_variables(self):
        if self.config:
            self.environment['MPICXX'] = self.config.mpi.mpicxx
            self.environment['MPICC'] = self.config.mpi.mpicc
            self.environment['MPIF90'] = self.config.mpi.mpif95
            return

    
    def set_compiler_variables(self):
        if self.config and not hasattr(self.config.compilers, 'found_fftw'):
            raise DistutilsError("configuration is not up to date, please reconfigure amuse by running 'configure'")
            
        if self.config:
            self.environment['CXX'] = self.config.compilers.cxx
            self.environment['CC'] = self.config.compilers.cc
            self.environment['FC'] = self.config.compilers.fc
            self.environment['CFLAGS'] = self.config.compilers.cc_flags
            self.environment['CXXFLAGS'] = self.config.compilers.cxx_flags
            self.environment['FFLAGS'] = self.config.compilers.fc_flags
            
            if self.config.compilers.found_fftw == 'yes':
                self.environment['FFTW_FLAGS'] = self.config.compilers.fftw_flags
                self.environment['FFTW_LIBS'] = self.config.compilers.fftw_libs
            
            
            if self.config.compilers.found_gsl == 'yes':
                self.environment['GSL_FLAGS'] = self.config.compilers.gsl_flags
                self.environment['GSL_LIBS'] = self.config.compilers.gsl_libs
                
            return

    def set_java_variables(self):
        if self.config and hasattr(self.config, 'java') and hasattr(self.config.java, 'is_enabled') and self.config.java.is_enabled:
            self.environment['JAVA'] = self.config.java.java
            self.environment['JAVAC'] = self.config.java.javac
            self.environment['JAR'] = self.config.java.jar
        else:
            self.environment['JAVA'] = ''
            self.environment['JAVAC'] = ''
            self.environment['JAR'] = ''
        return

    def set_openmp_flags(self):
        if self.config and hasattr(self.config, 'openmp'):
            self.environment['OPENMP_FCFLAGS'] = self.config.openmp.fcflags
            self.environment['OPENMP_CFLAGS'] = self.config.openmp.cflags
        else:
            self.environment['OPENMP_FCFLAGS'] = ''
            self.environment['OPENMP_CFLAGS'] = ''
            
    def set_libdir_variables(self):
        for varname in ('SAPPORO_LIBDIRS', 'GRAPE6_LIBDIRS'):
            if varname in self.environment:
                continue
                
            if varname in os.environ:
                self.environment[varname] = os.environ[varname]
            else:
                self.environment_notset[varname] ='-L<directory>'
        
        if 'SAPPORO_LIBDIRS' in self.environment:
            self.environment['SAPPORO_LIBS'] = '-L{0} -lsapporo'.format(
                self.environment['SAPPORO_LIBDIRS']
            )
        else:
            if self.config and hasattr(self.config.cuda, 'sapporo_version'):
                if self.config.cuda.sapporo_version == '2':
                    self.environment['SAPPORO_LIBS'] = '-L{0}/lib/sapporo-2 -lsapporo {1}'.format(
                        os.path.abspath(os.getcwd()),
                        self.config.openmp.cflags
                    )
                else:
                    self.environment['SAPPORO_LIBS'] = '-L{0}/lib/sapporo_light -lsapporo'.format(
                        os.path.abspath(os.getcwd())
                    )
            else:
                self.environment['SAPPORO_LIBS'] = '-L{0}/lib/sapporo_light -lsapporo'.format(
                    os.path.abspath(os.getcwd())
                )
            self.environment['BOOSTLIBS'] = ''
     
    def set_libs_variables(self):
        for varname, libname in []:
            if varname in self.environment:
                continue
                
            if varname in os.environ:
                self.environment[varname] = os.environ[varname]
            else:
                self.environment_notset[varname] ='-L<directory> -l{0}'.format(libname)

    def copy_config_to_build_dir(self):
        configpath=os.path.abspath(os.getcwd())
        topath=os.path.join(self.build_lib, "amuse")
        self.copy_file(os.path.join(configpath,"config.mk"), topath) 
     
    def copy_build_prereq_to_build_dir(self):
        if not os.path.exists(self.build_temp):
            self.mkpath(self.build_temp)

        if supportrc["framework_install"]:
            configpath=os.path.abspath(os.getcwd())
            self.copy_file(os.path.join(configpath,"config.mk"), self.build_temp) 
            self.copy_file(os.path.join(configpath,"build.py"), self.build_temp) 
        #~ self.copy_tree(os.path.join(configpath,"support"), os.path.join(self.build_temp,"support") )
        #~ self.copy_tree(os.path.join(configpath,"src"), os.path.join(self.build_temp,"src") )
        path=os.path.join(self.build_temp,"src")
        if not os.path.exists(path) and not os.path.islink(path):
            os.symlink(os.path.relpath(self.build_lib,self.build_temp), path)
        
    def copy_codes_to_build_dir(self):

        for dir in self.makefile_paths(self.codes_src_dir):
            reldir = os.path.relpath(dir, self.codes_src_dir)
            self.copy_tree(
                dir, 
                os.path.join(self.codes_dir, reldir)
            )
        
    def copy_lib_to_build_dir(self):
        for dir in self.makefile_paths(self.lib_src_dir):
            reldir = os.path.relpath(dir, self.lib_src_dir)
            self.copy_tree(
                dir, 
                os.path.join(self.lib_dir, reldir)
            )


    def copy_worker_codes_to_build_dir(self):
        if sys.platform == 'win32':
            worker_code_re = re.compile(r'([a-zA-Z0-9]+_)?worker(_[a-zA-Z0-9]+)?(.exe)?')
        else:
            worker_code_re = re.compile(r'([a-zA-Z0-9]+_)?worker(_[a-zA-Z0-9]+)?')
            worker_so_re = re.compile(r'([a-zA-Z0-9]+_)?cython(_[a-zA-Z0-9]+)?.so')
            
        
        lib_binbuilddir = os.path.join(self.build_lib, supportrc["package_name"], '_workers')
        if not os.path.exists(lib_binbuilddir):
            self.mkpath(lib_binbuilddir)
            
        for srcdir in self.makefile_paths(self.codes_src_dir):
            reldir = os.path.relpath(srcdir, self.codes_src_dir)
            temp_builddir = os.path.join(self.codes_dir, reldir)
            
            self.announce("will copy worker: {0}".format(srcdir), level = log.INFO)
            lib_builddir = os.path.join(self.build_lib, os.path.relpath(srcdir, os.path.join(self.amuse_src_dir, '..')))
            
            
            shortname = reldir.lower()
            self.announce(shortname, level = log.INFO)
            
            for name in os.listdir(temp_builddir):
                path = os.path.join(temp_builddir, name)
                stat = os.stat(path)
                
                if os.path.isfile(path):
                    if worker_so_re.match(name):
                        topath = os.path.join(lib_builddir, name)
                        self.copy_file(path, topath)
                        continue


                #self.announce("will copy worker: {0}".format(name), level = log.INFO)
                if os.path.isfile(path) and os.access(path, os.X_OK):
                    if worker_code_re.match(name):
                        topath = os.path.join(lib_binbuilddir, name)
                        self.copy_file(path, topath)
                    elif not name.endswith('.py'):
                        self.announce("will not copy executable: {0}, it does not match the worker pattern".format(name), level = log.WARN)
            
            # also copy file or dir named data
            path=os.path.join(temp_builddir,'data')
            topath = os.path.join(lib_builddir, 'data')
            if os.path.isfile(path):
                self.copy_file(path, topath)                
            if os.path.isdir(path):
                self.copy_tree(path, topath)                
            
                                    
    def subdirs_in_path(self,path):
        if not os.path.exists(path):
            return

        names = sorted(os.listdir(path))
        for name in names:
            if name.startswith('.'):
                continue
                
            path_ = os.path.join(path, name)
            if os.path.isdir(path_):
                yield path_

    def makefile_paths(self,path):
        for x in self.subdirs_in_path(path):
            for name in ('makefile', 'Makefile'):
                makefile_path = os.path.join(x, name)
                if os.path.exists(makefile_path):
                    yield x
                    break

    def update_environment_from_cfgfile(self):
        if os.path.exists('amuse.cfg'):
            config = configparser.ConfigParser()
            config.read(['amuse.cfg'])
            for name, value in config.items('environment'):
                if isinstance(value, str) and value:
                    varname = name.upper()
                    self.environment[varname] = value
                    if varname in self.environment_notset:
                        del self.environment_notset[varname]
                        
    def save_cfgfile_if_not_exists(self):
        if not os.path.exists('amuse.cfg'):
            config = configparser.RawConfigParser()
            config.add_section('environment')
            for name, value in self.environment.items():
                config.set('environment', name, value)
                
            for name, value in self.environment_notset.items():
                config.set('environment', name, '')
            
            with open('amuse.cfg', 'w') as f:
                config.write(f)

    
    
        
    def get_special_targets(self, name, directory, environment):
        process = Popen(['make','-qp', '-C', directory], env = environment, stdout = PIPE, stderr = PIPE)
        stdoutstring, stderrstring = process.communicate()
        if sys.hexversion > 0x03000000:
            stdoutstring = str(stdoutstring, 'utf-8')
        lines = stdoutstring.splitlines()
        result = []
        for line in lines:
            if line.startswith('muse_worker_gpu:'):
                result.append(('muse_worker_gpu', 'GPU',))
            elif line.startswith('muse_worker_grape:'):
                result.append(('muse_worker_grape', 'GRAPE6',))
            elif line.startswith('muse_worker_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    targetname = line[len('muse_worker_'):index_of_the_colon]
                    if '%' not in targetname: result.append((line[:index_of_the_colon], targetname,))
            elif line.startswith('worker_code_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    targetname = line[len('worker_code_'):index_of_the_colon]
                    if '%' not in targetname: result.append((line[:index_of_the_colon], targetname,))
            elif line.startswith(name + '_worker_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    targetname = line[len(name + '_worker_'):index_of_the_colon]
                    if '%' not in targetname: result.append((line[:index_of_the_colon], targetname,))
                    
        return result
    
    def call(self, arguments, buildlogfile = None, **keyword_arguments):
        stringio = StringIO()
         
        self.announce(' '.join(arguments), log.DEBUG)
        
        process = Popen(
            arguments, 
            stdout = PIPE,
            stderr = STDOUT,
            **keyword_arguments
        )
        
        while True:
            line = process.stdout.readline()
            if len(line) == 0:
                break
            
            if not buildlogfile is None:
                buildlogfile.write(line)
            self.announce(line[:-1], log.DEBUG)
            if sys.hexversion > 0x03000000:
                stringio.write(str(line, 'utf-8'))
            else:
                stringio.write(line)
            
        result = process.wait()
        content = stringio.getvalue()
        
        stringio.close()
        return result, content
        
    def build_environment(self):
        environment=self.environment.copy()
        environment.update(os.environ)
        path=os.path.join(environment["MUSE_PACKAGE_DIR"],"src")
        if environment["MUSE_PACKAGE_DIR"]!=environment["AMUSE_DIR"]:
            path=path+":"+os.path.join(environment["AMUSE_DIR"],"src")
        path=path+':'+environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"]=path
        return environment


    def do_clean(self):
        environment = self.build_environment()        

        for x in self.makefile_paths(self.lib_dir):
            self.announce("cleaning libary " + x)
            self.call(['make','-C', x, 'clean'], env=environment)
           
        for x in self.makefile_paths(self.codes_dir):
            if os.path.exists(x):
                self.announce("cleaning " + x)
                self.call(['make','-C', x, 'clean'], env=environment)

    def do_distclean(self):
        environment = self.build_environment()        

        for x in self.makefile_paths(self.lib_dir):
            self.announce("cleaning libary:" + x)
            self.call(['make','-C', x, 'distclean'], env=environment)
            
        for x in self.makefile_paths(self.codes_dir):
            self.announce("cleaning community code:" + x)
            self.call(['make','-C', x, 'distclean'], env=environment)


    
class SplitOutput(object) :
    def __init__(self, file1, file2) :
        self.file1 = file1
        self.file2 = file2

    def __del__(self) :
        self.close()

    def close(self):
        self.file1.close()
        self.file2.close()

    def write(self, text) :
        self.file1.write(text)
        self.file2.write(text)

    def flush(self) :
        self.file1.flush()
        self.file2.flush()
    
        
class BuildCodes(CodeCommand):

    description = "build interfaces to codes"

    user_options = list(CodeCommand.user_options)
    user_options.append( ('clean=', 'c', "clean code",), )

    def initialize_options(self):
        CodeCommand.initialize_options(self)
        self.clean = 'no'
        
    def finalize_options (self):
        CodeCommand.finalize_options(self)
        self.must_clean = self.clean == 'yes'
        self.must_dist_clean = self.clean == 'dist'
    
    def run_make_on_directory(self, codename, directory, target, environment):
        buildlog = os.path.abspath("build.log")
        
        with open(buildlog, "a") as output:
            output.write('*'*100)
            output.write('\n')
            output.write('Building code: {0}, target: {1}, in directory: {2}\n'.format(codename, target, directory))
            output.write('*'*100)
            output.write('\n')
            output.flush()
        
        with open(buildlog, "ab") as output:
            result, resultcontent = self.call(
                ['make','-C', directory, target], 
                output,
                env = environment
            )
        
        with open(buildlog, "a") as output:
            output.write('*'*100)
            output.write('\n')
            
        return result, resultcontent
    
    def is_download_needed(self, string):
        for line in string.splitlines():
            if 'DOWNLOAD_CODES' in line:
                return True
        return False
        
    def is_cuda_needed(self, string):
        for line in string.splitlines():
            if 'CUDA_TK variable is not set' in line:
                return True
            if 'CUDA_SDK variable is not set' in line:
                return True
        return False
        
    def are_python_imports_needed(self, string):
        for line in string.splitlines():
            if 'Python imports not available' in line:
                return True
        return False
    
    def run (self):
        if self.must_clean:
            self.do_clean()
        if self.must_dist_clean:
            self.do_distclean()

        not_build = list()
        is_download_needed = list()
        is_cuda_needed = list()
        not_build_special = {}
        are_python_imports_needed = list()
        build = list()
        lib_build = list()
        lib_not_build = list()
        environment = self.build_environment()
        
        buildlog = 'build.log'
        
        self.announce("building libraries and community codes", level = log.INFO)
        self.announce("build, for logging, see '{0}'".format(buildlog), level = log.INFO)
                
        with open(buildlog, "w") as output:
            output.write('*'*100)
            output.write('\n')
            output.write('Building libraries and codes\n')
            output.write('*'*100)
            output.write('\n')
        
        if not self.lib_dir == self.lib_src_dir:
            self.copy_build_prereq_to_build_dir()
            self.copy_lib_to_build_dir()
            if sys.hexversion > 0x03000000:
                run_2to3_on_build_dirs(self.makefile_paths(self.lib_src_dir), self.lib_dir,self.lib_src_dir)
                          
        for x in self.makefile_paths(self.lib_dir):
            
            shortname = x[len(self.lib_dir) + 1:] + '-library'
            starttime = datetime.datetime.now()
            self.announce("[{1:%H:%M:%S}] building {0}".format(shortname, starttime), level =  log.INFO)
            returncode, outputlog = self.run_make_on_directory(shortname, x, 'all', environment)
            
            endtime = datetime.datetime.now()
            if returncode == 2:
                self.announce("[{2:%H:%M:%S}] building {0}, failed, see {1!r} for error log".format(shortname, buildlog, endtime), level =  log.DEBUG)
                if self.is_download_needed(outputlog):
                    is_download_needed.append(x[len(self.lib_dir) + 1:])
                elif self.is_cuda_needed(outputlog):
                    is_cuda_needed.append(x[len(self.lib_dir) + 1:])
                else:
                    lib_not_build.append(shortname)
            else:
                self.announce("[{1:%H:%M:%S}] building {0}, succeeded".format(shortname, endtime), level =  log.DEBUG)
                lib_build.append(shortname)
            
        if not self.codes_dir == self.codes_src_dir:
            self.copy_codes_to_build_dir()
            if sys.hexversion > 0x03000000:
                run_2to3_on_build_dirs(self.makefile_paths(self.codes_src_dir), self.codes_dir,self.codes_src_dir)
        
        #environment.update(self.environment)
        makefile_paths = list(self.makefile_paths(self.codes_dir))

        build_to_special_targets = {}
        
        for x in makefile_paths:
            shortname = x[len(self.codes_dir) + 1:].lower()
            starttime = datetime.datetime.now()
            # For binary builds we do not want
            # to distribute mesa, it will make the
            # download size from about 100mb size 
            # to > 1Gb size.
            #
            # Could we remove some of the data files from mesa?
            #
            if not self.inplace and shortname == 'mesa':
                self.announce("[{1:%H:%M:%S}] skipping {0}".format(shortname, starttime), level =  log.INFO)
                continue 
                
            self.announce("[{1:%H:%M:%S}] building {0}".format(shortname, starttime), level =  log.INFO)
            returncode, outputlog = self.run_make_on_directory(shortname, x, 'all', environment)
            endtime = datetime.datetime.now()
            if returncode > 0:
                self.announce("[{2:%H:%M:%S}] building {0}, failed, see {1!r} for error log".format(shortname, buildlog, endtime), level =  log.DEBUG)
                if self.is_download_needed(outputlog):
                    is_download_needed.append(shortname)
                elif self.is_cuda_needed(outputlog):
                    is_cuda_needed.append(shortname)
                elif self.are_python_imports_needed(outputlog):
                    are_python_imports_needed.append(shortname)
                else:
                    not_build.append(shortname)
                    
                if self.is_mpi_enabled():
                    continue
            else:
                build.append(shortname)
                is_built = True
                self.announce("[{1:%H:%M:%S}] building {0}, succeeded".format(shortname, endtime), level =  log.DEBUG)
            
            if not self.variant:
                continue
                
            special_targets = self.get_special_targets(shortname, x, environment)
            for target,target_name in special_targets:
                starttime = datetime.datetime.now()
                self.announce("[{2:%H:%M:%S}] building {0} - {1}".format(shortname, target_name, starttime), level =  log.DEBUG)
                returncode, outputlog = self.run_make_on_directory(shortname, x, target, environment)
                endtime = datetime.datetime.now()
                if returncode > 0:
                    specials_list = not_build_special.setdefault(shortname,[])
                    specials_list.append(target_name)
                    self.announce("[{3:%H:%M:%S}] building {0} - {1}, failed, see {2!r} for error log".format(shortname, target_name, buildlog,endtime), level =  log.DEBUG)
                else:
                    build_to_special_targets.setdefault(shortname, list()).append(target_name)
                    self.announce("[{2:%H:%M:%S}] building {0} - {1}, succeeded".format(shortname, target_name, endtime), level =  log.DEBUG)
                
        
        if not self.codes_dir == self.codes_src_dir:
            if supportrc["framework_install"]:
                self.copy_config_to_build_dir()
            self.copy_worker_codes_to_build_dir()
            
        with open(buildlog, "a") as output:
            output.write('*'*80)
            output.write('\n')
            output.write('Building finished\n')
            output.write('*'*80)
            output.write('\n')
            
        self.announce("Environment variables")
        self.announce("="*80)
        sorted_keys = sorted(self.environment.keys())
        for x in sorted_keys:
            self.announce("%s\t%s" % (x , self.environment[x] ))
        
        if not self.is_mpi_enabled():
            print(build_to_special_targets)
            all_build = set(build)
            not_build_copy = []
            for x in not_build:
                if x in build_to_special_targets:
                    if not x in all_build:
                        build.append(x)
                        all_build.add(x)
                else:
                    not_build_copy.append(x)
            not_build = not_build_copy
                
        
        if not_build or not_build_special or is_download_needed or is_cuda_needed or are_python_imports_needed:
            if not_build:
                level = log.WARN
            else:
                level = log.INFO
            if not_build:
                self.announce("Community codes not built (because of errors):",  level = level)
                self.announce("="*80,  level = level)
                for x in not_build:
                    self.announce(' * {0}'.format(x), level =  level)
            if not_build_special:
                self.announce("Optional builds failed, need special libraries:",  level = level)
                for x in sorted(not_build_special.keys()):
                    self.announce(' * {0} - {1}'.format(x, ', '.join(not_build_special[x])), level = level)
            if is_cuda_needed:
                self.announce("Optional builds failed, need CUDA/GPU libraries:",  level = level)
                for x in is_cuda_needed:
                    self.announce(' * {0}'.format(x), level = level)
            if are_python_imports_needed:
                self.announce("Optional builds failed, need additional python packages:",  level = level)
                for x in are_python_imports_needed:
                    self.announce(' * {0}'.format(x), level = level)
            if is_download_needed:
                self.announce("Optional builds failed, need separate download",  level = level)
                for x in is_download_needed:
                    self.announce(' * {0} , make {0}.code DOWNLOAD_CODES=1'.format(x), level = level)

            self.announce("="*80,  level = level)
        
        if build:
            level = log.INFO
            self.announce("Community codes built",  level = level)
            self.announce("="*80,  level = level)
            for x in build:
                if x in build_to_special_targets:
                    y = build_to_special_targets[x]
                    self.announce('* {0} ({1})'.format(x,','.join(y)),  level = level)
                else:
                    self.announce('* {0}'.format(x),  level = level)
            self.announce("="*80,  level = level)
        
        level = log.INFO
        self.announce(
            "{0} out of {1} codes built, {2} out of {3} libraries built".format(
                len(build), 
                len(build) + len(not_build), 
                len(lib_build), 
                len(lib_build) + len(lib_not_build)
            ),  
            level = level
        )
        
        if self.config and (not hasattr(self.config, 'java') or not hasattr(self.config.java, 'is_enabled')):
            self.announce(
                "Your configuration is out of date, please rerun configure",
                level = level
            )
        
 
class BuildLibraries(BuildCodes):

    description = "build just the supporting libraries"

    def subdirs_in_path(self,path):
        # bit hackish way to filter out non lib stuff
        if path not in [self.lib_dir, self.lib_src_dir]:
            return 
            
        if not os.path.exists(path):
            return

        names = sorted(os.listdir(path))
        for name in names:
            if name.startswith('.'):
                continue
                
            path_ = os.path.join(path, name)
            if os.path.isdir(path_):
                yield path_

class ConfigureCodes(CodeCommand):

    description = "run configure for amuse"

    def run (self):

        if os.path.exists('config.mk') or self.config:
            self.announce("Already configured, not running configure", level = 2)
            return
        environment = self.build_environment()
        self.announce("Running configure for AMUSE", level = 2)
        result,content=self.call(['./configure'], env=environment, shell=True)
        if not os.path.exists('config.mk'):
            self.announce("config.mk not generated; output of configure:", level=2)
            self.announce(content, level=2)
            raise Exception("configure failed")
        with open("config.mk") as infile:
            self.announce("configure generated config.mk", level=2)
            self.announce("="*80, level=2)
            for line in infile:
                self.announce(line[:-1], level=2)
            self.announce("="*80, level=2)

        
class CleanCodes(CodeCommand):

    description = "clean build products in codes"

    def run (self):
                  
        self.announce("Cleaning libraries and community codes", level = 2)
        self.do_clean()
        
 
class DistCleanCodes(CodeCommand):

    description = "clean for distribution"

    def run (self):      
      
        self.announce("Cleaning for distribution, libraries and community codes", level = 2)
        self.do_distclean()
        
class BuildOneCode(BuildCodes):  
    description = "build one code"

    user_options = list(BuildCodes.user_options)
    user_options.append( ('code-name=', 'n', "name of the code",), )

    def initialize_options(self):
        BuildCodes.initialize_options(self)
        self.code_name = None
        
    def finalize_options (self):
        BuildCodes.finalize_options(self)
        if self.code_name is None:
            raise Exception("no code was specified")

    def subdirs_in_path(self,path):
        if not os.path.exists(path):
            return
            
        names = os.listdir(path)
        for name in names:
            if name.startswith('.'):
                continue
            if not name.lower().startswith(self.code_name.lower()):
                continue
            path_ = os.path.join(path, name)
            if os.path.isdir(path_):
                yield path_

    def run(self):
        if not self.inplace:
            self.run_command("build_py")
        
        BuildCodes.run(self)

class Clean(clean):

    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

class Install(install):

    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

        install.run(self)

def setup_commands():
    mapping_from_command_name_to_command_class = {
        'build_codes': BuildCodes,
        'build_code': BuildOneCode,
        'clean_codes': CleanCodes,
        'dist_clean': DistCleanCodes,
        'clean_python': clean,
        'clean': Clean,
        'install': install,
        'build_libraries': BuildLibraries,
        'install_libraries': InstallLibraries
    }
    
    if sys.hexversion > 0x03000000:
        mapping_from_command_name_to_command_class['build_py'] = build_py_2to3
    
    build.sub_commands.append(('build_codes', None))
    Clean.sub_commands.append(('clean_codes', None))
    Clean.sub_commands.append(('clean_python', None))
    Install.sub_commands.append(('install_libraries', None))
    
    if supportrc["framework_install"]:
        mapping_from_command_name_to_command_class.update(
            {
                'configure_codes': ConfigureCodes,
                'generate_install_ini': GenerateInstallIni,
                'build_latex': build_latex,
                'generate_main': generate_main,
                'tests': run_tests,
            }
        )
        build.sub_commands.insert(0, ('configure_codes', None))
        Install.sub_commands.insert(0, ('generate_install_ini', None))



    
    return mapping_from_command_name_to_command_class
