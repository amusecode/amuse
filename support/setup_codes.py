__revision__ = "$Id:$"

import sys, os, re, subprocess
import ConfigParser
import os.path
import datetime

from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.util import convert_path
from distutils import log
from distutils import spawn
from distutils.errors import DistutilsError
from subprocess import call, Popen, PIPE, STDOUT
from numpy.distutils import fcompiler
import StringIO
# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')

try:
    from . import config
    is_configured = hasattr(config, 'compilers')
except ImportError:
    is_configured = False
    
    
class CodeCommand(Command):
    user_options = [
        ('build-lib=', 'b',
         "directory for compiled extension modules"),
        ('build-temp=', 't',
         "directory for temporary files (build by-products)"),
        ('inplace', 'i',
         "ignore build-lib and put compiled extensions into the source " +
         "directory alongside your pure Python modules"),
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
        ('code-dir=', 'd', "directory containing codes"),
        ('lib-dir=', 'l', "directory containing libraries to build"),
    ]

    boolean_options = ['force', 'inplace', 'debug', 'variant']

    def initialize_options (self):
        self.codes_dir = None
        self.lib_dir = None
        self.amuse_src_dir =  os.path.join('src','amuse')
        self.environment = {}
        self.environment_notset = {}
        self.found_cuda = False
        self.found_sapporo = False
        self.variant = False
        self.inplace = 1
        
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
        
        if self.codes_dir is None:
            if self.inplace:
                self.codes_dir = os.path.join(self.amuse_src_dir,'community')
            else:
                self.codes_dir = os.path.join(self.build_lib, 'amuse', 'community')
                
        if self.lib_dir is None:
            self.lib_dir = 'lib'
        
        
        self.set_cuda_variables()
        self.set_mpi_variables()
        self.set_compiler_variables()
        
        
        self.set_fortran_variables()
        
        self.environment['F90'] = self.environment['FORTRAN']
        self.environment['FC'] = self.environment['FORTRAN']
        self.set_java_variables()
        self.set_libdir_variables()
        self.set_libs_variables()
        self.save_cfgfile_if_not_exists()
    
    
    def set_fortran_variables(self):
        if 'FORTRAN' in self.environment:
            return
            
        if 'FORTRAN' in os.environ:
            self.environment['FORTRAN'] = os.environ['FORTRAN']
            return
            
        if is_configured:
            self.environment['FORTRAN'] = config.compilers.fc
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
            
        compiler = fcompiler.new_fcompiler(requiref90=True)
        fortran_executable = compiler.executables['compiler_f90'][0]
        self.environment['FORTRAN'] = fortran_executable
    
    
    
    def is_mpi_enabled(self):
        if is_configured and hasattr(config.mpi, 'is_enabled'):
            return config.mpi.is_enabled
        else:
            return True
    
    def set_cuda_variables(self):
        all_found = True
        if is_configured and config.cuda.is_enabled:
            self.found_cuda = True
            self.environment['CUDA_LIBDIRS'] = '-L'+config.cuda.toolkit_path+'/lib' + ' -L'+config.cuda.toolkit_path+'/lib64'
            self.environment['CUDA_TK'] = config.cuda.toolkit_path
            self.environment['CUDA_SDK'] = config.cuda.sdk_path
            if hasattr(config.cuda, 'cuda_libs'):
                self.environment['CUDA_LIBS'] = config.cuda.cuda_libs
            else:
                raise DistutilsError("configuration is not up to date for cuda, please reconfigure amuse by running 'configure --enable-cuda'")
               
            return
            
        if is_configured and not config.cuda.is_enabled:
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
        if is_configured:
            self.environment['MPICXX'] = config.mpi.mpicxx
            self.environment['MPICC'] = config.mpi.mpicc
            self.environment['MPIF90'] = config.mpi.mpif95
            return

    
    def set_compiler_variables(self):
        if is_configured and not hasattr(config.compilers, 'found_fftw'):
            raise DistutilsError("configuration is not up to date, please reconfigure amuse by running 'configure'")
            
        if is_configured:
            self.environment['CXX'] = config.compilers.cxx
            self.environment['CC'] = config.compilers.cc
            self.environment['FC'] = config.compilers.fc
            self.environment['CFLAGS'] = config.compilers.cc_flags
            self.environment['CXXFLAGS'] = config.compilers.cxx_flags
            self.environment['FFLAGS'] = config.compilers.fc_flags
            
            if config.compilers.found_fftw == 'yes':
                self.environment['FFTW_FLAGS'] = config.compilers.fftw_flags
                self.environment['FFTW_LIBS'] = config.compilers.fftw_libs
            
            
            if config.compilers.found_gsl == 'yes':
                self.environment['GSL_FLAGS'] = config.compilers.gsl_flags
                self.environment['GSL_LIBS'] = config.compilers.gsl_libs
                
            return

    def set_java_variables(self):
        if is_configured:
            self.environment['JNI_INCLUDES'] = config.java.jni_includes
            self.environment['JDK'] = config.java.jdk
        return

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
            if is_configured and hasattr(config.cuda, 'sapporo_version'):
                if config.cuda.sapporo_version == '2':
                    self.environment['SAPPORO_LIBS'] = '-L{0}/lib/sapporo-2 -lsapporo'.format(
                        os.path.abspath(os.getcwd())
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
     
    
            
    def subdirs_in_codes_dir(self):
        names = os.listdir(self.codes_dir)
        for name in names:
            if name.startswith('.'):
                continue
            path = os.path.join(self.codes_dir, name)
            if os.path.isdir(path):
                yield path
                
    def subdirs_in_lib_dir(self):
        names = os.listdir(self.lib_dir)
        for name in names:
            if name.startswith('.'):
                continue
            path = os.path.join(self.lib_dir, name)
            if os.path.isdir(path):
                yield path
            
    def makefile_libpaths(self):
        for x in self.subdirs_in_lib_dir():
            for name in ('makefile', 'Makefile'):
                makefile_path = os.path.join(x, name)
                if os.path.exists(makefile_path):
                    yield x
                    break
                    
    def makefile_paths(self):
        for x in self.subdirs_in_codes_dir():
            for name in ('makefile', 'Makefile'):
                makefile_path = os.path.join(x, name)
                if os.path.exists(makefile_path):
                    yield x
                    break

    def update_environment_from_cfgfile(self):
        if os.path.exists('amuse.cfg'):
            config = ConfigParser.ConfigParser()
            config.read(['amuse.cfg'])
            for name, value in config.items('environment'):
                if isinstance(value, str) and value:
                    varname = name.upper()
                    self.environment[varname] = value
                    if varname in self.environment_notset:
                        del self.environment_notset[varname]
                        
    def save_cfgfile_if_not_exists(self):
        if not os.path.exists('amuse.cfg'):
            config = ConfigParser.RawConfigParser()
            config.add_section('environment')
            for name, value in self.environment.iteritems():
                config.set('environment', name, value)
                
            for name, value in self.environment_notset.iteritems():
                config.set('environment', name, '')
            
            with open('amuse.cfg', 'wb') as f:
                config.write(f)

    
    def get_special_targets(self, name, directory, environment):
        process = Popen(['make','-qp', '-C', directory], env = environment, stdout = PIPE, stderr = PIPE)
        stdoutstring, stderrstring = process.communicate()
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
                    result.append((line[:index_of_the_colon], targetname,))
            elif line.startswith('worker_code_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    targetname = line[len('worker_code_'):index_of_the_colon]
                    result.append((line[:index_of_the_colon], targetname,))
            elif line.startswith(name + '_worker_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    targetname = line[len(name + '_worker_'):index_of_the_colon]
                    result.append((line[:index_of_the_colon], targetname,))
                    
        if not self.is_mpi_enabled():
            return [x for x in result if x[-1].endswith('sockets')]
        return result
    
    def call(self, arguments, buildlogfile = None, **keyword_arguments):
        stringio = StringIO.StringIO()
         
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
            stringio.write(line)
            
        result = process.wait()
        content = stringio.getvalue()
        
        stringio.close()
        return result, content
    
    
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
    
    
    def run_make_on_directory(self, codename, directory, target, environment):
        buildlog = os.path.abspath("build.log")
        
        with open(buildlog, "a") as output:
            output.write('*'*100)
            output.write('\n')
            output.write('Building code: {0}, target: {1}, in directory: {2}\n'.format(codename, target, directory))
            output.write('*'*100)
            output.write('\n')
            output.flush()
        
        with open(buildlog, "a") as output:
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
        
    def run (self):
        not_build = list()
        is_download_needed = list()
        is_cuda_needed = list()
        not_build_special = list()
        build = list()
        environment = self.environment
        environment.update(os.environ)
        
        buildlog = 'build.log'
        
        self.announce("building libraries and community codes", level = log.INFO)
        self.announce("build, for logging, see '{0}'".format(buildlog), level = log.INFO)
        
        with open(buildlog, "w") as output:
            output.write('*'*100)
            output.write('\n')
            output.write('Building libraries and codes\n')
            output.write('*'*100)
            output.write('\n')
        
        
        for x in self.makefile_libpaths():
            
            shortname = x[len(self.lib_dir) + 1:] + '-library'
            starttime = datetime.datetime.now()
            self.announce("[{1:%H:%M:%S}] building {0}".format(shortname, starttime), level =  log.INFO)
            returncode, outputlog = self.run_make_on_directory(shortname, x, 'all', environment)
            
            endtime = datetime.datetime.now()
            if returncode == 2:
                self.announce("[{2:%H:%M:%S}] building {0}, failed, see {1!r} for error log".format(shortname, buildlog, endtime), level =  log.DEBUG)
                if self.is_download_needed(outputlog):
                    is_download_needed.append(shortname)
                elif self.is_cuda_needed(outputlog):
                    is_cuda_needed.append(shortname)
                else:
                    not_build.append(shortname)
            else:
                self.announce("[{1:%H:%M:%S}] building {0}, succeeded".format(shortname, endtime), level =  log.DEBUG)
                build.append(shortname)
            
        #environment.update(self.environment)
        makefile_paths = list(self.makefile_paths())
        build_to_special_targets = {}
        
        for x in makefile_paths:
            shortname = x[len(self.codes_dir) + 1:].lower()
            starttime = datetime.datetime.now()
            self.announce("[{1:%H:%M:%S}] building {0}".format(shortname, starttime), level =  log.INFO)
            returncode, outputlog = self.run_make_on_directory(shortname, x, 'all', environment)
            endtime = datetime.datetime.now()
            if returncode > 0:
                self.announce("[{2:%H:%M:%S}] building {0}, failed, see {1!r} for error log".format(shortname, buildlog, endtime), level =  log.DEBUG)
                if self.is_download_needed(outputlog):
                    is_download_needed.append(shortname)
                elif self.is_cuda_needed(outputlog):
                    is_cuda_needed.append(shortname)
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
                    not_build_special.append(shortname + " - " + target_name)
                    self.announce("[{3:%H:%M:%S}] building {0} - {1}, failed, see {2!r} for error log".format(shortname, target_name, buildlog,endtime), level =  log.DEBUG)
                else:
                    build_to_special_targets.setdefault(shortname, list()).append(target_name)
                    self.announce("[{2:%H:%M:%S}] building {0} - {1}, succeeded".format(shortname, target_name, endtime), level =  log.DEBUG)
                
        
        with open(buildlog, "a") as output:
            output.write('*'*100)
            output.write('\n')
            output.write('Building finished\n')
            output.write('*'*100)
            output.write('\n')
            
        self.announce("Environment variables")
        self.announce("="*100)
        sorted_keys = sorted(self.environment.keys())
        for x in sorted_keys:
            self.announce("%s\t%s" % (x , self.environment[x] ))
        
        if not self.is_mpi_enabled():
            print build_to_special_targets
            all_build = set(build)
            for i, x in enumerate(list(not_build)):
                print i, x, x in build_to_special_targets
                if x in build_to_special_targets:
                    del not_build[i]
                    if not x in all_build:
                        build.append(x)
                        all_build.add(x)
                    
                
        
        if not_build or not_build_special or is_download_needed or is_cuda_needed:
            if not_build:
                level = log.WARN
            else:
                level = log.INFO
            
            self.announce("Community codes not built (because of errors):",  level = level)
            self.announce("="*100,  level = level)
            for x in not_build:
                self.announce(' * {0}'.format(x), level =  level)
            for x in not_build_special:
                self.announce(' * {0} ** optional, needs special libraries or hardware to compile **'.format(x), level = level)
            for x in is_cuda_needed:
                self.announce(' * {0} ** needs CUDA library, please configure with --enable-cuda **'.format(x), level = level)
            for x in is_download_needed:
                self.announce(' * {0} ** needs to be download, use make {0}.code DOWNLOAD_CODES=1 to download and build **'.format(x), level = level)

            self.announce("="*100,  level = level)
        
        if build:
            level = log.INFO
            self.announce("Community codes built",  level = level)
            self.announce("="*100,  level = level)
            for x in build:
                if x in build_to_special_targets:
                    y = build_to_special_targets[x]
                    self.announce('* {0} ({1})'.format(x,','.join(y)),  level = level)
                else:
                    self.announce('* {0}'.format(x),  level = level)
            self.announce("="*100,  level = level)
        
        
 
class CleanCodes(CodeCommand):

    description = "clean build products in codes"

    def run (self):
        buildlog = "make-clean.log"
            
        self.announce("Cleaning libraries and community codes", level = 2)
        with open(buildlog, "aw") as output:
            for x in self.makefile_libpaths():
                self.announce("cleaning libary " + x)
                call(['make','-C', x, 'clean'], stdout = output, stderr = output)
               
                
            for x in self.makefile_paths():
                self.announce("cleaning " + x)
                call(['make','-C', x, 'clean'], stdout = output, stderr = output)
 
class DistCleanCodes(CodeCommand):

    description = "clean for distribution"

    def run (self):
        buildlog = "make-clean.log"
            
        self.announce("Cleaning for distribution, libraries and community codes", level = 2)
        with open(buildlog, "aw") as output:
            for x in self.makefile_libpaths():
                self.announce("cleaning libary:" + x)
                call(['make','-C', x, 'distclean'], stdout = output, stderr = output)
                
            for x in self.makefile_paths():
                self.announce("cleaning community code:" + x)
                call(['make','-C', x, 'distclean'], stdout = output, stderr = output)
        
class BuildOneCode(CodeCommand):  
    description = "build one code"
    user_options = list(CodeCommand.user_options)
    user_options.append( ('code-name=', 'n', "name of the code",), )
    
    
    def initialize_options(self):
        CodeCommand.initialize_options(self)
        self.code_name = None
        
    
    def finalize_options (self):
        CodeCommand.finalize_options(self)
        if self.code_name is None:
            raise Exception("no code was specified")
    
    
    def subdirs_in_codes_dir(self):
        names = os.listdir(self.codes_dir)
        for name in names:
            if name.startswith('.'):
                continue
            if not name.lower().startswith(self.code_name.lower()):
                continue
            path = os.path.join(self.codes_dir, name)
            if os.path.isdir(path):
                yield path
        
    def subdirs_in_lib_dir(self):
        names = os.listdir(self.lib_dir)
        for name in names:
            if name.startswith('.'):
                continue
            if not name.lower().startswith(self.code_name.lower()):
                continue
            path = os.path.join(self.lib_dir, name)
            if os.path.isdir(path):
                yield path
                
    
    
    
    def run (self):
        environment = self.environment
        environment.update(os.environ)
        
        results = []
        
        for x in self.makefile_paths():
            shortname = x[len(self.codes_dir) + 1:].lower()
            
            self.announce("cleaning " + x)
            self.call(['make','-C', x, 'clean'], env=environment)
                        
            if self.is_mpi_enabled():
                returncode, _ = self.call(['make','-C', x, 'all'], env = environment)
                results.append(('default',returncode,))
            
            special_targets = self.get_special_targets(shortname, x, environment)
            for target,target_name in special_targets:
                self.announce("building " + x + " version: " + target_name)
                returncode, _ = self.call(['make','-C', x, target], env = environment)
                results.append((target,returncode,))
        
        for x in self.makefile_libpaths():
            shortname = x[len(self.codes_dir) + 1:].lower()
            
            self.announce("cleaning " + x)
            self.call(['make','-C', x, 'clean'], env=environment)

            returncode, _ = self.call(['make','-C', x, 'all'], env = environment)
            results.append(('default',returncode,))
            
            special_targets = self.get_special_targets(shortname, x, environment)
            for target,target_name in special_targets:
                self.announce("building " + x + " version: " + target_name)
                returncode, _ = self.call(['make','-C', x, target], env = environment)
                results.append((target,returncode,))
            
        
        for name, returncode in results:
            print name, "...", "failed" if returncode == 2 else "succeeded"
            

                
            
            
            
            
    
            
           
            
            
            
            
    
            
