__revision__ = "$Id:$"

import sys, os, re, subprocess
import ConfigParser
import os.path

from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.util import convert_path
from distutils import log
from distutils import spawn
from subprocess import call, Popen, PIPE
from numpy.distutils import fcompiler
# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')

class LegacyCommand(Command):
    user_options = [
        ('legacy-dir=', 'd', "directory containing legacy codes"),
        ('lib-dir=', 'l', "directory containing libraries to build"),
    ]

    boolean_options = ['force']

    def initialize_options (self):
        self.legacy_dir = None
        self.lib_dir = None
        self.amuse_src_dir =  os.path.join('src','amuse')
        self.environment = {}
        self.environment_notset = {}
        self.found_cuda = False
        self.found_sapporo = False
        
    def finalize_options (self):
        if self.legacy_dir is None:
            self.legacy_dir = os.path.join(self.amuse_src_dir,'legacy')
        
        if self.lib_dir is None:
            self.lib_dir = 'lib'
        
        #self.update_environment_from_cfgfile()
        self.set_fortran_variables()
        
        self.environment['F90'] = self.environment['FORTRAN']
        self.environment['FC'] = self.environment['FORTRAN']
        
        self.set_cuda_variables()
        self.set_libdir_variables()
        self.set_libs_variables()
        self.save_cfgfile_if_not_exists()
    
    
    def set_fortran_variables(self):
        if 'FORTRAN' in self.environment:
            return
            
        if 'FORTRAN' in os.environ:
            self.environment['FORTRAN'] = os.environ['FORTRAN']
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
    
    
    def set_cuda_variables(self):
        all_found = True
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

    
    def set_libdir_variables(self):
        for varname in ('SAPPORO_LIBDIRS', 'GRAPE6_LIBDIRS'):
            if varname in self.environment:
                continue
                
            if varname in os.environ:
                self.environment[varname] = os.environ[varname]
            else:
                self.environment_notset[varname] ='-L<directory>'
     
    def set_libs_variables(self):
        for varname, libname in (('PGLIBS','pg5'),):
            if varname in self.environment:
                continue
                
            if varname in os.environ:
                self.environment[varname] = os.environ[varname]
            else:
                self.environment_notset[varname] ='-L<directory> -l{0}'.format(libname)
     
    
            
    def subdirs_in_legacy_dir(self):
        names = os.listdir(self.legacy_dir)
        for name in names:
            if name.startswith('.'):
                continue
            path = os.path.join(self.legacy_dir, name)
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
                    
        for x in self.subdirs_in_legacy_dir():
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

    
    def get_special_targets(self, directory, environment):
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
                    name = line[len('muse_worker_'):index_of_the_colon]
                    result.append((line[:index_of_the_colon], name,))
            elif line.startswith('worker_code_'):
                index_of_the_colon = line.index(':')
                if(index_of_the_colon > 0):
                    name = line[len('worker_code_'):index_of_the_colon]
                    result.append((line[:index_of_the_colon], name,))
        return result
        
class BuildLegacy(LegacyCommand):

    description = "build interfaces to legacy codes"
    
    def run_make_on_directory(self, codename, directory, target, environment):
        buildlog = "{0}-{1}-build.log".format(codename, target)
        with open(buildlog, "w") as output:
            process = Popen(['make','-C', directory, target], env = environment, stdout = output, stderr = output)
            return process.wait(), buildlog
    
    def run (self):
        not_build = list()
        not_build_special = list()
        build = list()
        environment = self.environment
        environment.update(os.environ)
        
        
        for x in self.makefile_libpaths():
            
            shortname = x[len(self.lib_dir) + 1:] + '-library'
            self.announce("building {0}".format(shortname), level =  log.INFO)
            returncode, buildlog = self.run_make_on_directory(shortname, x, 'all', environment)
            
            if returncode == 2:
                self.announce("building {0}, failed, see {1} for error log".format(shortname, buildlog), level =  log.INFO)
                not_build.append(shortname)
            else:
                self.announce("building {0}, succeeded".format(shortname), level =  log.INFO)
                build.append(shortname)
            
        #environment.update(self.environment)
        makefile_paths = list(self.makefile_paths())
        
        for x in makefile_paths:
            shortname = x[len(self.legacy_dir) + 1:]
            self.announce("building {0}".format(shortname), level =  log.INFO)
            returncode, buildlog = self.run_make_on_directory(shortname, x, 'all', environment)
            if returncode > 0:
                not_build.append(shortname)
                self.announce("building {0}, failed, see {1} for error log".format(shortname, buildlog), level =  log.INFO)
            else:
                build.append(shortname)
                self.announce("building {0}, succeeded".format(shortname), level =  log.INFO)
            
            special_targets = self.get_special_targets(x, environment)
            for target,target_name in special_targets:
                self.announce("building " + x + " version: " + target_name, level =  log.INFO)
                returncode, buildlog = self.run_make_on_directory(shortname, x, target, environment)
                if returncode > 0:
                    not_build_special.append(shortname + " - " + target_name)
                    self.announce("building {0} - {1}, failed, see {2!r} for error log".format(shortname, target_name, buildlog), level =  log.INFO)
                else:
                    build.append(shortname + " - " + target_name)
                    self.announce("building {0} - {1}, succeeded".format(shortname, target_name), level =  log.INFO)
                
        
        print
        #print "Environment variables"
        #print "====================="
        #sorted_keys = sorted(self.environment.keys())
        #for x in sorted_keys:
        #    print "%s\t%s" % (x , self.environment[x] )
        #print
        print "Environment variables not set"
        print "============================="
        sorted_keys = sorted(self.environment_notset.keys())
        for x in sorted_keys:
            print "%s\t%s" % (x, self.environment_notset[x] )
        print
        print
        if not_build or not_build_special:
            print
            print
            print "Legacy codes not built (because of errors):"
            print "==========================================="
            for x in not_build:
                print '*', x 
            for x in not_build_special:
                print '*', x, '** optional, needs special libraries or hardware to compile **'
        if build:
            print
            print
            print "Legacy codes built"
            print "=================="
            for x in build:
                print '*', x
        
        
 
class CleanLegacy(LegacyCommand):

    description = "clean build products in legacy codes"

    def run (self):
        for x in self.makefile_libpaths():
            self.announce("cleaning libary " + x)
            call(['make','-C', x, 'clean'])
            
        for x in self.makefile_paths():
            self.announce("cleaning " + x)
            call(['make','-C', x, 'clean'])
        
class BuildOneLegacyCode(LegacyCommand):  
    description = "build one code"
    user_options = list(LegacyCommand.user_options)
    user_options.append( ('code-name=', 'n', "name of the legacy code",), )
    
    
    def initialize_options(self):
        LegacyCommand.initialize_options(self)
        self.code_name = None
        
    
    def finalize_options (self):
        LegacyCommand.finalize_options(self)
        if self.code_name is None:
            raise Exception("no legacy code was specified")
    
    
    def subdirs_in_legacy_dir(self):
        names = os.listdir(self.legacy_dir)
        for name in names:
            if name.startswith('.'):
                continue
            if not name.lower().startswith(self.code_name.lower()):
                continue
            path = os.path.join(self.legacy_dir, name)
            if os.path.isdir(path):
                yield path
                
                
    def run (self):
        environment = self.environment
        environment.update(os.environ)
        
        results = []
        for x in self.makefile_paths():
            self.announce("cleaning " + x)
            call(['make','-C', x, 'clean'])
            returncode = call(['make','-C', x, 'all'], env = environment)
            results.append(('default',returncode,))
            
            special_targets = self.get_special_targets(x, environment)
            for target,target_name in special_targets:
                self.announce("building " + x + " version: " + target_name)
                returncode = call(['make','-C', x, target], env = environment)
                results.append((target,returncode,))
        
        for name, returncode in results:
            print name, "...", "failed" if returncode == 2 else "succeeded"
            

                
            
            
            
            
    
            
           
            
            
            
            
    
            
