#!/usr/bin/env python



import sys
import os.path
import os
import urllib
import subprocess
import shutil


IS_ON_OSX = sys.platform == 'darwin'

def late(function):
    class LateProperty(object):
        def __init__(self, initializer):
            self.initializer = initializer
        def __get__(self, instance, owner):
            if instance is None:
                return self
            value = self.initializer(instance)
            setattr(instance,self.initializer.__name__,value)
            return value
    return LateProperty(function)

    
        
class InstallPrerequisites(object):
    @late
    def prefix(self):
        return sys.prefix
            
    @late
    def static_prefix(self):
        return os.path.normpath(os.path.join(self.prefix, '..', 'static_libs'))
        
    @late
    def applications(self):
       return [
          #('openssl' , [], '0.9.8k' , 'openssl-', '.tar.gz', 'http://www.openssl.org/source/', self.openssl_build),
          (
            'numpy' ,                  #name to refer by
            [],                        #names of prerequisites (unused)
            '1.6.2' ,                  #version string
            'numpy-', '.tar.gz',       #pre- and postfix for filename
            'http://ignum.dl.sourceforge.net/sourceforge/numpy/', #download url, filename is appended
            self.numpy_build          #method to use for building
          ),
          (
            'nose', 
            [], 
            '1.1.2', 
            'nose-' , '.tar.gz', 
            'http://d.pypi.python.org/packages/source/n/nose/', 
            self.python_build
          ),
          (
            'distribute', 
            [], 
            '0.6.24', 
            'distribute-' , '.tar.gz', 
            'http://d.pypi.python.org/packages/source/d/distribute/', 
            self.python_build
          ),
          (
            'hdf' ,
            [],  
            '1.8.8',
            'hdf5-' , '.tar.gz' , 
            'http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.8/src/', 
            #'http://www.hdfgroup.org/ftp/HDF5/current/src/',
            self.hdf5_build
          ) ,
          (
            'h5py', 
            ['hdf'], 
            '2.1.0', 
            'h5py-' , '.tar.gz', 
            'http://h5py.googlecode.com/files/', self.h5py_build
          ) ,
          (
            'docutils', 
            [], 
            '0.9.1', 
            'docutils-','.tar.gz', 
            'http://downloads.sourceforge.net/project/docutils/docutils/0.9.1/', 
            self.python_build
          ) ,
          (
            'mpich2', 
            [], 
            '1.5', 
            'mpich2-', '.tar.gz', 
            'http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.5/', 
            self.mpich2_build
          ) ,
          (
            'mpi4py', 
            ['mpich2'], 
            '1.3', 
            'mpi4py-', '.tar.gz', 
            'http://mpi4py.googlecode.com/files/', 
            self.python_build
          ) ,
          #('openmpi', [], '1.3.3', 'openmpi-', '.tar.gz', 'http://www.open-mpi.org/software/ompi/v1.3/downloads/', self.openmpi_build) ,
          #('setuptools', [], '0.6c11', 'setuptools-', '-py2.6.egg', 'http://pypi.python.org/packages/2.6/s/setuptools/', self.setuptools_install) ,
          #http://pypi.python.org/packages/2.6/s/setuptools/setuptools-0.6c11-py2.6.egg#md5=bfa92100bd772d5a213eedd356d64086
          (
            'fftw3' ,                  #name to refer by
            [],                        #names of prerequisites (unused)
            '3.2.2' ,                  #version string
            'fftw-', '.tar.gz',        #pre- and postfix for filename
            'http://www.fftw.org/',    #download url, filename is appended
            self.fftw_build            #method to use for building
          ) ,
          (
            'gsl' ,                    #name to refer by
            [],                        #names of prerequisites (unused)
            '1.14' ,                   #version string
            'gsl-', '.tar.gz',         #pre- and postfix for filename
            'http://ftp.gnu.org/gnu/gsl/', #download url, filename is appended
            self.fftw_build            #method to use for building - same as for FFTW should work
          ) ,
          (
            'cmake' ,                   #name to refer by
            [],                         #names of prerequisites (unused)
            '2.8.10' ,                   #version string
            'cmake-', '.tar.gz',        #pre- and postfix for filename
            'http://www.cmake.org/files/v2.8/', #download url, filename is appended
            self.cmake_build             #method to use for building - same as for FFTW should work
          ) ,
          (
            'freetype' ,                   #name to refer by
            [],                         #names of prerequisites (unused)
            '2.4.10' ,                   #version string
            'freetype-', '.tar.gz',        #pre- and postfix for filename
            'http://download.savannah.gnu.org/releases/freetype/', #download url, filename is appended
            self.basic_build             #method to use for building - same as for FFTW should work
          ) ,
          (
            'zlib' ,                   #name to refer by
            [],                         #names of prerequisites (unused)
            '1.2.7' ,                   #version string
            'zlib-', '.tar.gz',        #pre- and postfix for filename
            'http://zlib.net/', #download url, filename is appended
            self.basic_build             #method to use for building - same as for FFTW should work
          ) ,
          (
            'png' ,                   #name to refer by
            [],                         #names of prerequisites (unused)
            '1.5.11' ,                   #version string
            'libpng-', '.tar.gz',        #pre- and postfix for filename
            'http://downloads.sourceforge.net/project/libpng/libpng15/older-releases/1.5.11/', #download url, filename is appended
            self.basic_build             #method to use for building - same as for FFTW should work
          ) ,
          #(
          #  'tcl' ,                   #name to refer by
          #  [],                         #names of prerequisites (unused)
          #  '8.5.11' ,                   #version string
          #  'tcl', '-src.tar.gz',        #pre- and postfix for filename
          #  'http://downloads.sourceforge.net/project/tcl/Tcl/8.5.11/', #download url, filename is appended
          #  self.tcl_build             #method to use for building - same as for FFTW should work
          #) ,
          #(
          #  'tk' ,                   #name to refer by
          #  [],                         #names of prerequisites (unused)
          #  '8.5.11' ,                   #version string
          #  'tk', '-src.tar.gz',        #pre- and postfix for filename
          #  'http://downloads.sourceforge.net/project/tcl/Tcl/8.5.11/', #download url, filename is appended
          #  self.tk_build             #method to use for building - same as for FFTW should work
          #) ,
          (
            'zmq' ,                   #name to refer by
            [],                         #names of prerequisites (unused)
            '3.2.1' ,                   #version string
            'zeromq-', '-rc2.tar.gz',        #pre- and postfix for filename
            'http://download.zeromq.org/', #download url, filename is appended
            self.basic_build             #method to use for building - same as for FFTW should work
          ) ,
          
        ]
        
    @late
    def temp_dir(self):
        return os.path.normpath(os.path.join(self.prefix,'..','install'))
    
    
    @late
    def fortran90_compiler(self):
        if 'FC' in os.environ:
            return os.environ['FC']
        else:
            return self.fortran_compiler
            
    @late
    def fortran77_compiler(self):
        if 'F77' in os.environ:
            return os.environ['F77']
        else:
            return None
    
    
    @late
    def fortran_compiler(self):
        if 'FC' in os.environ:
            return os.environ['FC']
        elif 'FORTRAN' in os.environ:
            return os.environ['FORTRAN']
        elif 'F77' in os.environ:
            return os.environ['F77']
        elif 'FORT' in os.environ:
            return os.environ['FORT']
        else:
            return None
            
    @late
    def use_hydra_process_manager(self):
        return False
        
    def setup_temp_dir(self):
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
    
    def run_application(self, args, cwd, env = None):
        print "starting " , ' '.join(args)
        process = subprocess.Popen(args, cwd=cwd, env = env)
        returncode = process.wait()
        if returncode != 0:
            commandline = ' '.join(args)
            raise Exception("Error when running <" + commandline + ">")
        print "finished " , ' '.join(args)
    
    def h5py_build(self, path):
        
        self.run_application([sys.executable,'setup.py','build','--hdf5='+self.prefix], cwd=path)
        self.run_application([sys.executable,'setup.py','install'], cwd=path)
        
    def setuptools_install(self, path):
        self.run_application(['sh',], cwd=path)
        
    def hdf5_build(self, path):
        commands = []
        commands.append([
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared', 
          '--enable-production',
          '--without-zlib',
          '--with-pthread=/usr', 
          '--enable-threadsafe',
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
    
    def python_build(self, path):
        self.run_application([sys.executable,'setup.py','build'], cwd=path)
        self.run_application([sys.executable,'setup.py','install'], cwd=path)
    
    def numpy_build(self, path):
        env = os.environ.copy()
        env['BLAS'] = 'None'
        env['LAPACK'] = 'None'
        env['ATLAS'] = 'None'
        self.run_application([sys.executable,'setup.py','build'], cwd=path, env=env)
        self.run_application([sys.executable,'setup.py','install'], cwd=path, env=env)
        
    def mercurial_build(self, path):
        self.run_application(['make','install','PREFIX='+self.prefix], cwd=path)
    
    def openmpi_build(self, path):
        commands = []
        commands.append([
          './configure','--prefix='+self.prefix,
          #'--enable-mpi-threads', 
          '--enable-cxx-exceptions',
          '--enable-debug',
          '--enable-orterun-prefix-by-default',
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
            
    def mpich2_build(self, path):
        commands = []
        env = os.environ.copy()
        env['CPP'] = '/usr/bin/cpp'
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared',
          '--enable-sharedlibs=gcc',
          '--enable-fc', 
          '--with-python='+sys.executable,
          '--with-device=ch3:sock',
        ]
        if self.use_hydra_process_manager:
            command.append('--with-pm=hydra:mpd:gforker')
        else:
            command.append('--with-pm=gforker')
        if not self.fortran90_compiler is None:
            command.append('FC=' + self.fortran90_compiler)
        
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
            
        self.check_mpich2_install(commands, path)
        
    def fftw_build(self, path):
        commands = []
        command = [
          './configure',
          '--prefix='+self.static_prefix,
          '--disable-shared',
          '--enable-threads'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
            
    
    def cmake_build(self, path):
        commands = []
        command = [
          './configure',
          '--prefix='+self.static_prefix
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
            
    def basic_build(self, path):
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
            
    def tcl_build(self, path):
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared',
          '--enable-threads'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        path = os.path.join(path, 'unix')
        for x in commands:
            self.run_application(x, path)
            
    def tk_build(self, path):            
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared',
          '--disable-xss',
          '--enable-threads'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        path = os.path.join(path, 'unix')
        for x in commands:
            self.run_application(x, path)
            
        
    def check_mpich2_install(self, commands, path):
        bin_directory = os.path.join(self.prefix, 'bin')
        mpif90_filename = os.path.join(bin_directory, 'mpif90')
        if not os.path.exists(mpif90_filename):
            print "-----------------------------------------------------------------"
            print "MPICH build incomplete, no fortran 90 support"
            print "-----------------------------------------------------------------"
            print "The 'mpif90' command was not build"
            print "This is usually caused by an incompatible C and fortran compiler"
            print "Please set the F90, F77 and CC environment variables"
            print 
            print "After changing the environment variables,"
            print "you can restart the install with:"
            print
            print "  ./install.py install mpich2 mpi4py"
            print
            print "You can rerun the build by hand, using:"
            print 
            print "  cd", path
            for command in commands:
                print
                if len(command) < 3:
                    print ' ', ' '.join(command)
                else:
                    print '   \\\n    '.join(command)
            
            sys.exit(1)
            
            
    def openssl_build(self, path):
        commands = []
        commands.append([
          './config','--prefix='+self.prefix,
          ' --openssldir='+self.prefix+'/openssl',
          '--shared'
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)

    def download_apps(self, names, skip):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            if skip and name in skip:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            url = url_prefix + app_file
            temp_app_file = os.path.join(self.temp_dir, app_file)
            if not os.path.exists(temp_app_file):
                print "Downloading ", app_file
                urllib.urlretrieve(url, os.path.join(self.temp_dir, app_file))
                print "...Finished"
                
    
    def list_apps(self, names, skip):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if skip and name in skip:
                continue
            print name, " - dowloaded from", url_prefix
                
    def unpack_apps(self, names, skip):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            if skip and name in skip:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            url = url_prefix + app_file
            temp_app_file = os.path.join(self.temp_dir, app_file)
            temp_app_dir = os.path.join(self.temp_dir , app_dir)
            if os.path.exists(temp_app_dir):
                shutil.rmtree(temp_app_dir)
            
            print "Unpacking ", app_file
            try:
                self.run_application(['tar','-xzf',app_file], cwd=self.temp_dir)
            except:
                print "----------------------------------------------------------"
                print "Could not unpack source file of", name
                print "----------------------------------------------------------"
                print
                print "Download location may have changed"
                print "Please download the source file yourself, "
                print "or contact the AMUSE development team."
                print "http://castle.strw.leidenuniv/trac/amuse"
                print
                print "To download the file you can update the URL in"
                print "one of the following lines and run the command."
                print
                print "curl ", url, "-o", temp_app_file
                print
                print "or" 
                print
                print "wget ", url, "-O", temp_app_file
                print
                print "Note: The name of the output file must not be changed (after the -o or -O parameter)"
                sys.exit(1)
            print "...Finished"
            
    def build_apps(self, names, skip):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            if skip and name in skip:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            temp_app_dir = os.path.join(self.temp_dir , app_dir)
            if not os.path.exists(temp_app_dir):
                if prefix.endswith('-'):
                    app_dir = prefix[:-1]
                else:
                    app_dir = prefix
                temp_app_dir = os.path.join(self.temp_dir , app_dir)
                if not os.path.exists(temp_app_dir):
                    app_file = prefix + version + suffix
                    if app_file.endswith('.tar.gz'):
                        app_dir = app_file[:-len('.tar.gz')]
                    elif app_file.endswith('.tar.bz2'):
                        app_dir = app_file[:-len('.tar.bz2')]
                    else:
                        app_dir, ignore = os.path.os.path.splitext(app_file)
                        
                    temp_app_dir = os.path.join(self.temp_dir , app_dir)
                    if not os.path.exists(temp_app_dir):
                        print "Package was not correctly unpacked: ", app_file
                        return
    
            print "Building ", app_file
            function(temp_app_dir)
            print "...Finished"
            
class InstallPrerequisitesOnOSX(InstallPrerequisites):   
    @late
    def prefix(self):
        path = os.path.split(sys.executable)[0]
        if 'Resources' in path:
            return path[:path.index('Resources')]
        else:
            return sys.prefix
            
    @late
    def temp_dir(self):
        return os.path.normpath(os.path.abspath('install'))
        
    @late
    def static_prefix(self):
        return os.path.normpath(os.path.abspath('static_libs'))
    

    def mpich2_build(self, path):
        
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-fc',
          '--enable-shared',
          '--with-python='+sys.executable,
          '--with-pm=gforker',
          '--enable-sharedlibs=osx-gcc',
          '--with-device=ch3:sock',
        ]
        if self.use_hydra_process_manager:
            command.append('--with-pm=hydra:mpd:gforker')
        else:
            command.append('--with-pm=gforker')
            
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
        
        self.check_mpich2_install(commands, path)
        
    def tcl_build(self, path):
        if True:
            return
            
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared',
          '--enable-threads'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        path = os.path.join(path, 'unix')
        for x in commands:
            self.run_application(x, path)
            
    def tk_build(self, path):
        if True:
            return
            
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared',
          '--disable-xss',
          '--enable-threads'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        path = os.path.join(path, 'unix')
        for x in commands:
            self.run_application(x, path)
            
            
     
if IS_ON_OSX:
    INSTALL = InstallPrerequisitesOnOSX()
else:
    INSTALL = InstallPrerequisites()
           
def download(names, skip):
    INSTALL.download_apps(names, skip)

def install(names, skip):
    INSTALL.download_apps(names, skip)
    INSTALL.unpack_apps(names, skip)
    INSTALL.build_apps(names, skip)
    
def list(names, skip):
    INSTALL.list_apps(names, skip)

_commands = {
    'download' : download,
    'install' : install,
    'list' : list,
}

if __name__ == '__main__':
    print "Files are installed in: ", INSTALL.prefix
    print "Files are downloaded to: ", INSTALL.temp_dir
    print ""
    
    if INSTALL.fortran90_compiler is None or INSTALL.fortran77_compiler is None:
        print """No fortran 90 compiler environment variable set.
A FORTRAN 90 compiler is needed for MPI and several modules, 
please set FC and F77 first by (bash, replace gfortran with your preferred
compiler):

export FC=gfortran
export F77=gfortran

or (csh):

setenv FC gfortran 
setenv F77 gfortran 

"""
        sys.exit(1)
    else:
        print "Fortran 90 compiler used will be: ", INSTALL.fortran90_compiler
        print "Fortran 77 compiler used will be: ", INSTALL.fortran77_compiler
    
    print ""
    INSTALL.setup_temp_dir()
    do = []
    names = []
    flag = False
    skip = []
    
    for x in sys.argv:
        if x in _commands.keys():
            do.append(x)
            flag = True
        else:
            if x == '--hydra':
                INSTALL.use_hydra_process_manager = True
            if flag:
                if x.startswith('no-'):
                    skip.append(x[3:])
                else:
                    names.append(x)
    for x in do:
        _commands[x](names, skip)
    
    if len(do) == 0:
        print "Usage: install.py download|install|list [package names]"
        print ""
        print "download  download the packages to the download directory"
        print "install   unpack and install the packages to the prefix directory"
        print ""
        print "you can also install download or install individual packages"
        print "please specify a list of packages to install"
        print ""
        print "to install all prerequisites do:"
        print ""
        print "./install.py install"
        print ""
        print "to get a list of all packages:"
        print ""
        print "./install.py list"
        print ""
        print "to install mpich2 with the hydra process manager do:"
        print ""
        print "./install.py --hydra install"
        print ""
        
        
        sys.exit(1)
    

        
