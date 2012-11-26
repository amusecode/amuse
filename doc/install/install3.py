#!/usr/bin/env python



import sys
import os.path
import os
import urllib.request, urllib.parse, urllib.error
import subprocess
import shutil


IS_ON_OSX = sys.platform == 'darwin'
PYTHON = sys.executable

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
        path = os.path.split(sys.executable)[0]
        if 'Framework' in path:
            return path[:path.index('Framework')]
        else:
            return path[:path.index('bin')-1]
    @late
    def applications(self):
       return [
       #numpy-1.7.0b2.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fnumpy%2Ffiles%2FNumPy%2F1.7.0b2%2F&ts=1353670630&use_mirror=garr
          #('openssl' , [], '0.9.8k' , 'openssl-', '.tar.gz', 'http://www.openssl.org/source/', self.openssl_build),
          (
            'numpy' ,                  #name to refer by
            [],                        #names of prerequisites (unused)
            '1.7.0b2' ,                  #version string
            'numpy-', '.tar.gz',       #pre- and postfix for filename
            #'http://pypi.python.org/packages/source/n/numpy/', #download url, filename is appended
            'http://downloads.sourceforge.net/project/numpy/NumPy/1.7.0b2/',
            self.numpy_build          #method to use for building
          ),
          (
            'nose', 
            [], 
            '1.2.1', 
            'nose-' , '.tar.gz', 
            'http://pypi.python.org/packages/source/n/nose/', 
            self.python_build
          ),
          (
            'hdf' ,
            [],  
            '1.8.8',
            'hdf5-' , '.tar.gz' , 
            'http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.8/src/',
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
            'http://pypi.python.org/packages/source/d/docutils/', 
            self.python_build
          ) ,
          (
            'mpich2', 
            [], 
            '1.5', 
            'mpich2-', '.tar.gz', 
            'http://www.mpich.org/static/tarballs/1.5/', 
            self.mpich2_build
          ) ,
          (
            'mpi4py', 
            ['mpich2'], 
            '1.2.2', 
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
            '2.8.8' ,                   #version string
            'cmake-', '.tar.gz',        #pre- and postfix for filename
            'http://www.cmake.org/files/v2.8/', #download url, filename is appended
            self.cmake_build             #method to use for building
          ) ,
          (
            'gmp',                      #name to refer by
            [],                         #names of prerequisites (unused)
            '5.0.3' ,                   #version string
            'gmp-', '.tar.bz2',        #pre- and postfix for filename
            'ftp://ftp.gmplib.org/pub/gmp-5.0.3/', #download url, filename is appended
            self.gmp_build             #method to use for building
          ) ,
          ( # NOTE: When library version is changed, url to 'allpatches' in self.mpfr_build must be changed too!
            'mpfr' ,                    #name to refer by
            ['gmp'],                    #names of prerequisites
            '3.1.1' ,                   #version string
            'mpfr-', '.tar.gz',         #pre- and postfix for filename
            'http://mpfr.loria.fr/mpfr-3.1.1/', #download url, filename is appended
            self.mpfr_build             #method to use for building
          ) ,
        ]
        
    @late
    def temp_dir(self):
        return os.path.join(self.prefix,'install','_temp')
    
    
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
           
    @late
    def use_gforker_process_manager(self):
        return False
        
    def setup_temp_dir(self):
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
    
    def run_application(self, args, cwd, env = None):
        print("starting " , ' '.join(args))
        process = subprocess.Popen(args, cwd=cwd, env = env)
        returncode = process.wait()
        if returncode != 0:
            commandline = ' '.join(args)
            raise Exception("Error when running <" + commandline + ">")
        print("finished " , ' '.join(args))
    
    def h5py_build(self, path):
        self.run_application([
            'patch', 
            '-p1', 
            '-i', 
            os.path.abspath('h5py_imports_python33.patch')
        ], cwd=path)
       
        self.run_application([PYTHON,'setup.py','build','--hdf5='+self.prefix], cwd=path)
        self.run_application([PYTHON,'setup.py','install'], cwd=path)
        
    def setuptools_install(self, path):
        self.run_application(['sh',], cwd=path)
        
    def hdf5_build(self, path):
        commands = []
        commands.append([
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared', 
          '--enable-production',
          '--with-pthread=/usr', 
          '--enable-threadsafe'
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
    
    def python_build(self, path):
        self.run_application([PYTHON,'setup.py','build'], cwd=path)
        self.run_application([PYTHON,'setup.py','install'], cwd=path)
    
    def numpy_build(self, path):
        env = os.environ.copy()
        env['BLAS'] = 'None'
        env['LAPACK'] = 'None'
        env['ATLAS'] = 'None'
        self.run_application([PYTHON,'setup.py','build'], cwd=path, env=env)
        self.run_application([PYTHON,'setup.py','install'], cwd=path, env=env)
        
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
        elif self.use_gforker_process_manager:
            command.append('--with-pm=gforker:hydra:mpd')
        else:
            command.append('--with-pm=mpd:hydra:gforker')
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
          '--prefix='+self.prefix,
          '--enable-shared',
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
          '--prefix='+self.prefix
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
    
    def gmp_build(self, path):
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-shared'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'check'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
    
    def mpfr_build(self, path):
        #temp_patch_file = os.path.join(self.temp_dir, "mpfr-allpatches")
        #if not os.path.exists(temp_patch_file):
        #    print "Downloading mpfr-allpatches"
        #    urllib.urlretrieve("http://www.mpfr.org/mpfr-3.1.1/allpatches", temp_patch_file)
        #    print "...Finished"
        
        commands = []
        #commands.append(['patch', '-N', '-Z', '-p1', '-i', temp_patch_file])
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--with-gmp='+self.prefix,
          '--enable-shared',
          '--enable-thread-safe'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'check'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)
    
    def check_mpich2_install(self, commands, path):
        bin_directory = os.path.join(self.prefix, 'bin')
        mpif90_filename = os.path.join(bin_directory, 'mpif90')
        if not os.path.exists(mpif90_filename):
            print("-----------------------------------------------------------------")
            print("MPICH build incomplete, no fortran 90 support")
            print("-----------------------------------------------------------------")
            print("The 'mpif90' command was not build")
            print("This is usually caused by an incompatible C and fortran compiler")
            print("Please set the F90, F77 and CC environment variables")
            print() 
            print("After changing the environment variables,")
            print("you can restart the install with:")
            print()
            print("  ./install.py install mpich2 mpi4py")
            print()
            print("You can rerun the build by hand, using:")
            print() 
            print("  cd", path)
            for command in commands:
                print()
                if len(command) < 3:
                    print(' ', ' '.join(command))
                else:
                    print('   \\\n    '.join(command))
            
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
                print("Downloading ", app_file, url)
                urllib.request.urlretrieve(url, os.path.join(self.temp_dir, app_file))
                print("...Finished")
                
    
    def list_apps(self, names, skip):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if skip and name in skip:
                continue
            print(name, " - dowloaded from", url_prefix)
                
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
            
            print("Unpacking ", app_file)
            try:
                self.run_application(['tar','-xf',app_file], cwd=self.temp_dir)
            except:
                print("----------------------------------------------------------")
                print("Could not unpack source file of", name)
                print("----------------------------------------------------------")
                print()
                print("Download location may have changed")
                print("Please download the source file yourself, ")
                print("or contact the AMUSE development team.")
                print("http://castle.strw.leidenuniv/trac/amuse")
                print()
                print("To download the file you can update the URL in")
                print("one of the following lines and run the command.")
                print()
                print("curl ", url, "-o", temp_app_file)
                print()
                print("or") 
                print()
                print("wget ", url, "-O", temp_app_file)
                print()
                print("Note: The name of the output file must not be changed (after the -o or -O parameter)")
                sys.exit(1)
            print("...Finished")
            
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
                        print("Package was not correctly unpacked: ", app_file)
                        return
    
            print("Building ", app_file)
            function(temp_app_dir)
            print("...Finished")
            
class InstallPrerequisitesOnOSX(InstallPrerequisites):


   
          
    def mpich2_build(self, path):
        
        commands = []
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--enable-fc',
          '--enable-shared',
          '--with-python='+sys.executable,
          '--with-pm=mpd',
          '--enable-sharedlibs=osx-gcc',
          '--with-device=ch3:sock',
        ]
        if self.use_hydra_process_manager:
            command.append('--with-pm=hydra:mpd:gforker')
        elif self.use_gforker_process_manager:
            command.append('--with-pm=gforker:hydra:mpd')
        else:
            command.append('--with-pm=mpd:hydra:gforker')
            
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
        
        self.check_mpich2_install(commands, path)
        
    def mpfr_build(self, path):
        #temp_patch_file = os.path.join(self.temp_dir, "mpfr-allpatches")
        #if not os.path.exists(temp_patch_file):
        #    print "Downloading mpfr-allpatches"
        #    urllib.urlretrieve("http://www.mpfr.org/mpfr-3.1.1/allpatches", temp_patch_file)
        #    print "...Finished"
        
        commands = []
        #commands.append(['patch', '-N', '-Z', '-p1', '-i', temp_patch_file])
        command = [
          './configure',
          '--prefix='+self.prefix,
          '--with-gmp='+self.prefix,
          '--enable-shared'
        ]
        commands.append(command)
        commands.append(['make'])
        commands.append(['make', 'check'])
        commands.append(['make', 'install'])
        
        for x in commands:
            self.run_application(x, path)

class InstallMatplotlib(InstallPrerequisites):
    
    @late
    def applications(self):
        return (
               (
                'freetype' ,                   #name to refer by
                [],                         #names of prerequisites (unused)
                '2.4.9' ,                   #version string
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
              ),
              (
                'matplotlib' ,                   #name to refer by
                [],                         #names of prerequisites (unused)
                '1.1.0' ,                   #version string
                'matplotlib-', '.tar.gz',        #pre- and postfix for filename
                ' http://pypi.python.org/packages/source/m/matplotlib/', #download url, filename is appended
                self.matplotlib_build             #method to use for building - same as for FFTW should work
              ),
        )
        
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
            
    def matplotlib_build(self, path):
        env = os.environ.copy()
        env['CFLAGS'] ="-I{0}/include -I{0}/include/freetype2".format(self.prefix)
        env['LDFLAGS'] = "-L{0}/lib".format(self.prefix)
        self.run_application([PYTHON,'setup.py','build'], cwd=path, env = env)
        self.run_application([PYTHON,'setup.py','install'], cwd=path, env = env)
      
     
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
    
def listpackages(names, skip):
    INSTALL.list_apps(names, skip)

_commands = {
    'download' : download,
    'install' : install,
    'list' : listpackages,
}

if __name__ == '__main__':
    print("")
    
    if INSTALL.fortran90_compiler is None or INSTALL.fortran77_compiler is None:
        print("""No fortran 90 compiler environment variable set.
A FORTRAN 90 compiler is needed for MPI and several modules, 
please set FC and F77 first by (bash, replace gfortran with your preferred
compiler):

export FC=gfortran
export F77=gfortran

or (csh):

setenv FC gfortran 
setenv F77 gfortran 

""")
        sys.exit(1)
    else:
        print("Fortran 90 compiler used will be: ", INSTALL.fortran90_compiler)
        print("Fortran 77 compiler used will be: ", INSTALL.fortran77_compiler)
    
    print("")
    do = []
    names = []
    flag = False
    skip = []
    
    for x in sys.argv:
        if x in list(_commands.keys()):
            do.append(x)
            flag = True
        else:
            if x == '--matplotlib':
                INSTALL = InstallMatplotlib()
                print("----------------------------------------------------------")
                print("This feature is optional and experimental!")
                print("----------------------------------------------------------")
                print()
                print("Will download and install matplotlib")
                print("plus it most common prerequisites (zlib, png and freetype)")
            elif x == '--hydra':
                INSTALL.use_hydra_process_manager = True
            elif x == '--gforker':
                INSTALL.use_gforker_process_manager = True
            elif x.startswith('--prefix='):
                INSTALL.prefix = x[len('--prefix='):]
            elif flag:
                if x.startswith('no-'):
                    skip.append(x[3:])
                else:
                    names.append(x)

    print("Files are installed in: ", INSTALL.prefix)
    print("Files are downloaded to: ", INSTALL.temp_dir)

    INSTALL.setup_temp_dir()

    for x in do:
        _commands[x](names, skip)
    
    if len(do) == 0:
        print("Usage: install.py download|install|list [package names]")
        print("")
        print("download  download the packages to the download directory")
        print("install   unpack and install the packages to the prefix directory")
        print("")
        print("you can also install download or install individual packages")
        print("please specify a list of packages to install")
        print("")
        print("to install all prerequisites do:")
        print("")
        print("./install.py install")
        print("")
        print("to get a list of all packages:")
        print("")
        print("./install.py list")
        print("")
        print("hydra is a modern process manager with more options and")
        print("faster distributed process start-up")
        print("to install mpich2 with the hydra process manager do:")
        print("")
        print("./install.py --hydra install")
        print("")
        print("the gforker process manager is easier to run (no daemon)")
        print("but only works on the local machine") 
        print("to install mpich2 with the gforker process manager do:")
        print("")
        print("./install.py --gforker install")
        print("")
        
        
        sys.exit(1)
    

        
