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
        path = os.path.split(sys.executable)[0]
        if 'Framework' in path:
            return path[:path.index('Framework')]
        else:
            return path[:path.index('bin')-1]
    @late
    def applications(self):
       return [
          ('openssl' , [], '0.9.8k' , 'openssl-', '.tar.gz', 'http://www.openssl.org/source/', self.openssl_build),
          ('numpy' , [], '1.3.0' , 'numpy-', '.tar.gz', 'http://superb-east.dl.sourceforge.net/sourceforge/numpy/', self.python_build),
          ('nose', [], '0.11.1', 'nose-' , '.tar.gz', 'http://somethingaboutorange.com/mrl/projects/nose/', self.python_build),
          ('hdf' , [],  '1.8.3' , 'hdf5-' , '.tar.gz' , 'http://www.hdfgroup.org/ftp/HDF5/current/src/', self.hdf5_build) ,
          ('h5py', ['hdf'], '1.2.0', 'h5py-' , '.tar.gz', 'http://h5py.googlecode.com/files/', self.h5py_build) ,
          ('mpich2', [], '1.1', 'mpich2-', '.tar.gz', 'http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.1/', self.mpich2_build) ,
          ('mpi4py', ['mpich2'], '1.0.0', 'mpi4py-', '.tar.gz', 'http://mpi4py.googlecode.com/files/', self.python_build) ,
        ]
    @late
    def temp_dir(self):
        return os.path.join(self.prefix,'install','_temp')
    def setup_temp_dir(self):
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
    def run_application(self, args, cwd):
        print "starting " , ' '.join(args)
        subprocess.Popen(args, cwd=cwd).wait()
        print "finished " , ' '.join(args)
    def h5py_build(self, path):
        self.run_application(['python','setup.py','configure', '--hdf5='+self.prefix, '--api=18'], cwd=path)
        self.python_build(path)
    def hdf5_build(path):
        commands = []
        commands.append([
          './configure','--prefix='+self.prefix,
          '--enable-shared', '--enable-production',
          '--with-pthread=/usr', '--enable-threadsafe'
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
    def python_build(self, path):
        self.run_application(['python','setup.py','build'], cwd=path)
        self.run_application(['python','setup.py','install'], cwd=path)
    def mercurial_build(self, path):
        self.run_application(['make','install','PREFIX='+self.prefix], cwd=path)
    def mpich2_build(self, path):
        commands = []
        commands.append([
          './configure','--prefix='+self.prefix,
          '--enable-sharedlibs=gcc',
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
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
    def download_apps(self, names):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            url = url_prefix + app_file
            temp_app_file = os.path.join(self.temp_dir, app_file)
            if not os.path.exists(temp_app_file):
                print "Downloading ", app_file
                urllib.urlretrieve(url, os.path.join(self.temp_dir, app_file))
                print "...Finished"
    def unpack_apps(self, names):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            url = url_prefix + app_file
            temp_app_file = os.path.join(self.temp_dir, app_file)
            temp_app_dir = os.path.join(self.temp_dir , app_dir)
            if os.path.exists(temp_app_dir):
                shutil.rmtree(temp_app_dir)
            
            print "Unpacking ", app_file
            self.run_application(['tar','-xf',app_file], cwd=self.temp_dir)
            print "...Finished"
    def build_apps(self, names):
        for (name, dependencies, version, prefix, suffix, url_prefix, function) in self.applications:
            if names and name not in names:
                continue
            app_file = prefix + version + suffix
            app_dir = prefix + version 
            temp_app_dir = os.path.join(self.temp_dir , app_dir)
            if not os.path.exists(temp_app_dir):
                print "Package was not correctly unpacked: ", app_file
                return
    
            print "Building ", app_file
            function(temp_app_dir)
            print "...Finished"
            
class InstallPrerequisitesOnOSX(InstallPrerequisites):

    def mpich2_build(self, path):
        commands = []
        commands.append([
          './configure','--prefix='+self.prefix,
          '--enable-sharedlibs=osx-gcc',
        ])
        commands.append(['make'])
        commands.append(['make', 'install'])
        for x in commands:
            self.run_application(x, path)
     
if IS_ON_OSX:
    INSTALL = InstallPrerequisitesOnOSX()
else:
    INSTALL = InstallPrerequisites()
           
def download(names):
    INSTALL.download_apps(names)

def install(names):
    INSTALL.download_apps(names)
    INSTALL.unpack_apps(names)
    INSTALL.build_apps(names)

_commands = {
    'download' : download,
    'install' : install
}

if __name__ == '__main__':
    print "files are installed to: ", INSTALL.prefix
    print "files are downloaded to: ", INSTALL.temp_dir
    INSTALL.setup_temp_dir()
    do = []
    names = []
    flag = False
    for x in sys.argv:
        if x in _commands.keys():
            do.append(x)
            flag = True
        else:
            if flag:
                names.append(x)
    for x in do:
    	_commands[x](names)
