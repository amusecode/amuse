#!/usr/bin/env python



import sys
import os.path
import os
import urllib
import subprocess
import shutil

#temp_dir = os.path.join('_downloaded','_temp')
def get_prefix():
    path = os.path.split(sys.executable)[0]
    if 'Framework' in path:
        prefix = path[:path.index('Framework')]
    else:
        prefix = path[:path.index('bin')-1]
    return prefix


prefix = get_prefix()
print "prefix: <", prefix,">"
temp_dir = os.path.join(prefix,'install','_temp')
print "files are downloaded to: " , temp_dir

def setup_temp_dir():
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

def run_application(args, cwd):
    print "starting " , ' '.join(args)
    subprocess.Popen(args, cwd=cwd).wait()
    print "finished " , ' '.join(args)
    
   
    
def h5py_build(path):
    run_application(['python','setup.py','configure', '--hdf5='+prefix, '--api=18'], cwd=path)
    python_build(path)

def python_build(path):
    run_application(['python','setup.py','build'], cwd=path)
    run_application(['python','setup.py','install'], cwd=path)
    
def mercurial_build(path):
    run_application(['make','install','PREFIX='+prefix], cwd=path)

def hdf5_build(path):
    commands = []
    commands.append([
      './configure','--prefix='+prefix,
      '--enable-shared', '--enable-production',
      '--with-pthread=/usr', '--enable-threadsafe'
    ])
    commands.append(['make'])
    commands.append(['make', 'install'])
    for x in commands:
        run_application(x, path)

def mpich2_build(path):
    commands = []
    commands.append([
      './configure','--prefix='+prefix,
      '--enable-sharedlibs=gcc',
    ])
    commands.append(['make'])
    commands.append(['make', 'install'])
    for x in commands:
        run_application(x, path)

def openssl_build(path):
    commands = []
    commands.append([
      './config','--prefix='+prefix,
      ' --openssldir='+prefix+'/openssl',
      '--shared'
    ])
    commands.append(['make'])
    commands.append(['make', 'install'])
    for x in commands:
        run_application(x, path)







apps = [
  ('openssl' , [], '0.9.8k' , 'openssl-', '.tar.gz', 'http://www.openssl.org/source/', openssl_build),
  ('numpy' , [], '1.3.0' , 'numpy-', '.tar.gz', 'http://superb-east.dl.sourceforge.net/sourceforge/numpy/', python_build),
  ('nose', [], '0.11.1', 'nose-' , '.tar.gz', 'http://somethingaboutorange.com/mrl/projects/nose/', python_build),
  ('hdf' , [],  '1.8.3' , 'hdf5-' , '.tar.gz' , 'http://www.hdfgroup.org/ftp/HDF5/current/src/', hdf5_build) ,
  ('h5py', ['hdf'], '1.1.0', 'h5py-' , '.tar.gz', 'http://h5py.googlecode.com/files/', h5py_build) ,
  ('mpich2', [], '1.1', 'mpich2-', '.tar.gz', 'http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.1/', mpich2_build) ,
  ('mpi4py', ['mpich2'], '1.0.0', 'mpi4py-', '.tar.gz', 'http://mpi4py.googlecode.com/files/', python_build) ,
]

def download_apps(apps, names):
    for (name, dependencies, version, prefix, suffix, url_prefix, function) in apps:
        if names and name not in names:
            continue
        app_file = prefix + version + suffix
        app_dir = prefix + version 
        url = url_prefix + app_file
        temp_app_file = os.path.join(temp_dir, app_file)
        if not os.path.exists(temp_app_file):
            print "Downloading ", app_file
            urllib.urlretrieve(url, os.path.join(temp_dir, app_file))
            print "...Finished"


def unpack_apps(apps, names):
    for (name, dependencies, version, prefix, suffix, url_prefix, function) in apps:
        if names and name not in names:
            continue
        app_file = prefix + version + suffix
        app_dir = prefix + version 
        url = url_prefix + app_file
        temp_app_file = os.path.join(temp_dir, app_file)
        temp_app_dir = os.path.join(temp_dir , app_dir)
        if os.path.exists(temp_app_dir):
            shutil.rmtree(temp_app_dir)
        
        print "Unpacking ", app_file
        run_application(['tar','-xf',app_file], cwd=temp_dir)
        print "...Finished"

def build_apps(apps, names):
    for (name, dependencies, version, prefix, suffix, url_prefix, function) in apps:
        if names and name not in names:
            continue
        app_file = prefix + version + suffix
        app_dir = prefix + version 
        temp_app_dir = os.path.join(temp_dir , app_dir)
        if not os.path.exists(temp_app_dir):
            print "Package was not correctly unpacked: ", app_file
            return

        print "Building ", app_file
        function(temp_app_dir)
        print "...Finished"


hg_server = "koppoel.strw.leidenuniv.nl"
hg_path = "data/hg/Ratatouille"
hg_ssh = "ssh://"+hg_server+"//"+hg_path

def get_code():
    args = ['hg', 'clone', hg_ssh]
    run_application(args, '.')

def download(names):
    download_apps(apps, names)
    
def install(names):
    unpack_apps(apps, names)
    build_apps(apps, names)

def update_code():
    pass
    
_commands = {
    'download' : download,
    'install' : install,
    'get' : get_code,
    'update' : update_code
}

if __name__ == '__main__':
    setup_temp_dir()
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
