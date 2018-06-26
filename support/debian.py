from __future__ import print_function

import os
import shutil
import subprocess
try:  # Python 3
    from urllib.request import urlretrieve
except ImportError:  # Python 2
    from urllib import urlretrieve
import platform
import sys
from optparse import OptionParser


depends = \
[
    'bash (>= 2.05a-11)',
    'build-essential',
    'cmake',
    'gfortran',
    'python2.6 (>=2.6.0)',
    'python-numpy (>=1.3.0)',
    'python-nose (>=0.11)',
    'python-matplotlib (>=0.99)',
    'python-setuptools',
    'python-docutils',
    'python-h5py (>=1.2.1)',
    'libhdf5-serial-dev (>=1.6)',
    'hdf5-tools',
    'libfftw3-3', 
    'libfftw3-dev', 
    'libfftw3-doc',
    'libopenmpi-dev (>=1.4.1)',
    'openmpi-bin',
    'libgsl0ldbl',
    'libgsl0-dev',
    
]
control = \
"""
Package: amuse
Version: {0}
Section: base
Priority: optional
Architecture: all
Depends: {1} 
Maintainer: amuseteam <info@amusecode.org>
Description: Astrophysical Multipurpose Software Environment
 A software framework for large-scale simulations of dense 
 stellar systems, in which existing codes for dynamics, 
 stellar evolution, and hydrodynamics can be easily coupled. 
"""
postinst = \
"""#!/bin/sh -e
easy_install mpi4py
"""


class generate_debian_package(object):
    
    def __init__(self, version = 'svn', arch = None):
        self.version = version
        
        if self.version == 'svn':
            self.version = 'r{0}'.format(self.get_svn_revision())
            
        self.amuse_version = 'amuse-{0}'.format(self.version)
        self.debianversion = '0ubuntu'
        
        if arch is None:
            self.architecture = 'i386' if platform.architecture()[0] == '32bit' else 'amd64'
        else:
            self.architecture = arch
            
        self.package_name = 'amuse_{0}-{1}-{2}'.format(self.version, self.debianversion, self.architecture)
        self.package_path = os.path.abspath(self.package_name)
    
    def get_svn_revision(self):
        stdoutstring, stderrstring = subprocess.Popen(
            ['svn','info'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        ).communicate()
        lines = stdoutstring.splitlines()
        revision = '0000'
        for x in lines:
            if x.startswith('Revision:'):
                label, revision = x.split(': ')
        return revision
        
        
    def run(self):
        self.setup_deb_builddir()
        self.makescripts()
        self.install_in_deb_builddir()
        self.package()
        self.cleanup()
        
        print("generated debian package: {0}.deb".format(self.package_name))

    def makescripts(self):
        #os.system('python setup.py generate_main --amuse-dir=/usr/share/{0}'.format(amuse_version))
        pass
        
    def setup_deb_builddir(self):
        if os.path.exists(self.package_path):
            shutil.rmtree(self.package_path)
        os.makedirs(self.package_path)
        
    def install_in_deb_builddir(self):
        debian_path = os.path.join(self.package_path, 'DEBIAN')
        os.makedirs(debian_path)
        
        
        dependency_string = ', '.join(depends)
        with open(os.path.join(debian_path, 'control'),'w') as f:
            f.write(control.format(self.version, dependency_string))

        #f = open(os.path.join(debian_path, 'postinst'),'w')
        #f.write(postinst)
        #f.close()
        #os.chmod(os.path.join(debian_path, 'postinst'), 0b11110  
        
        if not os.path.exists('build'):
            os.makedirs('build')
            
        mpi4pyfile = 'mpi4py-1.2.2.tar.gz'
        urlretrieve(
            'http://mpi4py.googlecode.com/files/{0}'.format(mpi4pyfile),
            mpi4pyfile
        )
        shutil.copyfile(mpi4pyfile, os.path.join('build', mpi4pyfile))
        
        subprocess.call(
            [
            'tar',
            '-xf',
            mpi4pyfile
            ]
            ,
            cwd=os.path.join('build')
        )
        subprocess.call(
            [
            sys.executable,
            'setup.py',
            'install',
            '--prefix=/usr',
            '--install-layout=deb',
            '--root={0}'.format(self.package_path),
            ]
            ,
            cwd=os.path.join('build','mpi4py-1.2.2')
        )
        
        subprocess.call(
            [
            sys.executable,
            'setup.py',
            'install',
            '--prefix=/usr',
            '--install-layout=deb',
            '--root={0}'.format(self.package_path),
            ]
        )
            
        #shutil.copytree('src', './{0}/usr/share/{1}/src'.format(package_name, amuse_version))
        #shutil.copytree('test', './{0}/usr/share/{1}/test'.format(package_name, amuse_version))
        #shutil.copytree('data', './{0}/usr/share/{1}/data'.format(package_name, amuse_version))
        #shutil.copytree('lib', './{0}/usr/share/{1}/lib'.format(package_name, amuse_version))
        #shutil.copytree('doc', './{0}/usr/share/doc/{1}/doc'.format(package_name, amuse_version))
        #self.touch('./{0}/usr/share/{1}/build.py'.format(package_name, amuse_version))

        #shutil.copy('amuse.sh', './{0}/usr/bin'.format(package_name))
        #shutil.copy('iamuse.sh', './{0}/usr/bin'.format(package_name))

        #os.chmod('./{0}/usr/bin/amuse.sh'.format(package_name), 0b111101101)
        #os.chmod('./{0}/usr/bin/iamuse.sh'.format(package_name), 0b111101101)

    def package(self):
        if os.path.exists(self.package_path):
            print("creating debian package..")
            subprocess.call(
                [
                'fakeroot',
                'dpkg-deb',
                '--build',
                self.package_name
                ]
            )
            
    def cleanup(self):
        
        #if os.path.exists(self.package_path):
        #    shutil.rmtree(self.package_path)
            
        #os.system('python setup.py generate_main')
        pass
          
    def touch(self, filename):
        if not os.path.exists(filename):
            open(filename, 'w').close()
            
    
def main(version = 'svn', arch = None):
    command = generate_debian_package(version, arch)
    command.run()
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-v", "--version", 
        default = 'svn',
        dest="version",
        help="version of the debian package to create, defaults to svn based",
        type="string"
    )    
    result.add_option(
        "-a", "--arch", 
        dest="arch",
        default = None,
        help="architecture to build (i386 or amd64)",
        type="string"
    )    
    return result
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)
    
   
