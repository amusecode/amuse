import os
import shutil
import subprocess
import urllib
import platform
import sys


version = '5.1'
amuse_version = 'amuse-{0}'.format(version)
debianversion = '0ubuntu'
architecture = 'i386' if platform.architecture()[0] == '32bit' else 'amd64'

#f=os.popen('uname -m')
#architecture = f.readlines()[0]
#print architecture
#f.close()

package_name = 'amuse_{0}-{1}-{2}'.format(version, debianversion, architecture)
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
    def __init__(self):
        self.makescripts()
        self.makedirs()
        self.package()
        self.cleanup()

    def makescripts(self):
        os.system('python setup.py generate_main --amuse-dir=/usr/share/{0}'.format(amuse_version))

    def makedirs(self):
        package_path = os.path.abspath(package_name)
        
        if os.path.exists(package_path):
            shutil.rmtree(package_path)

        os.makedirs(package_path)
        
        dependency_string = ', '.join(depends)
        
        debian_path = os.path.join(package_path, 'DEBIAN')
        
        os.makedirs(debian_path)
        f = open(os.path.join(debian_path, 'control'),'w')
        f.write(control.format(version, dependency_string))
        f.close()

        #f = open(os.path.join(debian_path, 'postinst'),'w')
        #f.write(postinst)
        #f.close()
        #os.chmod(os.path.join(debian_path, 'postinst'), 0b11110  
        
        if not os.path.exists('build'):
            os.makedirs('builld')
            
        mpi4pyfile = 'mpi4py-1.2.2.tar.gz'
        urllib.urlretrieve('http://mpi4py.googlecode.com/files/{0}'.format(mpi4pyfile))
        shutil.copyfile(mpi4pyfile, os.path.join('build', mpi4pyfile))
        
        subprocess.call([
            'tar',
            '-xf',
            mpi4pyfile
            ]
            ,
            cwd=os.path.join('build')
        )
        subprocess.call([
            sys.executable,
            'setup.py',
            'install',
            '--prefix=/usr',
            '--root={0}'.format(package_path),
            ]
            ,
            cwd=os.path.join('build','mpi4py-1.2.2')
        )
        
        subprocess.call([
            sys.executable,
            'setup.py',
            'install',
            '--prefix=/usr',
            '--root={0}'.format(package_path),
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
        if os.path.exists('./'+package_name):
            print "trying to create debian package.."
            os.system('fakeroot dpkg-deb --build {0}'.format(package_name))
            
    def cleanup(self):
        package_path = os.path.abspath(package_name)
        
        if os.path.exists(package_path):
            shutil.rmtree(package_path)
            
        os.system('python setup.py generate_main')
            
    def __repr__(self):
        s = "generated debian package: {0}.deb".format(package_name)
        return s

    def touch(self, filename):
        if not os.path.exists(filename):
            open(filename, 'w').close()
            
    
if __name__ == '__main__':
    instance = generate_debian_package()
   
