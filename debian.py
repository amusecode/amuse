#make debian binary package...
"""
These are the packages AMUSE needs:

    * Python (version >= 2.6)
    * Numpy (version >= 1.3.0)
    * HDF (version 1.6.5 - 1.8.3)
    * h5py (version >= 1.2.0)
    * MPI (OpenMPI or MPICH2)
    * mpi4py (version >= 1.0)
    * nose (version >= 0.11)
    * FFTW (version >= 3.0)
"""

import os
import shutil

class generate_debian_package(object):
    def __init__(self):
        pass

#f=os.popen('uname -m')
#architecture = f.readlines()[0]
#print architecture
#f.close()
version = '2.2-1'
debianversion = 'ubuntu'
architecture = 'amd64'
package_name = 'amuse_{0}-{1}-{2}'.format(version, debianversion, architecture)
depends = \
[
    'bash (>= 2.05a-11)',
    'python2.6 (>=2.6.0)',
    'python-numpy (>=1.3.0)',
    'python-nose (>=0.11)',
    'python-matplotlib (>=0.99)',
    'python-setuptools',
    'python-docutils',
    'libhdf5-serial-dev (>=1.6)',
    'hdf5-tools',
    'libfftw3-3', 
    'libfftw3-dev', 
    'libfftw3-doc',
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
dirtree = \
[
package_name,
'./{0}/DEBIAN'.format(package_name),
'./{0}/usr'.format(package_name),
'./{0}/usr/bin'.format(package_name),
'./{0}/usr/share'.format(package_name),
'./{0}/usr/share/amuse-2.2'.format(package_name),
'./{0}/usr/share/doc'.format(package_name),
'./{0}/usr/share/doc/amuse-2.2'.format(package_name)
]

def makeamuse():
    os.system('make clean')
    os.system('make')

def makedirs():
    if os.path.exists('./'+package_name):
        #os.system('rm -r ./amusedeb')
        shutil.rmtree(package_name)
    
    for d in dirtree:
        os.mkdir(d)

    dependency_string = ', '.join(depends)

    f = open('./{0}/DEBIAN/control'.format(package_name),'w')
    f.write(control.format(version, dependency_string))
    f.close()

    shutil.copytree('src', './{0}/usr/share/amuse-2.2/src'.format(package_name))
    shutil.copytree('test', './{0}/usr/share/amuse-2.2/test'.format(package_name))

    shutil.copytree('doc', './{0}/usr/share/doc/amuse-2.2/doc'.format(package_name))
    shutil.copy('amuse.sh', './{0}/usr/bin'.format(package_name))
    shutil.copy('iamuse.sh', './{0}/usr/bin'.format(package_name))
    os.chmod('./{0}/usr/bin/amuse.sh'.format(package_name), 0b111101101)
    os.chmod('./{0}/usr/bin/iamuse.sh'.format(package_name), 0b111101101)

def package():
    if os.path.exists('./'+package_name):
        print "trying to create debian package.."
        os.system('fakeroot dpkg-deb --build {0}'.format(package_name))
        #shutil.copy('{0}.deb'.format(package_name),'./dist')
    
if __name__ == '__main__':
    #makeamuse()
    makedirs()
    package()
