import os
import shutil

version = '2.2-1'
amuse_version = 'amuse-{0}'.format(version)
debianversion = 'ubuntu'
architecture = 'amd64'

#f=os.popen('uname -m')
#architecture = f.readlines()[0]
#print architecture
#f.close()

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
    'python-h5py (>=1.2.1)',
    'libhdf5-serial-dev (>=1.6)',
    'hdf5-tools',
    'libfftw3-3', 
    'libfftw3-dev', 
    'libfftw3-doc',
    'mpich2 (>=1.2.1)'
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

dirtree = \
[
package_name,
'./{0}/DEBIAN'.format(package_name),
'./{0}/usr'.format(package_name),
'./{0}/usr/bin'.format(package_name),
'./{0}/usr/share'.format(package_name),
'./{0}/usr/share/{1}'.format(package_name, amuse_version),
'./{0}/usr/share/doc'.format(package_name),
'./{0}/usr/share/doc/{1}'.format(package_name, amuse_version)
]

class generate_debian_package(object):
    def __init__(self):
        self.makescripts()
        self.makedirs()
        self.package()
        self.cleanup()

    def makescripts(self):
        os.system('python setup.py generate_main --amuse-dir=/usr/share/{0}'.format(amuse_version))

    def makedirs(self):
        if os.path.exists('./'+package_name):
            shutil.rmtree(package_name)

        for d in dirtree:
            os.mkdir(d)

        dependency_string = ', '.join(depends)

        f = open('./{0}/DEBIAN/control'.format(package_name),'w')
        f.write(control.format(version, dependency_string))
        f.close()

        f = open('./{0}/DEBIAN/postinst'.format(package_name),'w')
        f.write(postinst)
        f.close()
        os.chmod('./{0}/DEBIAN/postinst'.format(package_name), 0b111101101)

        shutil.copytree('src', './{0}/usr/share/{1}/src'.format(package_name, amuse_version))
        shutil.copytree('test', './{0}/usr/share/{1}/test'.format(package_name, amuse_version))
        shutil.copytree('data', './{0}/usr/share/{1}/data'.format(package_name, amuse_version))
        shutil.copytree('lib', './{0}/usr/share/{1}/lib'.format(package_name, amuse_version))
        shutil.copytree('doc', './{0}/usr/share/doc/{1}/doc'.format(package_name, amuse_version))
        self.touch('./{0}/usr/share/{1}/build.py'.format(package_name, amuse_version))

        shutil.copy('amuse.sh', './{0}/usr/bin'.format(package_name))
        shutil.copy('iamuse.sh', './{0}/usr/bin'.format(package_name))

        os.chmod('./{0}/usr/bin/amuse.sh'.format(package_name), 0b111101101)
        os.chmod('./{0}/usr/bin/iamuse.sh'.format(package_name), 0b111101101)

    def package(self):
        if os.path.exists('./'+package_name):
            print "trying to create debian package.."
            os.system('fakeroot dpkg-deb --build {0}'.format(package_name))
            
    def cleanup(self):
        shutil.rmtree(package_name)
        os.system('python setup.py generate_main')
            
    def __repr__(self):
        s = "generated debian package: {0}.deb".format(package_name)
        return s

    def touch(self, filename):
        if not os.path.exists(filename):
            open(filename, 'w').close()
            
    
if __name__ == '__main__':
    instance = generate_debian_package()
    print instance
