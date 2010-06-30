#make debian binary package...
import os
import shutil


class generate_debian_package(object):
    pass

control = \
"""
Package: amuse
Version: 2.2
Section: user
Priority: optional
Architecture: i386
Depends: python2.6
Maintainer: amuse team <info@amusecode.org>
Description: AMUSE
"""

dirtree = \
[
'amusedeb',
'./amusedeb/DEBIAN',
'./amusedeb/usr',
'./amusedeb/usr/share',
'./amusedeb/usr/share/amuse-2.2',
'./amusedeb/usr/share/doc',
'./amusedeb/usr/share/doc/amuse-2.2'
]

def makedirs():
    if os.path.exists('./amusedeb'):
        #os.system('rm -r ./amusedeb')
        shutil.rmtree('amusedeb')
    
    for d in dirtree:
        os.mkdir(d)

    f = open('./amusedeb/DEBIAN/control','w')
    f.write(control)
    f.close()

    shutil.copytree('src', './amusedeb/usr/share/amuse-2.2/src')
    shutil.copytree('test', './amusedeb/usr/share/amuse-2.2/test')

    shutil.copytree('doc', './amusedeb/usr/share/doc/amuse-2.2/doc')
    shutil.copy('amuse.sh', './amusedeb/usr/share/amuse-2.2/src')

def package():
    if os.path.exists('./amusedeb'):
        os.system('fakeroot dpkg-deb amusedeb')
        shutil.copy('amusedeb.deb','./dist')

if __name__ == '__main__':
    makedirs()
    package()
