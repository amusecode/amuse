import sys
import os

from setuptools import setup

import support
support.use("system")
from support.setup_codes import setup_commands

name = 'amuse-distributed'
major_version = 12
minor_version = 0
micro_version = 0
patch_version = ""  # e.g. ".rc1"
post_version = ""  # e.g. ".post1"
version = "%i.%i.%i%s%s" % (
    major_version, minor_version, micro_version, patch_version,
    post_version
)
main_version = "%i.%i.0" % (major_version, minor_version)
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'nose>=0.11.1',
    'mpi4py>=1.1.0',
    'h5py>=1.1.0',
    'amuse-framework>=%s' % (main_version),
]
description = 'The Astrophysical Multipurpose Software Environment - Distributed AMUSE'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"
classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: Apache Software License',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: C',
    'Programming Language :: C++',
    'Programming Language :: Fortran',
    'Topic :: Scientific/Engineering :: Astronomy',
]

extensions = []

all_data_files = []

packages = ['amuse.community.distributed']

package_data = {
}

mapping_from_command_name_to_command_class=setup_commands()

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={'amuse.community.distributed': 'src/amuse/community/distributed'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)
