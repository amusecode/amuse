import sys
import os

from distutils.command.build import build
from distutils.command.clean import clean
from distutils.command.install import install
from distutils.util import convert_path
from setuptools import setup, find_packages

from support.generate_main import generate_main
from support.build_latex import build_latex
from support.setup_codes import (
    BuildCodes, CleanCodes, DistCleanCodes, BuildOneCode, BuildLibraries,
    ConfigureCodes, GenerateInstallIni, InstallLibraries,
)
from support.run_tests import run_tests
from support.misc import find_data_files

if sys.hexversion > 0x03000000:
    from distutils.command.build_py import build_py_2to3

extensions = []

class Clean(clean):

    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

class Install(install):

    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

        install.run(self)

mapping_from_command_name_to_command_class = {
    'build_latex':build_latex,
    'build_codes':BuildCodes,
    'build_libraries':BuildLibraries,
    'build_code':BuildOneCode,
    'configure_codes':ConfigureCodes,
    'clean_codes':CleanCodes,
    'dist_clean':DistCleanCodes,
    'clean_python':clean,
    'clean': Clean,
    'tests':run_tests,
    'generate_main': generate_main,
    'generate_install_ini': GenerateInstallIni,
    'install':install,
    'install_libraries':InstallLibraries
}

if sys.hexversion > 0x03000000:
    mapping_from_command_name_to_command_class['build_py'] = build_py_2to3

build.sub_commands.insert(0, ('configure_codes', None))
build.sub_commands.append(('build_codes', None))
Clean.sub_commands.append(('clean_codes', None))
Clean.sub_commands.append(('clean_python', None))

Install.sub_commands.insert(0, ('generate_install_ini', None))
Install.sub_commands.append(('install_libraries', None))

all_data_files = find_data_files('data', 'share/amuse/data', '*', recursive=True)
# all_data_files.extend(find_data_files('support', 'share/amuse/support', '*', recursive=False))
# all_data_files.extend(find_data_files('support3', 'share/amuse/support3', '*', recursive = False))
# all_data_files.extend(find_data_files('lib', 'share/amuse/lib', '*.h', '*.a', '*.mod', '*.inc', '*.so', '*.dylib', recursive = True))
all_data_files.append(('share/amuse', ['./config.mk', './build.py']))

packages = find_packages('src')
packages.extend(['amuse.test.suite.' + x for x in find_packages('test')])
packages.extend(['amuse.examples.' + x for x in find_packages('examples')])

package_data = {
    'amuse.rfi.tools': ['*.template'],
    'amuse.test.suite.core_tests': [
        '*.txt', '*.dyn', '*.ini',
        '*.nemo',
        '*.dat', 'gadget_snapshot'
    ],
    'amuse.test.suite.codes_tests': [
        '*.txt', 'test_sphray_data*'
    ],
    'amuse.test.suite.ticket_tests': [
        '*.out'
    ],
    'amuse': [
        '*rc'
    ]
}

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='amuse',
    version="11",
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={'': 'src', 'amuse.test.suite' :'test', 'amuse.examples' : 'examples'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    scripts=[ "amusifier" ],
    classifiers=[
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
    ],
    url='http://www.amusecode.org/',
    author_email='info@amusecode.org',
    author='The Amuse Team',
    description='The Astrophysical Multipurpose Software Environment',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'docutils>=0.6',
        'numpy>=1.2.2',
        'nose>=0.11.1',
        'mpi4py>=1.1.0',
        'h5py>=1.1.0',
    ]
)
