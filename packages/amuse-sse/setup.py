import sys
import os

from distutils.command.build import build
from distutils.command.clean import clean
from distutils.command.install import install
from distutils.util import convert_path
from setuptools import setup

from support.setup_codes import (
    BuildCodes, CleanCodes, DistCleanCodes, BuildOneCode, BuildLibraries,
    ConfigureCodes, GenerateInstallIni, InstallLibraries,
)

if sys.hexversion > 0x03000000:
    from distutils.command.build_py import build_py_2to3

name = 'amuse-sse'
version = "12.0a5"
author = 'The AMUSE team'
author_email = 'info@amusecode.org'
license = "Apache License 2.0"
url = 'http://www.amusecode.org/'
install_requires = [
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'nose>=0.11.1',
    'mpi4py>=1.1.0',
    'h5py>=1.1.0',
    'amuse-framework>=12.0a5',
]
description = 'The Astrophysical Multipurpose Software Environment'
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
    'build_codes': BuildCodes,
    'build_libraries': BuildLibraries,
    'build_code': BuildOneCode,
    'configure_codes': ConfigureCodes,
    'clean_codes': CleanCodes,
    'dist_clean': DistCleanCodes,
    'clean_python': clean,
    'clean': Clean,
    'generate_install_ini': GenerateInstallIni,
    'install': install,
    'install_libraries': InstallLibraries
}

if sys.hexversion > 0x03000000:
    mapping_from_command_name_to_command_class['build_py'] = build_py_2to3

build.sub_commands.insert(0, ('configure_codes', None))
build.sub_commands.append(('build_codes', None))
Clean.sub_commands.append(('clean_codes', None))
Clean.sub_commands.append(('clean_python', None))


def find_packages(where='.', exclude=()):
    """Return a list all Python packages found within directory 'where'

    'where' should be supplied as a "cross-platform" (i.e. URL-style) path; it
    will be converted to the appropriate local path syntax.  'exclude' is a
    sequence of package names to exclude; '*' can be used as a wildcard in the
    names, such that 'foo.*' will exclude all subpackages of 'foo' (but not
    'foo' itself).
    """
    out = []
    stack = [(convert_path(where), '')]
    while stack:
        where, prefix = stack.pop(0)
        for name in os.listdir(where):
            fn = os.path.join(where, name)
            if (
                    '.' not in name and os.path.isdir(fn) and
                    os.path.isfile(os.path.join(fn, '__init__.py'))
            ):
                out.append(prefix+name)
                stack.append((fn, prefix+name+'.'))
    for pat in list(exclude)+['ez_setup', 'distribute_setup']:
        from fnmatch import fnmatchcase
        out = [item for item in out if not fnmatchcase(item, pat)]
    return out


all_data_files = []

packages = ['amuse.community.sse']

package_data = {
}

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={'amuse.community.sse': 'src/amuse/community/sse'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)
