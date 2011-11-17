from distribute_setup import use_setuptools
use_setuptools()

#from distutils.core import setup
from setuptools import setup, find_packages
from distutils.command.build import build
from distutils.command.clean import clean
from setuptools.command.install_lib import install_lib
from distutils.cmd import Command
from distutils.extension import Extension

from support.generate_main import generate_main
from support.build_latex import build_latex
from support.setup_codes import BuildCodes, CleanCodes, DistCleanCodes, BuildOneCode
from support.setup_codes import GenerateInstallIni
from support.run_tests import run_tests

import os
import fnmatch

#include_dirs.append(sysconfig.get_python_inc())

extensions = []

class Clean(clean):
    
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)
            
class InstallLib(install_lib):
    
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)
            
        install_lib.run(self)
        
            
            
mapping_from_command_name_to_command_class = {
    'build_latex':build_latex, 
    'build_codes':BuildCodes,
    'code':BuildOneCode,
    'clean_codes':CleanCodes,
    'dist_clean':DistCleanCodes,
    'clean_python':clean,
    'clean': Clean,
    'tests':run_tests, 
    'generate_main': generate_main,
    'generate_install_ini': GenerateInstallIni,
    'install_lib':InstallLib
}


Clean.sub_commands.append(('clean_codes',None))
Clean.sub_commands.append(('clean_python',None))

InstallLib.sub_commands.append( ('generate_install_ini',None) )

def find_data_files(srcdir, *wildcards, **kw):
    """
    get a list of all files under the srcdir matching wildcards,
    returned in a format to be used for install_data
    """
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname:
            return
        names = []
        lst, wildcards = arg
        for wc in wildcards:
            wc_name = os.path.normpath(os.path.join(dirname, wc))
            for f in files:
                filename = os.path.normpath(os.path.join(dirname, f))

                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            lst.append( (dirname, names ) )

    file_list = []
    recursive = kw.get('recursive', True)
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards))
    else:
        walk_helper((file_list, wildcards),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
    return file_list

setup(
    name = 'amuse',
    version = '5.1',
    cmdclass = mapping_from_command_name_to_command_class,
    ext_modules = extensions,
    package_dir = {'': 'src'},
    packages =  find_packages('src'),
    data_files= find_data_files('data', '*', recursive = True),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: GNU General Public License (GPL)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: C',
        'Programming Language :: C++',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    include_package_data = True,
    url = 'http://www.amusecode.org/',
    author_email = 'info@amusecode.org',
    author = 'The Amuse Team',
    zip_safe = False
)

