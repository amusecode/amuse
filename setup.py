from distutils.core import setup
from distutils.command.build import build
from distutils.command.clean import clean
from distutils.command.install import install
from distutils.util import convert_path
from distutils.cmd import Command
from distutils.extension import Extension

from support.generate_main import generate_main
from support.build_latex import build_latex
from support.setup_codes import BuildCodes, CleanCodes, DistCleanCodes, BuildOneCode
from support.setup_codes import GenerateInstallIni
from support.run_tests import run_tests

import os
import fnmatch
import re
import glob

#include_dirs.append(sysconfig.get_python_inc())

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
    'build_code':BuildOneCode,
    'clean_codes':CleanCodes,
    'dist_clean':DistCleanCodes,
    'clean_python':clean,
    'clean': Clean,
    'tests':run_tests, 
    'generate_main': generate_main,
    'generate_install_ini': GenerateInstallIni,
    'install':install
}

build.sub_commands.append(('build_codes', None))
Clean.sub_commands.append(('clean_codes',None))
Clean.sub_commands.append(('clean_python',None))

Install.sub_commands.insert(0, ('generate_install_ini',None) )

def find_packages(where='.', exclude=()):
    """Return a list all Python packages found within directory 'where'

    'where' should be supplied as a "cross-platform" (i.e. URL-style) path; it
    will be converted to the appropriate local path syntax.  'exclude' is a
    sequence of package names to exclude; '*' can be used as a wildcard in the
    names, such that 'foo.*' will exclude all subpackages of 'foo' (but not
    'foo' itself).
    """
    out = []
    stack=[(convert_path(where), '')]
    while stack:
        where,prefix = stack.pop(0)
        for name in os.listdir(where):
            fn = os.path.join(where,name)
            if ('.' not in name and os.path.isdir(fn) and
                os.path.isfile(os.path.join(fn,'__init__.py'))
            ):
                out.append(prefix+name); stack.append((fn,prefix+name+'.'))
    for pat in list(exclude)+['ez_setup', 'distribute_setup']:
        from fnmatch import fnmatchcase
        out = [item for item in out if not fnmatchcase(item,pat)]
    return out

def find_data_files(srcdir, destdir, *wildcards, **kw):
    """
    get a list of all files under the srcdir matching wildcards,
    returned in a format to be used for install_data
    """
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname:
            return
        names = []
        lst, wildcards, dirnameconverter, destdir = arg
        for wc in wildcards:
            wc_name = os.path.normpath(os.path.join(dirname, wc))
            for f in files:
                filename = os.path.normpath(os.path.join(dirname, f))

                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            destdirname = dirnameconverter.sub(destdir, dirname)
            lst.append( (destdirname, names ) )

    file_list = []
    recursive = kw.get('recursive', True)
    converter = re.compile('^({0})'.format(srcdir))
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards, converter, destdir))
    else:
        walk_helper((file_list, wildcards, converter, destdir),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(os.path.join(srcdir, '*'))])
    return file_list

all_data_files = find_data_files('data', 'share/amuse/data', '*', recursive = True)
all_data_files.extend(find_data_files('support', 'share/amuse/support', '*', recursive = False))
all_data_files.extend(find_data_files('lib', 'share/amuse/lib', '*.h', '*.a', '*.mod', '*.inc', recursive = True))
all_data_files.append(('share/amuse',['./config.mk','./build.py']))

packages = find_packages('src')
packages.extend(['amuse.test.suite.' + x for x in find_packages('test')])
packages.extend(['amuse.examples.' + x for x in find_packages('examples')])

package_data = {
    'amuse.rfi': ['*.template'],
    'amuse.test.suite.core_tests': [
        '*.txt', '*.dyn', '*.ini', 
        '*.nemo',
        '*.dat', 'gadget_snapshot'
    ],
    'amuse.test.suite.ticket_tests': [
        '*.out'
    ],
    'amuse': [
        '*rc'
    ]
}

setup(
    name = 'amuse',
    version = '6.0',
    cmdclass = mapping_from_command_name_to_command_class,
    ext_modules = extensions,
    package_dir = {'': 'src', 'amuse.test.suite' :'test', 'amuse.examples' : 'examples'},
    packages =  packages,
    package_data = package_data,
    data_files = all_data_files,
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
    url = 'http://www.amusecode.org/',
    author_email = 'info@amusecode.org',
    author = 'The Amuse Team',
)

