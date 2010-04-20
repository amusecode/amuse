from distutils.core import setup
from distutils.command.build import build
from distutils.command.clean import clean
from distutils.cmd import Command
from distutils.extension import Extension

from support.generate_main import generate_main
from support.build_latex import build_latex
from support.setup_legacy import BuildLegacy, CleanLegacy, BuildOneLegacyCode
from support.run_tests import run_tests
from support.config import config

#include_dirs.append(sysconfig.get_python_inc())

extensions = []

class Clean(clean):
    
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)
            
mapping_from_command_name_to_command_class = {
    'build_latex':build_latex, 
    'build_legacy':BuildLegacy,
    'code':BuildOneLegacyCode,
    'clean_legacy':CleanLegacy,
    'clean_python':clean,
    'clean': Clean,
    'tests':run_tests, 
    'config': config ,
    'generate_main': generate_main,
}
   
build.sub_commands.append(('build_legacy',None))
build.sub_commands.append(('generate_main',None))
Clean.sub_commands.append(('clean_legacy',None))
Clean.sub_commands.append(('clean_python',None))
 
setup(
    name = 'amuse',
    version = '2.2',
    cmdclass = mapping_from_command_name_to_command_class,
    ext_modules = extensions,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: C',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    url = 'http://www.amusecode.org/',
    author_email = 'info@amusecode.org',
    author = 'The Amuse Team',
)

