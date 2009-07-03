from distutils.core import setup
from distutils.command.build import build
from distutils.extension import Extension

from support.build_latex import build_latex
from support.run_tests import run_tests
from support.config import config

#include_dirs.append(sysconfig.get_python_inc())

extensions = []

mapping_from_command_name_to_command_class = {
    'build_latex':build_latex, 
    'tests':run_tests, 
    'config': config ,
}
   
build.sub_commands.append(('build_latex',None))
 
setup(
    name = 'amuse',
    version = '1.0',
    cmdclass = mapping_from_command_name_to_command_class,
    ext_modules = extensions,
)
	
