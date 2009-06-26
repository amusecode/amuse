from distutils.core import setup
from distutils.command.build import build
from support.build_latex import build_latex
from support.run_tests import run_tests

from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy.distutils.core import Extension as NumpyExtension
from numpy.distutils.core import setup as NumpySetup

from numpy.distutils.misc_util import get_numpy_include_dirs

include_dirs = []
include_dirs.extend(get_numpy_include_dirs())
include_dirs.append('src/amuse/codes/gravity/nbody1h/code')
#include_dirs.append(sysconfig.get_python_inc())

import glob

build.sub_commands.append(('build_latex',None))

bhtc_cpp_files = glob.glob('src/amuse/codes/gravity/bhtree/*/*.C')
bhtc_c_files = glob.glob('src/amuse/codes/gravity/bhtree/*/*.c')
bhtc_pyx_files = glob.glob('src/amuse/codes/gravity/bhtree/*.pyx')
files = []
files.extend(bhtc_pyx_files)
files.extend(bhtc_cpp_files)
files.extend(bhtc_c_files)
print files

nbody1_f_files = glob.glob('src/amuse/codes/gravity/nbody1h/code/*.f')
nbody1_F_files = glob.glob('src/amuse/codes/gravity/nbody1h/code/*.F')
nbody1_interface_files = glob.glob('src/amuse/codes/gravity/nbody1h/*.f')
files2 = []
#files2.extend(nbody1_interface_files)
files2.extend(nbody1_f_files)
files2.extend(nbody1_F_files)
extensions = [ #Extension('src.amuse.codes.gravity.bhtree.interface', files, include_dirs=include_dirs),
Extension('src.amuse.codes.gravity.nbody1h.lib1', nbody1_interface_files, include_dirs=include_dirs, extra_objects=['libnbody1h-1.a'], language="f90")]


NumpySetup(
	name = 'amuse',
	version = '1.0',
	cmdclass = {
	'build_latex':build_latex, 
	'tests':run_tests, 
	#'build_ext': build_ext,
	},
	
	libraries = [('src/amuse/codes/gravity/nbodyh1/nbody1h-1',  {'sources':files2, 'language':'f90'})],
	ext_modules = extensions,
)
#latex_documents = ['doc/install/amuse-scope.tex','doc/instal/amuse-overview.tex'],
	
