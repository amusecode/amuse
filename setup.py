from distutils.core import setup
from distutils.command.build import build
from support.build_latex import build_latex
from support.run_tests import run_tests

build.sub_commands.append(('build_latex',None))

setup(
	name = 'amuse',
	version = '1.0',
	cmdclass = {'build_latex':build_latex, 'tests':run_tests}
)
#latex_documents = ['doc/install/amuse-scope.tex','doc/instal/amuse-overview.tex'],
	
