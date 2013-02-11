__revision__ = "$Id:$"

import sys, os, re, subprocess,stat

from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.util import convert_path
from distutils import log

try:
    from . import config
    is_configured = hasattr(config, 'compilers')
except ImportError:
    is_configured = False

class generate_main(Command):

    description = "generate shell script to run amuse"

    user_options = [
        ('amuse-dir=', 'd', "root directory of the amuse project"),
    ]

    def initialize_options (self):
        self.amuse_dir = None

    def finalize_options (self):
        if self.amuse_dir is None:
            self.amuse_dir =os.path.dirname(os.path.dirname(__file__))
            #self.amuse_dir ='/usr/share/amuse-2.2'

    def get_source_files(self):
        return self.latex_documents

    def run (self):
        test_directory = os.path.join(self.amuse_dir, 'test')
        src_directory = os.path.join(self.amuse_dir, 'src')

        with open('amuse.sh','w') as script_file:
            script_file.write('#!/bin/sh')
            script_file.write('\n\n')
            script_file.write('export PYTHONPATH=${PYTHONPATH}')
            for x in [test_directory, src_directory]:
                script_file.write(':')
                script_file.write(x)

            script_file.write('\n')
            script_file.write('export AMUSE_DIR=')
            script_file.write(self.amuse_dir)
            script_file.write('\n')
            script_file.write(sys.executable)
            script_file.write(' "$@"\n')
        os.chmod('amuse.sh', stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)

        with open('iamuse.sh','w') as script_file:
            script_file.write('#!/usr/bin/env python')
            script_file.write('\n\n')
            script_file.write('import IPython.Shell\n')
            script_file.write('import sys\n')
            for x in [test_directory, src_directory]:
                script_file.write("sys.path.append('{0}')\n".format(x))
            script_file.write('amuse_root_dir = "')
            script_file.write(self.amuse_dir)
            script_file.write('"\n')
            script_file.write('IPython.Shell.start().mainloop()\n')
        os.chmod('iamuse.sh', stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)

        if is_configured and hasattr(config, 'java') and hasattr(config.java, 'is_enabled'):
            with open('ibis-deploy.sh','w') as script_file:
                script_file.write('#!/bin/sh')
                script_file.write('\n#Deploy support process. Only works if the Ibis library has been build\n\n')

                script_file.write('export AMUSE_DIR=')
                script_file.write(self.amuse_dir)
                script_file.write('\n')

                script_file.write('export IBIS_LIB_DIR=')
                script_file.write(self.amuse_dir)
                script_file.write('/lib/ibis')
                script_file.write('\n')

                script_file.write('export JAVA=')
                script_file.write(config.java.java)
                script_file.write('\n')

                script_file.write('\n')

		script_file.write('exec ${JAVA}')
                script_file.write(' -Xmx500M')
                script_file.write(' -classpath ${IBIS_LIB_DIR}:${IBIS_LIB_DIR}/lib/*:${IBIS_LIB_DIR}/deploy/lib/*')
                script_file.write(' -Dgat.adaptor.path=${IBIS_LIB_DIR}/deploy/lib/adaptors')
                script_file.write(' -Djava.library.path=${IBIS_LIB_DIR}/deploy/lib/native_libraries')
                script_file.write(' -Dibis.deploy.home=${IBIS_LIB_DIR}/deploy')
                script_file.write(' -Damuse.home=${AMUSE_DIR}')
                script_file.write(' ibis.amuse.Daemon')
                script_file.write(' "$@"\n')

            os.chmod('ibis-deploy.sh', stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)


