__revision__ = "$Id:$"

import sys, os, re, subprocess
from stat import ST_MODE
from distutils import sysconfig
from distutils.core import Command
from distutils.dep_util import newer
from distutils.util import convert_path
from distutils import log

# check if Python is called on the first line with this expression
first_line_re = re.compile('^#!.*python[0-9.]*([ \t].*)?$')

class build_latex(Command):

    description = "build latex documents (convert to pdf)"

    user_options = [
        ('build-dir', 'd', "directory to copy pdf to"),
        ('build-temp', 't', "directory to put temporary build by-products"),
        ('force', 'f', "forcibly build everything (ignore file timestamps")
        ]

    boolean_options = ['force']


    def initialize_options (self):
        self.build_dir = None
        self.build_temp = None
        self.latex_documents = None
        self.force = None
        self.outfiles = None

    def finalize_options (self):
        self.set_undefined_options('build',
                                   ('build_lib', 'build_dir'),
                                   ('build_temp', 'build_temp'),
                                   ('force', 'force'))
        self.latex_documents = []
        self.find_latex_document()

    def find_latex_document(self):
        for dir, dirs, files in os.walk('doc'):
            for x in files:
                if not dir.startswith('_') and x.endswith('.tex'):
                    self.latex_documents.append(os.path.join(dir,x))
    

    def get_source_files(self):
        return self.latex_documents

    def run (self):
        if not self.latex_documents:
            return
        self.make_pdf_files()


    def make_pdf_files (self):
        
        outfiles = []
        for x in self.get_source_files():
            x = convert_path(x)
            (dir,name) = os.path.split(x)
            pdfName = (name.rsplit('.',1))[0] + '.pdf'
            outdir = os.path.join(self.build_dir, dir)
            outfile = os.path.join(outdir, pdfName)
            tempdir = os.path.join(self.build_temp, dir)
            self.mkpath(tempdir)
            self.mkpath(outdir)
            outfiles.append(outfile)
    
            if not self.force and not newer(x, outfile):
                log.debug("not converting %s (up-to-date)", outfile)
                continue
            
            log.info("making pdf of %s", x)
            for i in range(0,2):
                process = subprocess.Popen(['pdflatex','-output-directory',tempdir, x], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
                process.communicate()
            tempfile = os.path.join(tempdir, pdfName)
            self.copy_file(tempfile,outfile)
            
            
            
            
            
    
            
