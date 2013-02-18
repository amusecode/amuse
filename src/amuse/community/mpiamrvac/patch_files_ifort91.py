#!/usr/bin/env python

import os
import sys
import subprocess
import re

PATCHESDIR = "patches"
QUILT_PC   = ".pc"

def execute_command_line(arguments, cwd = None):
    process = subprocess.Popen(arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd = cwd)
    stdoutstring, stderrstring = process.communicate()
    returncode = process.poll()
    return stdoutstring, stderrstring, returncode
    
def which(executablename):
    stdoutstring, stderrstring, returncode = execute_command_line(['which', executablename])
    if not returncode == 0:
        return None
    else:
        return stdoutstring
    
def run_patch(patchname, patchfile):
    arguments = ['patch', '-p1', '--backup', '--prefix={0}/{1}/'.format(QUILT_PC, patchname), '-E', '-i', patchfile]
    returncode = subprocess.call(arguments)
    if not returncode == 0:
        raise Exception("could not apply patch {0}".format(patchname))
        
def apply_patches_using_patch():
    run_patch('ifort91', 'ifort91.diff')
        
def main(undo_patches = False):
    print "applying ifort91 patches to source code"
    apply_patches_using_patch()
    
if __name__ == '__main__':
    main()
