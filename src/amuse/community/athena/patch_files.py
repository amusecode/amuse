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
    
def is_quilt_installed():
    if sys.platform == 'win32':
        return False
    path = which('quilt')
    if path is None:
        return False
        
    stdoutstring, stderrstring, returncode = execute_command_line(['quilt', '--version'])
    if not returncode == 0:
        return False
    
    version_re = re.compile(r'(\d).(\d\d)')
    match = version_re.match(stdoutstring)
    if not match:
        return False
    
    return True

def apply_patches_using_quilt():
    returncode = subprocess.call(['quilt', 'push', '-a'])
    if not returncode == 0:
        raise Exception("error in applying the patches, please apply by hand using quilt push")
        
def undo_patches_using_quilt():
    returncode = subprocess.call(['quilt', 'pop', '-a'])
    if not returncode == 0:
        raise Exception("error in undoing the patches, please undo by hand using quilt pop -a")

def run_patch(patchname, patchfile):
    arguments = ['patch', '-p1', '--backup', '--prefix={0}/{1}/'.format(QUILT_PC, patchname), '-E', '-i', patchfile]
    returncode = subprocess.call(arguments)
    if not returncode == 0:
        raise Exception("could not apply patch {0}".format(patchname))
        
def apply_patches_using_patch():
    with open("patches/series", "r") as f:
        lines = f.readlines()
    patches = [x.strip() for x in lines]
    patches = [x for x in patches if len(x) > 0]
    for patch in patches:
        path = os.path.join(PATCHESDIR, patch)
        run_patch(patch, path)
        
def main(undo_patches = False):
    print("checking if quilt is installed ... ")
    if not is_quilt_installed():
        print("... no")
        
        if undo_patches:
            print("quilt is not installed, cannot undo the patches")
            sys.exit(1)
        else:
            print("applying patches to source code")
            apply_patches_using_patch()
    else:
        print("... yes")
        
        if undo_patches:
            print("quilt is install, will try to undo the patches")
            undo_patches_using_quilt()
        else:
            print("applying patches to source code")
            apply_patches_using_quilt()
            print("all patches applied")
    
if __name__ == '__main__':
    main()
