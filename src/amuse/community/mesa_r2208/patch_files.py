#!/usr/bin/env python

import os
import sys
import subprocess
import re

PATCHESDIR = "patches"
QUILT_PC = ".pc"


def execute_command_line(arguments, cwd=None):
    with subprocess.Popen(
        arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd
    ) as process:
        stdoutstring, stderrstring = process.communicate()
        returncode = process.poll()
    return stdoutstring, stderrstring, returncode


def which(executablename):
    stdoutstring, stderrstring, returncode = execute_command_line(
        ['which', executablename])
    if returncode != 0:
        return None
    return stdoutstring


def is_quilt_installed():
    path = which('quilt')
    if path is None:
        return False

    stdoutstring, stderrstring, returncode = execute_command_line(
        ['quilt', '--version'])
    if returncode != 0:
        return False

    version_re = re.compile(r'(\d).(\d\d)')
    match = version_re.match(stdoutstring)
    if not match:
        return False

    return True


def apply_patches_using_quilt():
    returncode = subprocess.call(['quilt', 'push', '-a'])
    if returncode != 0:
        raise Exception(
            "error in applying the patches, please apply by hand using "
            "quilt push"
        )


def undo_patches_using_quilt():
    returncode = subprocess.call(['quilt', 'pop', '-a'])
    if returncode != 0:
        raise Exception(
            "error in undoing the patches, please undo by hand using "
            "quilt pop -a"
        )


def run_patch(patchname, patchfile):
    arguments = [
        'patch', '-p1', f'--prefix={QUILT_PC}/{patchname}/', '-E',
        '-i', patchfile
    ]
    print(arguments)
    returncode = subprocess.call(arguments)
    if returncode != 0:
        raise Exception(f"could not apply patch {patchname}")


def apply_patches_using_patch():
    with open("patches/series", "r", encoding="utf-8") as file:
        lines = file.readlines()
    patches = [x.strip() for x in lines]
    patches = [x for x in patches if len(x) > 0]
    for patch in patches:
        path = os.path.join(PATCHESDIR, patch)
        run_patch(patch, path)


def main(undo_patches=False):
    print("checking if quilt is installed ... ", end=' ')
    if not is_quilt_installed():
        print("no")

        if undo_patches:
            print("quilt is not installed, cannot undo the patches")
            sys.exit(1)
        else:
            print("applying patches to source code")
            apply_patches_using_patch()
    else:
        print("yes")

        if undo_patches:
            print("quilt is install, will try to undo the patches")
            undo_patches_using_quilt()
        else:
            print("applying patches to source code")
            apply_patches_using_quilt()
            print("all patches applied")


if __name__ == '__main__':
    main()
