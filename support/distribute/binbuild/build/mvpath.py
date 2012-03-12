import os
import re
import os.path
import subprocess
from optparse import OptionParser

shared_library_re=re.compile(r'.*\.dylib.*')
def get_dylib_files(path='.'):
    for name in os.listdir(path):
        if not shared_library_re.match(name):
            continue

        fullname = os.path.join(path, name)
        if os.path.islink(fullname):
            continue
        
        yield os.path.normpath(os.path.abspath(fullname))

def get_dylib_id(name):
    print name
    outputstring = subprocess.check_output(['otool', '-D', name])
    lines = outputstring.strip().splitlines()
    return lines[1].strip()

def change_dylib_id(name, newid, dryrun = True):
    arguments = ['install_name_tool', '-id', newid, name]
    if dryrun == True:
        print ' '.join(arguments)
        return

    outputstring = subprocess.check_output(arguments)
    print outputstring

def main(path='.', rpath='@rpath', dryrun = True):
    for x in get_dylib_files(path):
        dylibid = get_dylib_id(x)	
        basename = os.path.basename(dylibid)
        change_dylib_id(x, os.path.join(rpath, basename), dryrun)


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-r", "--rpath", 
        default = '@rpath',
        dest="rpath",
        help="new path for dynamic libraries, defaults to @rpath ",
        type="string"
    )
    result.add_option(
        "-p", "--path", 
        default = '.',
        dest="path",
        help="path to scan for dynamic libraries",
        type="string"
    )
    result.add_option(
        "--dry-run", 
        default = False,
        dest="dryrun",
        help="if given will show the commands, not execute these",
        action="store_true"
    )

    return result
        
if __name__ == '__main__':
    options,argumenst = new_option_parser().parse_args()
    main(**options.__dict__)
