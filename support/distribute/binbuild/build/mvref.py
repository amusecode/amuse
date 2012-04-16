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


def get_bin_files(path='.'):
    for dirname, subdirs, names in os.walk(path):
        for name in names:
            fullname = os.path.join(dirname, name)
            if os.path.islink(fullname):
                continue
                    
            outputstring = subprocess.check_output(['file', fullname])
            if not outputstring.find('Mach-O') >= 0:
                continue

            if (
                not outputstring.find('shared library') >= 0 and
                not outputstring.find('executable') >= 0 and
                not outputstring.find('bundle') >= 0
               ):
                continue

            yield os.path.normpath(os.path.abspath(fullname))


def get_bin_references(name):
    outputstring = subprocess.check_output(['otool', '-L', name])
    lines = outputstring.strip().splitlines()
    for x in lines[1:]:
        x = x.strip()
        index = x.find('(comp')
        if index > 0:
            x = x[:index]
            x = x.strip()
        yield x

def get_dylib_id(name):
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

def change_dylib_ref(name, oldid, newid, dryrun = True):
    arguments = ['install_name_tool', '-change', oldid, newid, name]
    if dryrun == True:
        print ' '.join(arguments)
        return
    
    try:
        outputstring = subprocess.check_output(arguments)
    except subprocess.CalledProcessError as ex:
        print "error while running:", ' '.join(arguments)
        print "output:", ex.output
    print outputstring

def add_rpath(name, rpath,  dryrun = True):
    arguments = ['install_name_tool', '-add_rpath', rpath, name]
    if dryrun == True:
        print ' '.join(arguments)
        return
        
    print ' '.join(arguments)
    try:
        outputstring = subprocess.check_output(arguments)
        print outputstring
    except subprocess.CalledProcessError as ex:
        print "error in addind rpath:"
        print ex
        print "will ignore"

def mainold(path='../lib', binpath='.', rpath='', dryrun = True):
    
    abslibpath = os.path.abspath(path)

    if len(rpath) == 0:
        absbinpath = os.path.abspath(binpath)
        relpath = os.path.relpath(abslibpath, absbinpath)
        print relpath
        rpath = os.path.join('@loader_path',relpath)
         
    for name in get_bin_files(binpath):
        must_add_rpath = False
        for x in get_bin_references(name):
            dir, basename = os.path.split(x)
            if dir == abslibpath:
                must_add_rpath = True
                change_dylib_ref(name, x, os.path.join('@rpath', basename), dryrun)

        if must_add_rpath:
            add_rpath(name, rpath, dryrun)


    #for x in get_dylib_files(path):
    #    dylibid = get_dylib_id(x)	
    #    basename = os.path.basename(dylibid)
    #    change_dylib_id(x, os.path.join('@rpath', basename), dryrun)
        

def main(libpath='../lib', binpath='.', rpath='', dryrun = True):
    
    abslibpath = os.path.abspath(libpath)

    for name in get_bin_files(binpath):
        must_add_rpath = False
        
        print "file:", name
        for x in get_bin_references(name):
            
            dir, basename = os.path.split(x)
            print dir
            print basename
            print abslibpath
            if dir == abslibpath or dir == '@rpath':
                must_add_rpath = True
                change_dylib_ref(name, x, os.path.join('@rpath', basename), dryrun)
        
        if must_add_rpath:
            if len(rpath) == 0:
                absbinpath = os.path.abspath(name)
                relpath = os.path.relpath(abslibpath, os.path.dirname(absbinpath))
                #print "BBBB"
                #print relpath
                #print binpath
                #print abslibpath
                #print os.path.dirname(absbinpath)
                newrpath = os.path.join('@loader_path',relpath)
            else:
                newrpath = rpath
            print "newrpath", newrpath, dryrun
            add_rpath(name, newrpath, dryrun)


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-r", "--rpath", 
        default = '',
        dest="rpath",
        help="new path for dynamic libraries, default will create relative path between bin and libdir",
        type="string"
    )
    result.add_option(
        "-p", "--path", 
        default = '../lib',
        dest="libpath",
        help="path to scan for dynamic libraries",
        type="string"
    )
    result.add_option(
        "-b", "--bin-path", 
        default = '.',
        dest="binpath",
        help="path to scan for binaries referencing dynamic libraries",
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
