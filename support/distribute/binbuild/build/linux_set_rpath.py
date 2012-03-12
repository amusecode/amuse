import os
import re
import os.path
import subprocess
from optparse import OptionParser
 
def check_output(*popenargs, **kwargs):
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        raise subprocess.CalledProcessError(retcode, cmd, output=output)
    return output

shared_library_re=re.compile(r'.*\.so.*')
def get_so_files(path='.'):
    for name in os.listdir(path):
        if not shared_library_re.match(name):
            continue

        fullname = os.path.join(path, name)
        if os.path.islink(fullname):
            continue
        
        yield fullname

def get_bin_files(path='.'):
    for dirname, subdirs, names in os.walk(path):
        for name in names:
            fullname = os.path.join(dirname, name)
            if os.path.islink(fullname):
                continue
                    
            outputstring = check_output(['file', fullname])
            if not outputstring.find('ELF') >= 0:
                continue

            if (
                not outputstring.find('shared object') >= 0 and
                not outputstring.find('executable') >= 0
               ):
                continue

            yield os.path.normpath(os.path.abspath(fullname))

def get_rpaths(name):
    arguments = ['patchelf', '--print-rpath', name]
    outputstring = check_output(arguments)
    outputstring = outputstring.strip()
    return [x for x in outputstring.split(':') if len(x) > 0]

def add_rpaths(name, rpaths,  dryrun = True):
    existing_paths = get_rpaths(name)
    
    parts = list(rpaths)
    parts.extend(existing_paths)
    
    set_rpaths(name, parts, dryrun)
    
def set_rpaths(name, rpaths, dryrun = True):
    rpath = ':'.join(rpaths)
        
    arguments = ['patchelf', '--set-rpath', rpath, name]
    if dryrun == True:
        print ' '.join(arguments)
    else:
        outputstring = check_output(arguments)
        print outputstring

def main(libpath='../lib', binpath='.', rpath='', dryrun = True):
    
    abslibpath = os.path.abspath(libpath)

    for name in get_bin_files(binpath):
        must_add_rpath = False
        
        if len(rpath) == 0:
            absbinpath = os.path.abspath(name)
            relpath = os.path.relpath(abslibpath, os.path.dirname(absbinpath))
            print relpath, abslibpath
            if not relpath or relpath == '.':
                newrpath = '$ORIGIN'
            else:
                newrpath = os.path.join('$ORIGIN',relpath)
        else:
            newrpath = rpath
            
        print "file:", name
        if dryrun:
            currentpaths = get_rpaths(name)
            if len(currentpaths) > 0:
                print "CURRENT RPATHS(s):"
            for x in get_rpaths(name):
                print " ** ",x
        
        if newrpath == '$ORIGIN' and len(get_rpaths(name)) == 0:
            print " NOT SETTING RPATH FOR: ",name
            continue
        #, '$ORIGIN'
        set_rpaths(name, [newrpath], dryrun)
    
    #for x in get_dylib_files(path):
    #    basename = os.path.basename(dylibid)
    #    change_dylib_id(x, os.path.join('@rpath', basename), dryrun)


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
        help="path containing the shared libraries",
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
