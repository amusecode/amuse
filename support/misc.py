import sys
import re
import os
import fnmatch

if sys.hexversion > 0x03000000:
    from os import walk as py_walk
    def walk(top, callback, args):
        for root, dirs, files in py_walk(top):
            callback(args, root, files)
else:
    from os.path import walk

def find_data_files(srcdir, destdir, *wildcards, **kw):
    """
    get a list of all files under the srcdir matching wildcards,
    returned in a format to be used for install_data
    """
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname:
            return
        names = []
        lst, wildcards, dirnameconverter, destdir = arg
        for wc in wildcards:
            wc_name = os.path.normpath(os.path.join(dirname, wc))
            for f in files:
                filename = os.path.normpath(os.path.join(dirname, f))

                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            destdirname = dirnameconverter.sub(destdir, dirname)
            lst.append((destdirname, names))

    file_list = []
    recursive = kw.get('recursive', True)
    converter = re.compile('^({0})'.format(srcdir))

    if recursive:
        walk(srcdir, walk_helper, (file_list, wildcards, converter, destdir))
    else:
        walk_helper((file_list, wildcards, converter, destdir),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(os.path.join(srcdir, '*'))])
    return file_list
