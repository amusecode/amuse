import sys
import ctypes
import tempfile
import shutil
import os

import atexit


class _ModuleRegister(object):
    is_cleanup_registered = False
    files_to_cleanup = []

def find_shared_object_file(dirpath, base_libname):
    for path in sys.path:
        fullpath = os.path.join(path, dirpath)
        if os.path.exists(fullpath) and os.path.isdir(fullpath):
            full_libname = os.path.join(fullpath, base_libname)
            if os.path.exists(full_libname):
                return full_libname
    return base_libname
    
def find_module(modulename):
    parts = modulename.split('.')
    modulename = parts[-1]
    dirparts = parts[:-1]
    base_libname = modulename + '.so'
    if len(dirparts) > 0:
        dirpath = os.path.join(*dirparts)
    else:
        dirpath = ''
    libname = find_shared_object_file(dirpath, base_libname)
    if not os.path.exists(libname):
        raise Exception("cannot find the shared object file of the module '{0}'".format(modulename))
    return modulename, libname
    
def import_unique(modulename):
    modulename, libname = find_module(modulename)
    
    if modulename in sys.modules:
        prevmodule = sys.modules[modulename]
    else:
        prevmodule = None
    if not _ModuleRegister.is_cleanup_registered:
        _ModuleRegister.is_cleanup_registered = True
        atexit.register(cleanup)
        
    if not os.path.exists('__modules__'):
        os.mkdir('__modules__')
    try:
        with tempfile.NamedTemporaryFile(suffix=".so", dir='__modules__', delete=False) as target:
            with open(libname, "rb") as source:
                shutil.copyfileobj(source, target)
            target.flush()
            
            _ModuleRegister.files_to_cleanup.append(target.name)
            
            lib = ctypes.pydll.LoadLibrary(target.name)
            initfunc = getattr(lib, "init"+modulename)
            initfunc()
            return sys.modules[modulename]
    finally:
        if prevmodule is None and modulename in sys.modules:
            del sys.modules[modulename]
        else:
            sys.modules[modulename] = prevmodule
            
def cleanup():
    for filename in _ModuleRegister.files_to_cleanup:
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except Exception as ex:
                print "Could not delete file:",filename,", exception:",ex
