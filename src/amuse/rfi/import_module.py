import sys
import ctypes
import tempfile
import shutil
import os

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
    
    
    try:
        with tempfile.NamedTemporaryFile(suffix=".so") as target:
            with open(libname, "rb") as source:
                shutil.copyfileobj(source, target)
            target.flush()
            
            lib = ctypes.pydll.LoadLibrary(target.name)
            initfunc = getattr(lib, "init"+modulename)
            initfunc()
            return sys.modules[modulename]
    finally:
        if prevmodule is None:
            del sys.modules[modulename]
        else:
            sys.modules[modulename] = prevmodule
