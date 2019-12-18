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
            result = sys.modules[modulename]
            result.__ctypeslib__ = lib
            result.__ctypesfilename__ = target.name
            result.__cleanup__ = cleanup_module
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
                print("Could not delete file:",filename,", exception:",ex)

# this struct will be passed as a ponter,
# so we don't have to worry about the right layout
class dl_phdr_info(ctypes.Structure):
  _fields_ = [
    ('padding0', ctypes.c_void_p), # ignore it
    ('dlpi_name', ctypes.c_char_p),
                            # ignore the reset
  ]


# call back function, I changed c_void_p to c_char_p
callback_t = ctypes.CFUNCTYPE(ctypes.c_int,
                       ctypes.POINTER(dl_phdr_info), 
                       ctypes.POINTER(ctypes.c_size_t), ctypes.c_char_p)

dl_iterate_phdr = ctypes.CDLL('libc.so.6').dl_iterate_phdr
# I changed c_void_p to c_char_p
dl_iterate_phdr.argtypes = [callback_t, ctypes.c_char_p]
dl_iterate_phdr.restype = ctypes.c_int


count = [0]
def callback(info, size, data):
    # simple search
    print("CLEANUP:", info.contents.dlpi_name)
    count[0] += 1
    return 0
  
def cleanup_module(mod):
    #print "CLEANUP!!"
    #sys.stdout.flush()
    #print "CLEANUP:", mod, len(list(os.listdir('/proc/self/fd')))
    #count[0] = 0
    #dl_iterate_phdr(callback_t(callback), "")
    #print "CLEANUP:", count[0]
    sys.stdout.flush()
    
    if hasattr(mod, '__ctypeslib__') and not mod.__ctypeslib__ is None:
        lib = mod.__ctypeslib__
        dlclose = ctypes.cdll.LoadLibrary('libdl.so').dlclose
        dlclose.argtypes = [ctypes.c_void_p]
        dlclose.restype = ctypes.c_int
        errorcode =  dlclose(lib._handle)
        mod.__ctypeslib__ = None
        filename = mod.__ctypesfilename__
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except Exception as ex:
                print("CLEANUP Could not delete file:",filename,", exception:",ex)
        mod.__ctypesfilename__ = None
                
