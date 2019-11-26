import os
import ctypes
import numpy
from amuse.support import exceptions
import threading

# run a little code, to create import error 
# for numpy in pypy
keys=numpy.ascontiguousarray([1,2,3], dtype="uintp")
keys.ctypes

class cell(ctypes.Structure):
    _fields_=[("key", ctypes.c_size_t),
              ("value", ctypes.c_size_t)]

cell_pointer=ctypes.POINTER(cell)
c_size_t_pointer=ctypes.POINTER(ctypes.c_size_t)
c_int_pointer=ctypes.POINTER(ctypes.c_int)

class simple_hash(ctypes.Structure):
    _fields_=[("m_cells", cell_pointer),
              ("m_arraySize", ctypes.c_size_t),
              ("m_population", ctypes.c_size_t),
              ("m_zeroUsed",ctypes.c_bool),
              ("m_zeroCell",cell)]

from amuse.support import get_amuse_root_dir
librarypath=os.path.join(get_amuse_root_dir(),"lib","simple_hash","libsimplehash.so")

lib_simple_hash=ctypes.CDLL(librarypath)

class SimpleHash(object):
    def __init__(self):
        self._map=simple_hash()
        self._dummy=ctypes.c_size_t()
        self._lib=lib_simple_hash
        self._map_ref = ctypes.byref(self._map)
        if self._lib.init_hash(ctypes.byref(self._map),128)!=0:
          raise MemoryError("allocation of SimpleHash")
        self.lock=threading.Lock()

    def __del__(self):
        self._lib.end_hash(self._map_ref)

    def lookup(self,inkeys):
        N=len(inkeys)
        keys=numpy.ascontiguousarray(inkeys, dtype="uintp")

        values=numpy.ascontiguousarray(numpy.zeros(N),dtype="uintp")
        errors=numpy.ascontiguousarray(numpy.zeros(N),dtype="int32")
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)
        cerrors=errors.ctypes.data_as(c_int_pointer)
        with self.lock:
            err=self._lib.hash_lookups(self._map_ref,N,ckeys,cvalues, cerrors)
        if err != 0:
            has_errors = errors!=0
            missing_keys = keys[has_errors]
            no_errors = ~has_errors
            raise exceptions.KeysNotInStorageException(keys[no_errors], values[no_errors], missing_keys)
        return values

    def insert(self,keys,values=None):
        N=len(keys)
        keys=numpy.ascontiguousarray(keys, dtype="uintp")
        if values is None:
            values=numpy.arange(N,dtype="uintp")
        else:
            assert len(keys)==len(values)
        values=numpy.ascontiguousarray(values,dtype="uintp")
      
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)
        with self.lock:
            err=self._lib.hash_inserts(self._map_ref,N,ckeys,cvalues)
        if err!=0:
            raise Exception("simple hash insert error")

    def reindex(self, keys, values=None):
        with self.lock:
            self._lib.end_hash(self._map_ref)
            if self._lib.init_hash(self._map_ref,len(keys))!=0:
                raise MemoryError("allocation of SimpleHash")
        self.insert(keys, values)

    def key_present(self,key):
        with self.lock:
            return self._lib.hash_lookup(self._map_ref,ctypes.c_size_t(key),ctypes.byref(self._dummy))==0


    def keys_present(self,keys):
        pass

    def match(self,keys):
        N=len(keys)
        keys=numpy.ascontiguousarray(keys, dtype="uintp")
        values=numpy.ascontiguousarray(numpy.zeros(N),dtype="uintp")
        errors=numpy.ascontiguousarray(numpy.zeros(N),dtype="int32")
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)
        cerrors=errors.ctypes.data_as(c_int_pointer)
        with self.lock:
            err=self._lib.hash_lookups(self._map_ref,N,ckeys,cvalues, cerrors)
        state = errors != 0
        return keys[state], values[state], keys[~state]

if __name__=="__main__":

  N=1000
  sm=SimpleHash()

  keys=[1,3,5,6,32,2]
  
  sm.reindex(keys)
  
  print(sm.lookup([5,32])) 
  
  del sm
