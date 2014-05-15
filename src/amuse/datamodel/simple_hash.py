import os
import ctypes
import numpy
from amuse.support import exceptions

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
        if self._lib.init_hash(ctypes.byref(self._map),128)!=0:
          raise MemoryError("allocation of SimpleHash")
        
    def __del__(self):
        self._lib.end_hash(ctypes.byref(self._map))

    def lookup(self,keys):
        N=len(keys)
        keys=numpy.ascontiguousarray(keys, dtype="uintp")
        values=numpy.ascontiguousarray(numpy.zeros(N),dtype="uintp")
        errors=numpy.ascontiguousarray(numpy.zeros(N),dtype="int")
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)
        cerrors=errors.ctypes.data_as(c_int_pointer)
        err=self._lib.hash_lookups(ctypes.byref(self._map),N,ckeys,cvalues, cerrors)
        if err!=0:
            missing_keys = keys[errors != 0]
            raise exceptions.KeysNotInStorageException(keys[errors==0], values[errors==0], missing_keys)
        return values

    def insert(self,keys,values=None):
        N=len(keys)
        keys=numpy.ascontiguousarray(keys, dtype="uintp")
        if values is None:
            values=numpy.arange(N)
        else:
            assert len(keys)==len(values)
        values=numpy.ascontiguousarray(values,dtype="uintp")
      
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)

        self._lib.hash_inserts(ctypes.byref(self._map),N,ckeys,cvalues)

    def reindex(self, keys,values=None):

        self._lib.end_hash(ctypes.byref(self._map))
        if self._lib.init_hash(ctypes.byref(self._map),len(keys))!=0:
          raise MemoryError("allocation of SimpleHash")
        self.insert(keys,values)

    def key_present(self,key):
        return self._lib.hash_lookup(ctypes.byref(self._map),ctypes.c_size_t(key),ctypes.byref(self._dummy))==0

    def keys_present(self,keys):
        pass

if __name__=="__main__":

  N=1000
  sm=SimpleHash()

  keys=[1,3,5,6,32,2]
  
  sm.reindex(keys)
  
  print sm.lookup([5,32]) 
  
  del sm
