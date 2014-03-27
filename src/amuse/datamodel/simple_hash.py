import os
import ctypes
import numpy

class cell(ctypes.Structure):
    _fields_=[("key", ctypes.c_size_t),
              ("value", ctypes.c_size_t)]

cell_pointer=ctypes.POINTER(cell)
c_size_t_pointer=ctypes.POINTER(ctypes.c_size_t)

class simple_hash(ctypes.Structure):
    _fields_=[("m_cells", cell_pointer),
              ("m_arraySize", ctypes.c_size_t),
              ("m_population", ctypes.c_size_t),
              ("m_zeroUsed",ctypes.c_bool),
              ("m_zeroCell",cell)]

librarypath=os.path.join(os.getenv("AMUSE_DIR"),"lib","simple_hash","libsimplehash.so")

lib_simple_hash=ctypes.CDLL(librarypath)

class SimpleHash(object):
    def __init__(self):
        self._map=simple_hash()
        self._dummy=ctypes.c_size_t()
        self._lib=lib_simple_hash
        self._lib.init_hash(ctypes.byref(self._map),128)
        
    def __del__(self):
        self._lib.end_hash(ctypes.byref(self._map))

    def lookup(self,keys):
        N=len(keys)
        keys=numpy.array(keys, dtype="uintp")
        values=numpy.zeros(N,dtype="uintp")
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)
        err=self._lib.hash_lookups(ctypes.byref(self._map),N,ckeys,cvalues)
        if err!=0:
          raise Exception("map lookup error")
        return values

    def insert(self,keys,values=None):
        N=len(keys)
        keys=numpy.array(keys, dtype="uintp")
        if not values:
            values=numpy.arange(N,dtype="uintp")
        else:
            assert len(keys)==len(values)
      
        ckeys=keys.ctypes.data_as(c_size_t_pointer)
        cvalues=values.ctypes.data_as(c_size_t_pointer)

        self._lib.hash_inserts(ctypes.byref(self._map),N,ckeys,cvalues)

    def reindex(self, keys):

        self._lib.end_hash(ctypes.byref(self._map))
        self._lib.init_hash(ctypes.byref(self._map),len(keys))
        self.insert(keys)

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
