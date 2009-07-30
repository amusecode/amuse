import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system

class SSE_new(object): 
    def __init__(self):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker1')
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
        self.channel = core.MpiChannel(self.intercomm)
    
    def __del__(self):
        self.stop_worker()
        
    @core.legacy_function
    def stop_worker():
        function = RemoteFunction()  
        function.id = 0
        return function;
        
    @core.legacy_function   
    def initialize():
        function = RemoteFunction()  
        function.id = 1
        function.addParameter('z_in', dtype='d', direction=function.IN)
        function.addParameter('neta_in', dtype='d', direction=function.IN)
        function.addParameter('bwind_in', dtype='d', direction=function.IN)
        function.addParameter('hewind_in', dtype='d', direction=function.IN)
        function.addParameter('sigma_in', dtype='d', direction=function.IN)
        function.addParameter('ifflag_in', dtype='i', direction=function.IN)
        function.addParameter('wdflag_in', dtype='i', direction=function.IN)
        function.addParameter('bhflag_in', dtype='i', direction=function.IN)
        function.addParameter('nsflag_in', dtype='i', direction=function.IN)
        function.addParameter('mxns_in', dtype='i', direction=function.IN)
        function.addParameter('pts1_in', dtype='d', direction=function.IN)
        function.addParameter('pts2_in', dtype='d', direction=function.IN)
        function.addParameter('pts3_in', dtype='d', direction=function.IN)
        function.addParameter('status', dtype='i', direction=function.OUT)
        return function
        
    @core.legacy_function     
    def evolve():
        function = RemoteFunction()  
        function.id = 2
        function.name = 'evolve0'
        function.addParameter('kw', dtype='d', direction=function.INOUT)
        function.addParameter('mass', dtype='d', direction=function.INOUT)
        function.addParameter('mt', dtype='d', direction=function.INOUT)
        function.addParameter('r', dtype='d', direction=function.INOUT)
        function.addParameter('lum', dtype='d', direction=function.INOUT)
        function.addParameter('mc', dtype='d', direction=function.INOUT)
        function.addParameter('rc,', dtype='d', direction=function.INOUT)
        function.addParameter('menv', dtype='d', direction=function.INOUT)
        function.addParameter('renv', dtype='d', direction=function.INOUT)
        function.addParameter('ospin', dtype='d', direction=function.INOUT)
        function.addParameter('epoch', dtype='d', direction=function.INOUT)
        function.addParameter('tm', dtype='d', direction=function.INOUT)
        function.addParameter('tphys', dtype='d', direction=function.INOUT)
        function.addParameter('tphysf', dtype='d', direction=function.INOUT)
        return function
        
    @core.legacy_function      
    def get_time_step():
        function = RemoteFunction()      
        function.id = 3
        function.addParameter('kw', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age', dtype='d', direction=function.IN)
        function.addParameter('mt', dtype='d', direction=function.IN)
        function.addParameter('tm', dtype='d', direction=function.IN)
        function.addParameter('epoch', dtype='d', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.OUT)
        return function
        
class SSE(object): 
    class double_property(object):
        def __init__(self, tag):
            self.tag = tag
        def __get__(self, instance, owner):
            if instance is None:
                return self
            id, int_result, double_result = instance.do_call(self.tag)
            return double_result[0]
        def __set__(self, instance, value):
            if instance is None:
                return self
            instance.do_call(self.tag, doubles_arg=[value])
        
    class int_property(object):
        def __init__(self, tag):
            self.tag = tag
        def __get__(self, instance, owner):
            if instance is None:
                return self
            id, int_result, double_result = instance.do_call(self.tag, int_arg2 = 0)
            return int_result
        def __set__(self, instance, value):
            if instance is None:
                return self
            instance.do_call(self.tag, int_arg1 = value, int_arg2 = 1)
             
    def __init__(self):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker')
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
    def __del__(self):
        self.do_call(0)

    def do_call(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_arg=[], ints_arg=[]):
        self.send_request(tag,id, int_arg1, int_arg2, doubles_arg, ints_arg)
        return self.recieve_result()
    def recieve_result(self):
        header = numpy.empty(6,  dtype='i')
        self.intercomm.Recv([header, MPI.INT], source=0, tag=999)
        id = header[1]
        int_result = header[2]
        n_doubles = header[4]
        n_ints = header[5]
        if n_doubles > 0:
            doubles_result = numpy.empty(n_doubles,  dtype='d')
            self.intercomm.Recv([doubles_result, MPI.DOUBLE], source=0, tag=999)
        else:
            doubles_result = []
        if n_ints > 0:
            ints_result = numpy.empty(n_ints,  dtype='i')
            self.intercomm.Recv([ints_result, MPI.INT], source=0, tag=999)
        else:
            ints_result = []
        if header[0] < 0:
            raise Exception("Not a valid message!")
        return (id, int_result, doubles_result, ints_result)
    def send_request(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_arg=[], ints_arg=[]):
        header = numpy.array([tag, id, int_arg1, int_arg2, len(doubles_arg), len(ints_arg)], dtype='i')
        self.intercomm.Send([header, MPI.INT], dest=0, tag=0)
        if doubles_arg:
            doubles = numpy.array(doubles_arg, dtype='d')
            self.intercomm.Send([doubles, MPI.DOUBLE], dest=0, tag=0)
        if ints_arg:
            ints = numpy.array(ints_arg, dtype='i')
            self.intercomm.Send([ints, MPI.INT], dest=0, tag=0)
    
       
    def initialize(self, z_in,  neta_in,  bwind_in,
                     hewind_in,  sigma_in,  ifflag_in,
                     wdflag_in,  bhflag_in,  nsflag_in,  mxns_in,
                     pts1_in,  pts2_in,  pts3_in ):
                         
        doubles = [z_in,  neta_in,  bwind_in, hewind_in,  sigma_in,  pts1_in,  pts2_in,  pts3_in ]
        ints = [ifflag_in, wdflag_in, bhflag_in, nsflag_in, mxns_in]
        id, int_result, double_result, ints_result = self.do_call(1, doubles_arg=doubles, ints_arg=ints)
        return int_result
    def evolve(self, kw, mass,  mt,  r,  lum, mc,  rc,  menv,  renv, ospin, epoch, tm, tphys, tphysf):
        	     
        doubles = [mass,  mt,  r,  lum, mc,  rc,  menv,  renv, ospin, epoch, tm, tphys, tphysf]
        ints = [kw]
        id, int_result, double_result, ints_result = self.do_call(3, doubles_arg=doubles, ints_arg=ints)
        result = [ints_result[0]]
        result.extend(double_result)
        return result
    
    
    def evolve2(self):
        function = RemoteFunction()      
        function.name = "evolve"
        function.addParameter('kw', dtype='d')
        function.addParameter('mt', dtype='d')
        function.addParameter('r', dtype='d')
        function.addParameter('lum', dtype='d')
        function.addParameter('mc', dtype='d')
        function.addParameter('rc,', dtype='d')
        function.addParameter('menv', dtype='d')
        function.addParameter('renv', dtype='d')
        function.addParameter('ospin', dtype='d')
        function.addParameter('epoch', dtype='d')
        function.addParameter('tm', dtype='d')
        function.addParameter('tphys', dtype='d')
        function.addParameter('tphysf', dtype='d')
        return function
    
    def get_time_step2(self):
        function = RemoteFunction()      
        function.name = "get_time_step"
        function.addParameter('kw', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age', dtype='d', direction=function.IN)
        function.addParameter('mt', dtype='d', direction=function.IN)
        function.addParameter('tm', dtype='d', direction=function.IN)
        function.addParameter('epoch', dtype='d', direction=function.OUT)
        return function
    def get_time_step(self, kw, mass, age, mt, tm, epoch):
        	     
        doubles = [mass, age, mt, tm, epoch]
        ints = [kw]
        id, int_result, double_result, ints_result = self.do_call(2, doubles_arg=doubles, ints_arg=ints)
        return double_result[0]
    
        
    
  
