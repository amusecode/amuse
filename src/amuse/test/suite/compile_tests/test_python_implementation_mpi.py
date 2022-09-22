from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI

from amuse import datamodel
from amuse.rfi.core import PythonCodeInterface, legacy_function,LegacyFunctionSpecification

from amuse.support import exceptions

import numpy
from mpi4py import MPI

class ForTestingInterface(PythonCodeInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, implementation_factory = ForTestingImplementation, **options)
        
    @legacy_function
    def get_range():
        function = LegacyFunctionSpecification()  
        function.addParameter('imin', dtype='int32', direction=function.OUT)
        function.addParameter('imax', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_x():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function    

    @legacy_function
    def get_dens():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('dens', dtype='float64', direction=function.OUT)
        function.addParameter('N', dtype='int32', direction=function.LENGTH)

        function.result_type = 'int32'
        function.must_handle_array = True
        return function    

class ForTestingImplementation(object):
    
    def __init__(self):
        self.comm=MPI.COMM_WORLD
        self.myrank=self.comm.Get_rank()
        self.N=self.comm.Get_size()
        self.Ngrid=3*4*5
        n=self.Ngrid//self.N
        x = (numpy.arange(n)+self.myrank*n)/(1.*self.Ngrid)
        self.local_imin=self.myrank*n
        self.local_imax=(self.myrank+1)*n-1
        self.dens = x**2
        
    def get_range(self,imin,imax):
        imin.value=0
        imax.value=self.Ngrid-1
        return 0
        
    def get_x(self,index,x):
        x.value=index/(1.*self.Ngrid)
        return 0

    def get_dens(self,index,dens,N):
        a=(index>=self.local_imin)*(index<=self.local_imax)
        _dens=numpy.zeros(N)
        _dens[a]=self.dens[index[a]-self.local_imin]
        dens.value=numpy.zeros(N)
        _dens=self.comm.Reduce(_dens, dens.value, MPI.SUM,root=0)
        return 0

class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(**options), **options)

    def define_grids(self, object):        
        object.define_grid('grid',axes_names = ['x'],grid_class=datamodel.CartesianGrid)
        object.set_grid_range('grid', 'get_range')
        object.add_getter('grid', 'get_dens', names=('dens',))
        object.add_getter('grid', 'get_x', names=('x',))
            
class TestInterface(TestWithMPI):
    
    def ForTesting(self, **options):
        options["worker_dir"]=self.get_path_to_results()
        return ForTesting( **options)
    def ForTestingInterface(self, **options):
        options["worker_dir"]=self.get_path_to_results()
        return ForTestingInterface(**options)

    def test1(self):
        interface=self.ForTesting(redirection="none",number_of_workers=1)        
        x=interface.grid.x
        dens=interface.grid.dens
        self.assertEqual(x,numpy.arange(60)/60.)
        self.assertEqual(dens,x**2)
        interface.stop()
    
    def test2(self):
        for n in [3,5,6]:
            interface=self.ForTesting(redirection="none",number_of_workers=n)        
            x=interface.grid.x
            dens=interface.grid.dens
            self.assertEqual(x,numpy.arange(60)/60.)
            self.assertEqual(dens,x**2)
            interface.stop()

