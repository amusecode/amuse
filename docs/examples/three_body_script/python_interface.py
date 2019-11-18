from amuse.rfi.core import *
import cPickle as Pickle
from amuse.units import units

#from testmod import f
from bin_planet import binary_with_planet_run as f

def func(*arg,**kwarg):
  return f(*arg,**kwarg)

class CodeImplementation(object):
   def _func(self,argin,kwargin,argout):
     try:
       arg=Pickle.loads(argin)
       kwarg=Pickle.loads(kwargin)
       result=func(*arg,**kwarg)
       argout.value=Pickle.dumps(result,-1)
       return 0
     except Exception as ex:
       print ex
       argout.value=Pickle.dumps(" ",-1)
       return -1

class CodeInterface(PythonCodeInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, CodeImplementation, **options)
    
    @legacy_function
    def _func():
        function = LegacyFunctionSpecification()
        function.addParameter('argin', dtype='string', direction=function.IN)
        function.addParameter('kwargin', dtype='string', direction=function.IN)
        function.addParameter('argout', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    def func(self,*args,**kwargs):
        result,err=self._func(Pickle.dumps( args,-1),Pickle.dumps( kwargs,-1))
        return Pickle.loads(result[0]),err

    def async_func(self,*args,**kwargs):
        request=self._func.async(Pickle.dumps( args,-1),Pickle.dumps( kwargs,-1))
        def f(x):
          result,err=x()
          return Pickle.loads(result[0]),err
        request.add_result_handler( f )
        return request
