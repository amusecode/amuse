"""
JobServer: class for job farming using amuse communication channels

usage:

start jobserver
  jobserver=JobServer(hosts=<list of hostnames> [ ,channel_type="mpi", preamble="<commands>", retry_jobs=True/ False] )

submit job
  job=jobserver.submit_job(somework, (args,))

wait for one result (encounters all):
  jobserver.wait()
  job=jobserver.last_finished_job

wait for all to finish, loop over all:
  jobserver.waitall()
  for job in jobserver.finished_jobs:
    print job.result
    
it is essential that the function which are to be executed remotely are pickleable, i.e. they must not be 
derived from the main module. Easy way to achieve this is to import them from a seperate file.    

"""


from amuse.rfi.core import *
import cPickle as pickle
from amuse.rfi.channel import AsyncRequestsPool
import inspect
from collections import deque
import threading

def dump_and_encode(x):
  return pickle.dumps(x,0)
def decode_and_load(x):
  return pickle.loads(x)

class CodeImplementation(object):
   def __init__(self):
     self.scope={}
     self.scope['dump_and_encode']=dump_and_encode
     self.scope['decode_and_load']=decode_and_load
   
   def exec_(self,arg):
     try:
       exec arg in self.scope
       return 0
     except Exception as ex:  
       print ex
       return -1
   def _func(self,f,argin,kwargin,argout):
     try:
       self.scope.update(dict(f=f,argin=argin,kwargin=kwargin))
       exec "func=decode_and_load(f)" in self.scope
       exec "arg=decode_and_load(argin)" in self.scope
       exec "kwarg=decode_and_load(kwargin)" in self.scope
       exec "result=func(*arg,**kwarg)" in self.scope
       argout.value=eval("dump_and_encode(result)",self.scope)
       return 0
     except Exception as ex:
       argout.value=dump_and_encode(ex)
       return -1

class CodeInterface(PythonCodeInterface):    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, CodeImplementation, **options)
    
    @legacy_function
    def _func():
        function = LegacyFunctionSpecification()
        function.addParameter('func', dtype='string', direction=function.IN)
        function.addParameter('argin', dtype='string', direction=function.IN)
        function.addParameter('kwargin', dtype='string', direction=function.IN)
        function.addParameter('argout', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def exec_():
        function = LegacyFunctionSpecification()
        function.addParameter('arg', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    def func(self,f,*args,**kwargs):
        result,err=self._func( dump_and_encode(f),
                               dump_and_encode(args),
                               dump_and_encode(kwargs) )
        return decode_and_load(result[0]),err

    def async_func(self,f,*args,**kwargs):
        request=self._func.async(dump_and_encode(f),
                                 dump_and_encode(args),
                                 dump_and_encode(kwargs) )
        def f(x):
          result,err=x()
          return decode_and_load(result[0]),err
        request.add_result_handler( f )
        return request

class Job(object):
    def __init__(self, f, args, kwargs):
      self.f=f
      self.args=args
      self.kwargs=kwargs
      self._result=None
      self.request=None
      self.err=None

class JobServer(object):
    def __init__(self,hosts,channel_type="mpi",preamble=None, retry_jobs=True):
      self.hosts=hosts
      self.job_list=deque()
      self.idle_codes=[]
      self.channel_type=channel_type
      self.retry_jobs=retry_jobs
      self._finished_jobs=deque()
      print "connecting %i hosts"%len(hosts),
      threads=[]
      for host in hosts:
        kwargs=dict( channel_type=self.channel_type,hostname=host,
                     copy_worker_code=True,redirection="none" )
        threads.append( threading.Thread(target=self._startup,kwargs=kwargs) )
      for thread in threads:
        thread.start()
      print "... waiting"
      for thread in threads:
        thread.join()               
      if preamble is not None:
        for code in self.idle_codes:
          code.exec_(preamble)          
      self.pool=AsyncRequestsPool()
      print "\nAMUSE JobServer launched with", len(self.idle_codes),"threads"
    
    def _startup(self, *args,**kwargs):
      try: 
        code=CodeInterface(*args,**kwargs) 
      except Exception as ex:
        print
        print "startup failed on", kwargs['hostname']
        print ex
      else:
          self.idle_codes.append(code)      
    
    def submit_job(self,f,args=(),kwargs={}):
      job=Job(f,args,kwargs)
      self.job_list.append( job)
      if self.idle_codes: 
          self._add_job(self.job_list.popleft(), self.idle_codes.pop())        
      return job

    def wait(self):
      if self._finished_jobs:
        return True
      elif len(self.pool)==0:
        return False
      else:
        self.pool.wait()
        return True

    def waitall(self):
      while len(self.pool)>0:
        self.pool.wait()

    @property
    def last_finished_job(self):
      if self._finished_jobs:
        return self._finished_jobs.popleft()
      else:
        return None
    
    @property
    def finished_jobs(self):
       while self._finished_jobs:
         yield self._finished_jobs.popleft()
    
    def _finalize_job(self,request,job,code):
      try: 
        job.result,job.err=request.result()
      except Exception as ex:
        job.result,job.err=ex,-2
        if self.retry_jobs:
          self.job_list.append( job)
        del code
      else:
        if self.job_list:
          self._add_job( self.job_list.popleft(), code)
        else:
          self.idle_codes.append(code)  
      self._finished_jobs.append(job)
    
    def _add_job(self,job,code):
      job.request=code.async_func(job.f,*job.args,**job.kwargs)
      self.pool.add_request(job.request,self._finalize_job, [job,code])
    
    def __del__(self):
      self.waitall()
      if self.job_list:
        print "Warning: unfinished jobs in JobServer"
      for code in self.idle_codes:
        code.stop()
