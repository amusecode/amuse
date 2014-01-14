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

2 issues to be fixed:
 - blocking startup of hosts may prevent threads shutting down, leading to freeze at end of script
   (so manual kill necessary) 
 - thread function _startup contains references to JobServer, hence del jobserver will actually not be called
   until the situation of issue 1 is resolved (so the warning given there is useless) 

"""
from amuse.rfi.core import *
import cPickle as pickle
from amuse.rfi.channel import AsyncRequestsPool
import inspect
from collections import deque
import threading
from time import sleep

def dump_and_encode(x):
  return pickle.dumps(x,-1)
def decode_and_load(x):
  return pickle.loads(x.encode("latin-1"))

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
      self.result=None
      self.request=None
      self.err=None

class JobServer(object):
    def __init__(self,hosts=[],channel_type="mpi",preamble=None, retry_jobs=True, 
                   no_wait=True,verbose=True):
      self.hosts=[]
      self.job_list=deque()
      self.idle_codes=[]
      self.retry_jobs=retry_jobs
      self._finished_jobs=deque()
      self.preamble=preamble
      self.pool=AsyncRequestsPool()
      self.number_available_codes=0
      self.number_starting_codes=0
      self.no_wait=no_wait
      self.last_finished_job=None
      self.verbose=verbose
      if self.verbose:
          print "AMUSE JobServer launching"

      self.add_hosts(hosts=hosts,channel_type=channel_type)
      
   
    def add_hosts(self,hosts=[],channel_type="mpi"):
      self.hosts.append(hosts)
      if self.verbose:
        print "JobServer: connecting %i hosts"%len(hosts),
      if channel_type=="mpi":
        for host in hosts:
          self.number_starting_codes+=1
          self._startup( channel_type=channel_type,hostname=host,
                           copy_worker_code=True,redirection="none" )
      else:  
        threads=[]
        for host in hosts:
          kwargs=dict( channel_type=channel_type,hostname=host,
                         copy_worker_code=True,redirection="none" )
          threads.append( threading.Thread(target=self._startup,kwargs=kwargs) )
        for thread in threads:
          self.number_starting_codes+=1
          thread.daemon=True
          thread.start()
        if not self.no_wait:  
          if self.verbose:
            print "... waiting"
          for thread in threads:
            thread.join()
        else:
          if self.verbose:
            print "... waiting for first available host"
          while self.number_available_codes==0 and self.number_starting_codes>0:
            sleep(0.1)
      if self.no_wait:
        if self.verbose:
          print "JobServer: launched"
      else:    
        if self.verbose:
          print "JobServer: launched with", len(self.idle_codes),"hosts"
    
    def _startup(self, *args,**kwargs):
      try: 
        code=CodeInterface(*args,**kwargs) 
      except Exception as ex:
        self.number_starting_codes-=1
        print "JobServer: startup failed on", kwargs['hostname']
        print ex
      else:
        if self.preamble is not None:
          code.exec_(self.preamble)
           
        self.number_available_codes+=1
        if self.no_wait:
          if self.number_available_codes & (self.number_available_codes-1) ==0:
            if self.verbose:
              print "JobServer: hosts now available:",self.number_available_codes
          self.number_starting_codes-=1
          if self.number_starting_codes==0:
            if self.verbose:
              print "JobServer: hosts in total:", self.number_available_codes
        if self.job_list: 
          self._add_job(self.job_list.popleft(), code)
        else:
          self.idle_codes.append(code)   
  
    def exec_(self,arg):
      while self.number_starting_codes>0:
        sleep(0.1)
      self.waitall()  
      for code in self.idle_codes:
        code.exec_(arg)
    
    def submit_job(self,f,args=(),kwargs={}):
      if len(self.pool)==0 and not self.job_list:
        if self.verbose:
          print "JobServer: submitting first job on queue"
      job=Job(f,args,kwargs)
      self.job_list.append( job)
      if self.idle_codes: 
          self._add_job(self.job_list.popleft(), self.idle_codes.pop())        
      return job

    def wait(self):
      if self._finished_jobs:
        self.last_finished_job=self._finished_jobs.popleft()
        return True
      elif len(self.pool)==0 and not self.job_list:
        if self.verbose:
          print "JobServer: no more jobs on queue or running"        
        return False
      else:
        while len(self.pool)==0 and self.job_list:
          sleep(0.1)
          if self.number_available_codes>0:
            raise Exception("JobServer: this should not happen")    
          if self.number_starting_codes==0:
            raise Exception("JobServer: no codes available")
        self.pool.wait()
        self.last_finished_job=self._finished_jobs.popleft()
        return True

    def waitall(self):
      while len(self.pool)>0 or self.job_list:
        self.pool.wait()
    
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
        self.number_available_codes-=1
      else:
        if self.job_list:
          self._add_job( self.job_list.popleft(), code)
          if not self.job_list:
            if self.verbose:
              print "JobServer: last job dispatched"
        else:
          self.idle_codes.append(code)
      self._finished_jobs.append(job)
    
    def _add_job(self,job,code):
      job.request=code.async_func(job.f,*job.args,**job.kwargs)
      self.pool.add_request(job.request,self._finalize_job, [job,code])
    
    def __del__(self):
      self.waitall()
      if self.job_list:
        print "JobServer: Warning: unfinished jobs"
      for code in self.idle_codes:
        code.stop()
      if self.number_starting_codes>0:
        print "JobServer: Warning: some hosts startup threads possibly blocking"
