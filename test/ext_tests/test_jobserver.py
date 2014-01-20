from amuse.test import amusetest

from amuse.ext.job_server import RemoteCodeInterface,JobServer

from socket import gethostname

from amuse.ic.plummer import new_plummer_model

def example_parallel_jobs(N,Nhosts=1):
    
  jobserver=JobServer(hosts=[gethostname()]*Nhosts)
    
  for i in range(N):
    jobserver.submit_job(new_plummer_model, (i,))
  
  result=dict()
  while jobserver.wait():
    job=jobserver.last_finished_job
    result[job.args[0]]=job.result

  return result

def example_parallel_jobs2(N,Nhosts=1):
    
  jobserver=JobServer(hosts=[gethostname()]*Nhosts)
    
  for i in range(N):
    jobserver.submit_job(new_plummer_model, (i,))
  
  jobserver.waitall()
  
  result=dict()
  
  for job in jobserver.finished_jobs:
    result[job.args[0]]=job.result

  return result

class TestRemoteCode(amusetest.TestCase):

    def test1(self):
      remote=RemoteCodeInterface()
      
      var=123
      remote.assign("var",var)
      var_=remote.evaluate("var+1")
      self.assertEqual(var_,var+1)
      remote.execute("var=var*2")
      var_=remote.evaluate("var")
      self.assertEqual(var_,var*2)
      
    def test2(self):
      remote=RemoteCodeInterface()
      
      from math import sqrt
      result=remote.func(sqrt, 64)
      self.assertEqual(result,8)

    def test3(self):
      remote=RemoteCodeInterface()
      
      var=new_plummer_model(100)
      remote.assign("var",var)
      var_=remote.evaluate("var.mass")
      self.assertEqual(var_,var.mass)

    def test4(self):
      remote=RemoteCodeInterface(channel_type="sockets")
      
      var=123
      remote.assign("var",var)
      var_=remote.evaluate("var+1")
      self.assertEqual(var_,var+1)
      remote.execute("var=var*2")
      var_=remote.evaluate("var")
      self.assertEqual(var_,var*2)


class TestJobServer(amusetest.TestCase):

    def test1(self):
                  
      result=example_parallel_jobs(10)
      for arg,res in result.items():
        self.assertEqual(arg,len(res))

    def test2(self):
                  
      result=example_parallel_jobs(10,4)
      for arg,res in result.items():
        self.assertEqual(arg,len(res))

    def test3(self):
                  
      result=example_parallel_jobs2(10,4)
      for arg,res in result.items():
        self.assertEqual(arg,len(res))
