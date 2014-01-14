from amuse.test import amusetest

from amuse.ext.job_server import JobServer

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

class TestJobServer(amusetest.TestCase):

    def test1(self):
                  
      result=example_parallel_jobs(10)
      for arg,res in result.items():
        self.assertEqual(arg,len(res))

    def test2(self):
                  
      result=example_parallel_jobs(10,4)
      for arg,res in result.items():
        self.assertEqual(arg,len(res))
