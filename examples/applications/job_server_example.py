from amuse.ext.job_server import JobServer
from time import sleep

def somework(x):
  sleep(0.5)
  return x*x

def example_parallel_jobs_1(N):

  from socket import gethostname
  
  jobserver=JobServer(["localhost"]*2,channel_type="mpi")
  print
  for i in range(N):
    jobserver.submit_job(somework, (i,))
   
  i=0 
  while jobserver.wait():
    job=jobserver.last_finished_job
    print job.args[0],job.result
  
def example_parallel_jobs_2(N):

  from socket import gethostname
  
  jobserver=JobServer(["localhost"]*2,channel_type="mpi")
  print
  for i in range(N):
    jobserver.submit_job(somework, (i,)) 
     
  print "waiting"  
  jobserver.waitall()
  print "done"
  for job in jobserver.finished_jobs:
    print job.args[0],job.result  

    
if __name__=="__main__":
# this is needed in order for all the functions and classes to be 
# pickled with full module name
  from job_server_example import example_parallel_jobs_1,example_parallel_jobs_2  
  example_parallel_jobs_1(10)
  example_parallel_jobs_2(10)
