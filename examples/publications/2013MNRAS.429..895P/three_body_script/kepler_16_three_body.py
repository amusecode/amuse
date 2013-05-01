from python_interface import CodeInterface
from amuse.units import units
from amuse.rfi.channel import AsyncRequestsPool
import numpy

from matplotlib import pyplot

ncpu=16
allhosts=[ ("paddegat",4),
           ("koppoel",4),
           ("gaasp",4),
           ("biesbosch",4),
           ]

hosts=reduce(lambda x,y: x+[y[0]]*y[1],allhosts,[])
print hosts

m_planet=0.333 | units.MJupiter
r_planet=0.754 | units.RJupiter,

m1=0.6897| units.MSun      # primary mass
m2=0.20255| units.MSun      # secondary mass
r1=0.6489 | units.RSun      # primary radius
r2=0.22623| units.RSun      # secondary radius
ecc_binary=0.15944        # binary orbit eccentric$
P_binary=41.08| units.day  # binary orbit period
a_binary=0.2243 # units.AU

amin=2.
amax=10.

results=[]

N=40
M=40

aset= ( a_binary*(amin+j*(amax-amin)/(M+1)) for j in range(M+1) )
eset= [ j*0.5/N for j in range(N+1)]
current_a=0

pool=AsyncRequestsPool()

def finalize_job(request,i_ecc,a,code,host):
  print "done with", eset[i_ecc],a
  result,err=request.result()
  print result
  results.append((eset[i_ecc],a,result))

  if result[0]=="stable" or i_ecc==0:
#  if i_ecc==0:
    try:
      a=aset.next()
    except:
      a=None
    if a is not None:  
      i_ecc=N
      add_job(i_ecc,a,code,host)
  else:
    i_ecc-=1
    add_job(i_ecc,a,code,host)

def add_job(i_ecc,a,code=None,host=""):
    ecc=eset[i_ecc]
    if code is None:
      host=hosts.pop()
      code=CodeInterface(hostname=host)
    print "adding:",(ecc,a),' on: ',host
    request=code.async_func(
      m1=m1,m2=m2,m_planet=m_planet,
      r1=r1,r2=r2,r_planet=r_planet,
      ecc_binary=ecc_binary,P_binary=P_binary,ecc_planet=ecc,a_planet=a | units.AU,
      pangle_planet=numpy.pi, tend=1000000.| units.yr,hostname=host) 
    pool.add_request( request, finalize_job, [i_ecc,a,code,host])

for i in range(ncpu):
    try:
      a=aset.next()
    except:
      a=None
    if a is not None:  
      i_ecc=N
      add_job(i_ecc,a)

while len(pool)>0:

  print "waiting for job to finish.."
  pool.wait()

f=open('results_16_pi','w')
import cPickle
cPickle.dump(results,f)
f.close()
