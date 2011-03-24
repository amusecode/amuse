"""
 example of calculating Pythagorean problem using low level interface
 
"""

import numpy
import time

from amuse.community.phiGRAPE.interface import PhiGRAPEInterface as phi
from amuse.community.hermite0.interface import HermiteInterface as her
from amuse.community.huayno.interface import HuaynoInterface as hua
from amuse.community.bhtree.interface import BHTreeInterface as bht

from matplotlib import pyplot

def pyth(interface):

  mass=[3.,4.,5.]
  
  x=[1.,-2.,1]
  y=[3,-1.,-1.]
  z=[0.,0.,0]
  
  vx=[0.,0.,0.]
  vy=[0.,0.,0.]
  vz=[0.,0.,0.]
  
  radius=[0.,0.,0.]
  
  nb = interface()
  nb.initialize_code()

  ids,error = nb.new_particle(mass,radius,x,y,z,vx,vy,vz)
  if filter(lambda x: x != 0, error) != []: raise Exception
  return nb,ids

def run_pyth(interface,tend=100,dt=0.125,parameters=[]):
  
  nb,ids=pyth(interface)
  for p,val in parameters:
     (eval("nb.set_"+p))(val)    
  nb.commit_particles()

  m,r,x,y,z,vx,vy,vz,err=nb.get_state(ids)
  xx=[x]
  yy=[y]
  
  t=0.
  while(t<tend-dt/2):
    t=t+dt
    nb.evolve(t)
    m,r,x,y,z,vx,vy,vz,err=nb.get_state(ids)
    xx.append(x)
    yy.append(y)
  nb.stop()
     
  return numpy.array(xx), numpy.array(yy)

if __name__ in ['__main__', '__plot__']:
  codes_to_run=[ ('Hermite0, $\eta=0.03$', her,  [("dt_param",0.03)] ),
                 ('Hermite0, $\eta=0.01$', her,  [("dt_param",0.01)] ),
                 ('Hermite0, $\eta=0.003$', her,  [("dt_param",0.003)] ),
                 ('Hermite0, $\eta=0.001$', her,  [("dt_param",0.001)] ) ]
  N=(len(codes_to_run)-1)/2+1
  f=pyplot.figure(figsize=(8,4*N))
  
  for i,(label,code,parameters) in enumerate(codes_to_run):
    xx,yy=run_pyth(code,tend=100,dt=0.0625,parameters=parameters)
    subplot=f.add_subplot(N,2,i+1)
    subplot.plot(xx[:,0],yy[:,0],'r')
    subplot.plot(xx[:,1],yy[:,1],'b')
    subplot.plot(xx[:,2],yy[:,2],'g')
    subplot.set_title(label)
    subplot.set_xlim(-8,8)
    subplot.set_ylim(-6,6)
  pyplot.show()


  
