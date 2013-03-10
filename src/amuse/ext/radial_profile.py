import numpy

def radial_profile(r,dat,N=100):
  n=len(r)
  a=r.argsort()
  i=0
  r_a=[]
  dat_a=[]
  while i < n:
    ra=r[a[i:i+N]].sum()/min(n-i,N)
    da=dat[a[i:i+N]].sum()/min(n-i,N)
    r_a.append(ra)
    dat_a.append(da)
    i=i+N
  return numpy.array(r_a),numpy.array(dat_a)

def radial_density(r,mass,N=100,dim=3, start_at_zero=False):
  if dim==3:
    volfac=numpy.pi*4./3.
  elif dim==2:
    volfac=numpy.pi
  else:
    volfac=1
  
  n=len(r)
  a=r.argsort()
  i=0
  r_a=[]
  dens=[]
  oldrshell=r[a[0]]
  if start_at_zero:
    oldrshell=0.*r[0]  
  while i < n:
    i1=min(n,i+N)
    rshell=r[a[i1-1]]
    ra=r[a[i:i1]].sum()/(i1-i)
    da=mass[a[i:i1]].sum()/(rshell**dim-oldrshell**dim)
    oldrshell=rshell
    r_a.append(ra)
    dens.append(da)
    i=i1
    
  return numpy.array(r_a),numpy.array(dens)/volfac

if __name__=="__main__":
  from matplotlib import pyplot
  from amuse.ic.plummer import MakePlummerModel
  
  plum=MakePlummerModel(100000).result
  r=(plum.x**2+plum.y**2+plum.z**2)**0.5
  ra,dens=radial_density(r,plum.mass,100,start_at_zero=True)
  
  ascl=1/1.695
  
  ra=numpy.array(map(lambda x: x.number,ra))
  dens=numpy.array(map(lambda x: x.number,dens))
  pyplot.subplot(211)
  pyplot.loglog(ra,dens)
  pyplot.loglog(ra, 3./4./numpy.pi/ascl**3/(1+(ra**2/ascl**2))**(5./2))
#  pyplot.plot(ra,(dens-3./4./numpy.pi/ascl**3/(1+(ra**2/ascl**2))**(5./2))/dens,'r.')
  pyplot.subplot(212)
  pyplot.plot(plum.x.number,plum.y.number,'r.')
  pyplot.show()
