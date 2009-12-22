import copy,numpy
import numpy as math

# tbd: initial value for time_Radius solver

def stumpff_C(z):
  if(z==0): 
    return 1/2.
  if(z>0):
    sz=math.sqrt(z)  
    return (1-math.cos(sz))/z  
  if(z<0):
    sz=math.sqrt(-z)  
    return -(math.cosh(sz)-1)/z  

def stumpff_S(z):
  if(z==0):
    return 1/6.
  if(z>0):
    sz=math.sqrt(z)  
    return (sz-math.sin(sz))/sz**3  
  if(z<0):
    sz=math.sqrt(-z)  
    return (math.sinh(sz)-sz)/sz**3  
    
def stumpff_C_prime(z):
  return (stumpff_C(z)-3*stumpff_S(z))/2/z 
def stumpff_S_prime(z):
  return (1-stumpff_S(z)-2*stumpff_C(z))/2/z 

def lagrange_f(xi,r0,vr0,smu,alpha):
  return 1.-xi**2/r0*stumpff_C(alpha*xi**2)    
def lagrange_dfdxi(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
  return xi/r0*(z*stumpff_S(z)-1)    
def lagrange_g(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
  return r0*vr0*(xi/smu)**2*stumpff_C(z)-r0*xi*z/smu*stumpff_S(z)+r0*xi/smu
def lagrange_dgdxi(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
  return r0*vr0/smu*(xi/smu)*(1-z*stumpff_S(z))-z*r0/smu*stumpff_C(z)+r0/smu

def universal_kepler(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
  return r0*vr0*xi**2*stumpff_C(z)/smu+ \
         (1-alpha*r0)*xi**3*stumpff_S(z)+r0*xi

def universal_kepler_dxi(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
  return r0*vr0*xi*(1-alpha*xi**2*stumpff_S(z))/smu + \
         (1-alpha*r0)*xi**2*stumpff_C(z)+r0  

def universal_kepler_dxidxi(xi,r0,vr0,smu,alpha):
  z=alpha*xi**2
#  return r0*vr0/smu-alpha*r0*vr0*xi**2/smu*stumpff_C(z)+\
#         (1-alpha*r0)*xi*(1-z*stumpff_S(z))
  return -alpha*universal_kepler(xi,r0,vr0,smu,alpha)+r0*vr0/smu+xi
      
def newton(f,x0,fprime=None,args=(),tol=1.48e-8,maxiter=50):
  if fprime is None:
    print "provide fprime"
    return x0
  i=0
  x=x0
  x1=0  
  while (i<maxiter):
    fv=f(x,*args)
    dfv=fprime(x,*args)
    if(dfv==0):
      return x0,-2
    delta=-fv/dfv
    if(abs(delta)<tol): 
      return x+delta,0
    x=x+delta
    i=i+1
  return x,-1    

def universal_time_radius_solver(radius,mu,pos0,vel0,dt):
  r02=reduce(lambda x,y: x+ y**2,pos0,0)
  v02=reduce(lambda x,y: x+ y**2,vel0,0)
  v0r0=(reduce(lambda x,y: x+y, pos0*vel0,0))
  r0=math.sqrt(r02)
  v0=math.sqrt(v02)
  vr0=v0r0/r0
  h2=r02*v02-(v0r0)**2
  p=h2/mu
  alpha=2./r0-v0**2/mu
  rp=p/(1+math.sqrt(1-p*alpha))
  if(radius < rp):
    return dt,0 
  if(r0==radius):
    return 0.,1

  def f(xi):
    return universal_kepler_dxi(xi,r0,vr0,smu,alpha)-radius
  def df(xi):
    return universal_kepler_dxidxi(xi,r0,vr0,smu,alpha)
      
  smu=math.sqrt(mu)
  xi0=1/math.sqrt(alpha)
  xi,err=newton(f,xi0,fprime=df,tol=1.e-10)  
  dt_coll=universal_kepler(xi,r0,vr0,smu,alpha)/smu
  if(dt_coll > 0 and dt_coll < dt):
    return dt_coll,1
  else:
    return dt,0
    
def collision(radius,mu,pos0,vel0,dt):
  return universal_time_radius_solver(radius,mu,pos0,vel0,dt)

def universal_solver(mu,pos0,vel0,dt):
  smu=math.sqrt(mu)
  
  r0=math.sqrt(reduce(lambda x,y: x+ y**2,pos0,0))
  v0=math.sqrt(reduce(lambda x,y: x+ y**2,vel0,0))
  vr0=(reduce(lambda x,y: x+y, pos0*vel0,0))/r0
  alpha=2./r0-v0**2/mu
    
  xi0=smu*abs(alpha)*dt
  
  def f(xi):
    return universal_kepler(xi,r0,vr0,smu,alpha)-smu*dt
  def df(xi):
    return universal_kepler_dxi(xi,r0,vr0,smu,alpha)
  
  xi,err=newton(f,xi0,fprime=df,tol=1.e-10)  
  
  pos=pos0*lagrange_f(xi,r0,vr0,smu,alpha)+vel0*lagrange_g(xi,r0,vr0,smu,alpha)
  r=math.sqrt(reduce(lambda x,y: x+ y**2,pos,0))
  vel=pos0*smu/r*lagrange_dfdxi(xi,r0,vr0,smu,alpha)+\
      vel0*smu/r*lagrange_dgdxi(xi,r0,vr0,smu,alpha)
  return pos,vel

class twobody(object):
  __G=6.673e-11
  def __init__(self):
    self.particles=[]
    self.tnow=0.
  def new_particle(self,mass,radius,x,y,z,vx,vy,vz):
    if( len(self.particles)>=2):
      return 0,-1
    self.particles.append( {'mass': mass, 'radius' : radius, \
                            'x' : x , 'y' : y, 'z' : z, \
                           'vx' : vx ,'vy' : vy,'vz' : vz })  
    return len(self.particles)-1,0
  def get_state(self,pid):
    try:
      return copy.deepcopy(self.particles[pid]),0
    except:        
      return {},-1
  def get_time(self):
    return self.tnow,0
  def evolve(self,time_end):
      if(len(self.particles)!=1 and len(self.particles)!=2):
        return -1
      if(len(self.particles)==1):
        mu=self.__G*self.particles[0]['mass']
        radius=self.particles[0]['radius']        
        dpos_initial=numpy.array( [self.particles[0]['x'], \
                         self.particles[0]['y'], \
                         self.particles[0]['z']] )
        dvel_initial=numpy.array( [self.particles[0]['vx'], \
                         self.particles[0]['vy'], \
                         self.particles[0]['vz']] )
        dt,collisionflag=collision(radius,mu,dpos_initial,dvel_initial,time_end-self.tnow)
        dpos,dvel=universal_solver(mu,dpos_initial,dvel_initial,dt)
        self.particles[0]['x']=dpos[0]
        self.particles[0]['y']=dpos[1]
        self.particles[0]['z']=dpos[2]
        self.particles[0]['vx']=dvel[0]
        self.particles[0]['vy']=dvel[1]
        self.particles[0]['vz']=dvel[2]                         
      if(len(self.particles)==2):
        mu=self.__G*(self.particles[0]['mass']+self.particles[1]['mass'])
        radius=self.particles[0]['radius']+self.particles[1]['radius']        
        tm=(self.particles[0]['mass']+self.particles[1]['mass'])
        m0=self.particles[0]['mass']
        m1=self.particles[1]['mass']
        pos0=numpy.array( [self.particles[0]['x'], \
                           self.particles[0]['y'], \
                           self.particles[0]['z']] )
        vel0=numpy.array( [self.particles[0]['vx'], \
                           self.particles[0]['vy'], \
                           self.particles[0]['vz']] )
        pos1=numpy.array( [self.particles[1]['x'], \
                           self.particles[1]['y'], \
                           self.particles[1]['z']] )
        vel1=numpy.array( [self.particles[1]['vx'], \
                           self.particles[1]['vy'], \
                           self.particles[1]['vz']] )
        dpos_initial=pos0-pos1
        dvel_initial=vel0-vel1
        cmpos=(m0*pos0+m1*pos1)/tm
        cmvel=(m0*vel0+m1*vel1)/tm        
        dt,collisionflag=collision(radius,mu,dpos_initial,dvel_initial,time_end-self.tnow)
        dpos,dvel=universal_solver(mu,dpos_initial,dvel_initial,dt)
        cmpos=cmpos+(time_end-self.tnow)*cmvel
        f0=m1/tm
        f1=m0/tm
        self.particles[0]['x']=cmpos[0]+f0*dpos[0]
        self.particles[0]['y']=cmpos[1]+f0*dpos[1]
        self.particles[0]['z']=cmpos[2]+f0*dpos[2]
        self.particles[0]['vx']=cmvel[0]+f0*dvel[0]
        self.particles[0]['vy']=cmvel[1]+f0*dvel[1]
        self.particles[0]['vz']=cmvel[2]+f0*dvel[2]                         
        self.particles[1]['x']=cmpos[0]-f1*dpos[0]
        self.particles[1]['y']=cmpos[1]-f1*dpos[1]
        self.particles[1]['z']=cmpos[2]-f1*dpos[2]
        self.particles[1]['vx']=cmvel[0]-f1*dvel[0]
        self.particles[1]['vy']=cmvel[1]-f1*dvel[1]
        self.particles[1]['vz']=cmvel[2]-f1*dvel[2]                         
      self.tnow=self.tnow+dt
      return collisionflag

if __name__=='__main__':
  nb=twobody()
  nb.new_particle(5.9742e24,6.371e6,3.85e8,0.,0.,0.,0.,0.)
  err=nb.evolve(3600.*24*30)
  dt,err=nb.get_time()
  print dt/24/3600.
