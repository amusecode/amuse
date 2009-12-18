import copy,numpy
import numpy as math

# note: collision detection TBD

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

def universal_kepler(xi,r0,vr0,smu,alpha,dt):
  z=alpha*xi**2
  return r0*vr0*xi**2*stumpff_C(z)/smu+ \
         (1-alpha*r0)*xi**3*stumpff_S(z)+r0*xi-smu*dt

def universal_kepler_derivative(xi,r0,vr0,smu,alpha,dt):
  z=alpha*xi**2
  return r0*vr0*xi*(1-alpha*xi**2*stumpff_S(z))/smu + \
         (1-alpha*r0)*xi**2*stumpff_C(z)+r0  
    
def test_stumpff():
  print stumpff_C(0.),stumpff_C(0.001),stumpff_C(-0.001)
  print stumpff_S(0.),stumpff_S(0.001),stumpff_S(-0.001)

def newton(f,x0,fprime=None,args=None,tol=1.48e-8,maxiter=50):
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

def universal_solver(mu,pos0,vel0,dt):
  smu=math.sqrt(mu)
  
  r0=math.sqrt(reduce(lambda x,y: x+ y**2,pos0,0))
  v0=math.sqrt(reduce(lambda x,y: x+ y**2,vel0,0))
  vr0=(reduce(lambda x,y: x+y, pos0*vel0,0))/r0
  alpha=2./r0-v0**2/mu
    
  xi0=smu*abs(alpha)*dt
  
  f=universal_kepler
  df=universal_kepler_derivative
  
  xi,err=newton(f,xi0,fprime=df,args=(r0,vr0,smu,alpha,dt),tol=1.e-10)  

  lagrange_f=1.-xi**2/r0*stumpff_C(alpha*xi**2)    
  lagrange_g=dt-1/smu*xi**3*stumpff_S(alpha*xi**2)
  
  pos=pos0*lagrange_f+vel0*lagrange_g
  r=math.sqrt(reduce(lambda x,y: x+ y**2,pos,0))
  lagrange_df=math.sqrt(mu)/r/r0*(alpha*xi**3*stumpff_S(alpha*xi**2)-xi)
  lagrange_dg=1-xi**2/r*stumpff_C(alpha*xi**2)    
  vel=pos0*lagrange_df+vel0*lagrange_dg

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
  def evolve(self,time_end):
      if(len(self.particles)!=1 and len(self.particles)!=2):
        return -1
      if(len(self.particles)==1):
        mu=self.__G*self.particles[0]['mass']
        dpos_initial=numpy.array( [self.particles[0]['x'], \
                         self.particles[0]['y'], \
                         self.particles[0]['z']] )
        dvel_initial=numpy.array( [self.particles[0]['vx'], \
                         self.particles[0]['vy'], \
                         self.particles[0]['vz']] )
        dpos,dvel=universal_solver(mu,dpos_initial,dvel_initial,time_end-self.tnow)
        self.particles[0]['x']=dpos[0]
        self.particles[0]['y']=dpos[1]
        self.particles[0]['z']=dpos[2]
        self.particles[0]['vx']=dvel[0]
        self.particles[0]['vy']=dvel[1]
        self.particles[0]['vz']=dvel[2]                         
      if(len(self.particles)==2):
        mu=self.__G*(self.particles[0]['mass']+self.particles[1]['mass'])
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
        dpos,dvel=universal_solver(mu,dpos_initial,dvel_initial,time_end-self.tnow)
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
      self.tnow=time_end
      return 0

