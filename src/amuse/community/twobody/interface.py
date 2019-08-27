from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode
import copy,numpy
import numpy as math

from amuse.units.quantities import zero

# tbd: initial value for time_Radius solver

from amuse.datamodel import Particles
from amuse.rfi.core import *
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
    return (r0*vr0*xi**2*stumpff_C(z)/smu+  
           (1-alpha*r0)*xi**3*stumpff_S(z)+r0*xi)

def universal_kepler_dxi(xi,r0,vr0,smu,alpha):
    z=alpha*xi**2
    return (r0*vr0*xi*(1-alpha*xi**2*stumpff_S(z))/smu +  
           (1-alpha*r0)*xi**2*stumpff_C(z)+r0)

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

def laguerre(f,x0,fprime=None,fprimeprime=None,args=(),order=4,tol=1.e-14,maxiter=50):
    if fprime is None:
        print "provide fprime"
        return x0
    if fprimeprime is None:
        print "provide fprimeprime"
        return x0
    i=0
    x=x0
    while (i<maxiter):
        fv=f(x,*args)
        dfv=fprime(x,*args)
        ddfv=fprimeprime(x,*args)
        if(dfv==0 or ddfv==0):
            return x0,-2            
        delta=-order*fv/(dfv+math.sign(dfv)*
          math.abs((order-1)**2*dfv**2-order*(order-1)*fv*ddfv)**0.5)
        if(abs(delta)<tol): 
            return x+delta,0
        x=x+delta
        i=i+1
    return x,-1    


def universal_time_radius_solver(radius,mu,pos0,vel0,dt):
    r02=reduce(lambda x,y: x+ y**2,pos0,0)
    v02=reduce(lambda x,y: x+ y**2,vel0,0)
    v0r0=reduce(lambda x,y: x+y, pos0*vel0,0)
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
    if(alpha>0):
        ra=p/(1-math.sqrt(1-p*alpha))
        if(radius > ra):
            return dt,0   

    smu=math.sqrt(mu)
    xi0=1/math.sqrt(alpha)
        
    def f(xi):
        return universal_kepler_dxi(xi,r0,vr0,smu,alpha)-radius
    def df(xi):
        return universal_kepler_dxidxi(xi,r0,vr0,smu,alpha)
        
    xi,err=newton(f,xi0,fprime=df,tol=1.e-10)    

    dt_coll=universal_kepler(xi,r0,vr0,smu,alpha)/smu
    if(dt_coll > 0 and dt_coll < dt):
        return dt_coll,1
    else:
        return dt,0
    
def collision(radius,mu,pos0,vel0,dt):
    return universal_time_radius_solver(radius,mu,pos0,vel0,dt)

def universal_solver_newton(mu,pos0,vel0,dt):
    smu=math.sqrt(mu)
    
    r0=math.sqrt(reduce(lambda x,y: x+ y**2,pos0,0))
    v0=math.sqrt(reduce(lambda x,y: x+ y**2,vel0,0))
    vr0=(reduce(lambda x,y: x+y, pos0*vel0,0))/r0
    alpha=2./r0-v0**2/mu
      
    if(alpha >= 0):
        xi0=smu*alpha*dt
    else:
        xi0=math.sign(dt)/math.sqrt(-alpha)*math.log(1-2*mu*dt*alpha/((vr0*r0)+  
         math.sign(dt)*smu/math.sqrt(-alpha)*(1-r0*alpha)) )
# this last formula is 4.5.11 in bate et al., fundamentals of astrodynamics 
# with +1 in the logarithm
        dxi0=smu/r0*dt
        if(abs(alpha*dxi0**2)<1):
            xi0=dxi0

    
    def f(xi):
        return universal_kepler(xi,r0,vr0,smu,alpha)-smu*dt
    def df(xi):
        return universal_kepler_dxi(xi,r0,vr0,smu,alpha)
    
    xi,err=newton(f,xi0,fprime=df,tol=1.e-10)    
#    print dt,xi,xi0
    
    pos=pos0*lagrange_f(xi,r0,vr0,smu,alpha)+vel0*lagrange_g(xi,r0,vr0,smu,alpha)
    r=math.sqrt(reduce(lambda x,y: x+ y**2,pos,0))
    vel=pos0*smu/r*lagrange_dfdxi(xi,r0,vr0,smu,alpha)+ \
        vel0*smu/r*lagrange_dgdxi(xi,r0,vr0,smu,alpha)
    return pos,vel


def universal_solver(mu,pos0,vel0,dt):
    smu=math.sqrt(mu)
    
    r0=math.sqrt(reduce(lambda x,y: x+ y**2,pos0,0))
    v0=math.sqrt(reduce(lambda x,y: x+ y**2,vel0,0))
    vr0=(reduce(lambda x,y: x+y, pos0*vel0,0))/r0
    alpha=2./r0-v0**2/mu
      
    if(alpha >= 0):
        xi0=smu*alpha*dt
    else:
        xi0=math.sign(dt)/math.sqrt(-alpha)*math.log(1-2*mu*dt*alpha/((vr0*r0)+  
         math.sign(dt)*smu/math.sqrt(-alpha)*(1-r0*alpha)) )
# this last formula is 4.5.11 in bate et al., fundamentals of astrodynamics 
# with +1 in the logarithm
        dxi0=smu/r0*dt
        if(abs(alpha*dxi0**2)<1):
            xi0=dxi0

    def f(xi):
        return universal_kepler(xi,r0,vr0,smu,alpha)-smu*dt
    def df(xi):
        return universal_kepler_dxi(xi,r0,vr0,smu,alpha)
    def ddf(xi):
        return universal_kepler_dxidxi(xi,r0,vr0,smu,alpha)
    
    xi,err=laguerre(f,xi0,fprime=df,fprimeprime=ddf)    
#    print dt,xi,xi0
    
    pos=pos0*lagrange_f(xi,r0,vr0,smu,alpha)+vel0*lagrange_g(xi,r0,vr0,smu,alpha)
    r=math.sqrt(reduce(lambda x,y: x+ y**2,pos,0))
    vel=pos0*smu/r*lagrange_dfdxi(xi,r0,vr0,smu,alpha)+ \
        vel0*smu/r*lagrange_dgdxi(xi,r0,vr0,smu,alpha)
    return pos,vel


class TwoBodyImplementation(object):
    __G = 1.0
    
    def __init__(self):
        self.initialize_code()
    
    def initialize_code(self):
        self.particles=[]
        self.tnow=0.0
        self.begin_time_parameter = 0.0
        return 0
          
    def cleanup_code(self):
        self.particles=[]
        return 0
      
    def commit_parameters(self):
        if self.tnow == 0.0:
            self.tnow = self.begin_time_parameter
            
        return 0
      
    def commit_particles(self):
        return 0
      
    def synchronize_model(self):
        return 0
      
    def initialize(self):
        pass
      
        
    def new_particle(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        index_of_the_particle.value = 0
        if( len(self.particles)>=2):
            return -1
        self.particles.append( 
            {
            'mass': mass, 
            'radius' : radius, 
            'x' : x, 
            'y' : y, 
            'z' : z,
            'vx' : vx,
            'vy' : vy,
            'vz' : vz,
            }
        )
        
        index_of_the_particle.value = len(self.particles)-1
        return 0
      
    def set_state(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        try:
            particle = self.particles[index_of_the_particle]
            particle['mass'] = mass  
            particle['radius'] = radius 
            particle['x'] =x
            particle['y'] =y
            particle['z'] =z
            particle['vx'] =vx
            particle['vy'] =vy
            particle['vz'] =vz
            return 0
        except:        
            return -1

    def set_mass(self, index_of_the_particle, mass):
        try:
            particle = self.particles[index_of_the_particle]
            particle['mass'] = mass  
            return 0
        except:        
            return -1

    def set_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            particle = self.particles[index_of_the_particle]
            particle['vx'] =vx
            particle['vy'] =vy
            particle['vz'] =vz
            return 0
        except:        
            return -1

    def set_position(self, index_of_the_particle, x, y, z):
        try:
            particle = self.particles[index_of_the_particle]
            particle['x'] =x
            particle['y'] =y
            particle['z'] =z
            return 0
        except:        
            return -1

      
    def get_state(self, index_of_the_particle, mass, radius, x, y, z, vx, vy, vz):
        try:
            particle = self.particles[index_of_the_particle]
            mass.value = particle['mass']
            radius.value = particle['radius']
            x.value = particle['x']
            y.value = particle['y']
            z.value = particle['z']
            vx.value = particle['vx']
            vy.value = particle['vy']
            vz.value = particle['vz']
            return 0
        except:        
            return -1

    def get_mass(self, index_of_the_particle, mass):
        try:
            particle = self.particles[index_of_the_particle]
            mass.value = particle['mass']
            return 0
        except:        
            return -1

    def get_radius(self, index_of_the_particle, radius):
        try:
            particle = self.particles[index_of_the_particle]
            radius.value = particle['radius']
            return 0
        except:        
            return -1

    def get_position(self, index_of_the_particle, x, y, z):
        try:
            particle = self.particles[index_of_the_particle]
            x.value = particle['x']
            y.value = particle['y']
            z.value = particle['z']
            return 0
        except:        
            return -1

    def get_velocity(self, index_of_the_particle, vx, vy, vz):
        try:
            particle = self.particles[index_of_the_particle]
            vx.value = particle['vx']
            vy.value = particle['vy']
            vz.value = particle['vz']
            return 0
        except:        
            return -1
      
    def get_time(self, time):
        time.value = self.tnow
        return 0
      
    def get_kinetic_energy(self, kinetic_energy):
        if(len(self.particles)!=1 and len(self.particles)!=2):
            return -1
          
        if(len(self.particles)==1):
            return -2
                      
        if(len(self.particles)==2):
            vel0=numpy.array( [self.particles[0]['vx'],  
                               self.particles[0]['vy'],  
                               self.particles[0]['vz']] )
            vel1=numpy.array( [self.particles[1]['vx'],  
                               self.particles[1]['vy'],  
                               self.particles[1]['vz']] )
            
            v0=reduce(lambda x,y: x+ y**2,vel0,0)
            v1=reduce(lambda x,y: x+ y**2,vel1,0)

            mass0=self.particles[0]['mass']
            mass1=self.particles[1]['mass']
        
            kinetic_energy.value = 0.5*(mass0*v0+mass1*v1)
            
            return 0 
  
    def get_gravity_at_point(self, eps, x, y, z, ax, ay, az, npoints):
        if(len(self.particles)==1):
            return -2
        elif(len(self.particles)==2):
            ax.value=0.
            ay.value=0.
            az.value=0.
            for i in [0,1]:
                mass=self.particles[i]['mass']
                xx=self.particles[i]['x']
                yy=self.particles[i]['y']
                zz=self.particles[i]['z']
                dr2=((xx-x)**2+(yy-y)**2+(zz-z)**2+eps**2)
                
                ax.value+=self.__G*mass*(xx-x)/dr2**1.5
                ay.value+=self.__G*mass*(yy-y)/dr2**1.5
                az.value+=self.__G*mass*(zz-z)/dr2**1.5
            
            return 0
        else:
            return -1
    
    def get_potential_at_point(self, eps, x, y, z, phi, npoints):
        if(len(self.particles)==1):
            return -2
        elif(len(self.particles)==2):
            phi.value = zero
            for i in [0,1]:
                mass=self.particles[i]['mass']
                xx=self.particles[i]['x']
                yy=self.particles[i]['y']
                zz=self.particles[i]['z']
                dr2=((xx-x)**2+(yy-y)**2+(zz-z)**2+eps**2)
                
                phi.value += -self.__G*mass/dr2**0.5
            return 0
        else:
            return -1
    
    
    def get_begin_time(self, value_out):
        value_out.value = self.begin_time_parameter
        return 0
        
    def set_begin_time(self, value_in):
        self.begin_time_parameter = value_in
        return 0
        
    def get_potential_energy(self, potential_energy):
        if(len(self.particles)!=1 and len(self.particles)!=2):
            return -1
          
        if(len(self.particles)==1):
            return -2    
        if(len(self.particles)==2):
            pos0=numpy.array( [self.particles[0]['x'],  
                               self.particles[0]['y'],  
                               self.particles[0]['z']] )
            pos1=numpy.array( [self.particles[1]['x'],  
                               self.particles[1]['y'],  
                               self.particles[1]['z']] )
            dpos=pos0-pos1
            mass0=self.particles[0]['mass']
            mass1=self.particles[1]['mass']
            r=math.sqrt(reduce(lambda x,y: x+ y**2,dpos,0))
            potential_energy.value = -self.__G*mass0*mass1/r
            return 0    

    def evolve_model(self, time):
        time_end = time
        
        if(len(self.particles)!=1 and len(self.particles)!=2):
            return -1
          
        if(len(self.particles)==1):
            mu=self.__G*self.particles[0]['mass']
            radius=self.particles[0]['radius']        
            dpos_initial=numpy.array( [self.particles[0]['x'],  
                             self.particles[0]['y'],  
                             self.particles[0]['z']] )
            dvel_initial=numpy.array( [self.particles[0]['vx'],  
                             self.particles[0]['vy'],  
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
            pos0=numpy.array( [self.particles[0]['x'],  
                               self.particles[0]['y'],  
                               self.particles[0]['z']] )
            vel0=numpy.array( [self.particles[0]['vx'],  
                               self.particles[0]['vy'],  
                               self.particles[0]['vz']] )
            pos1=numpy.array( [self.particles[1]['x'],  
                               self.particles[1]['y'],  
                               self.particles[1]['z']] )
            vel1=numpy.array( [self.particles[1]['vx'],  
                               self.particles[1]['vy'],  
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



class TwoBodyInterface(PythonCodeInterface, GravitationalDynamicsInterface, GravityFieldInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, TwoBodyImplementation, 'twobody_worker', **options)
    
   

class TwoBody(GravitationalDynamics, GravityFieldCode):
    
    
    def __init__(self, convert_nbody = None, **options):
        nbody_interface = TwoBodyInterface(**options)
        
        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )     
        
    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
    
    def get_epsilon_squared(self):
        return zero
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_epsilon_squared",
            None, 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = zero
        )
        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
    
