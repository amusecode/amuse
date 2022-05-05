"""
Wrapper of the bridge-like integrator for amuse. To see how bridge works go to the file bridge.py in this folder.

This bridge makes the integration of eqs of motion in a right hand counterclockwise rotating system.

Usage:
from amuse.ext.composition_methods import *
from amuse.ext.rotating_bridge import Rotating_Bridge

system= Rotating_Bridge(omega, timestep= dt_bridge, verbose= False, method= method)
system.add_system(cluster, (MW,), False)
system.add_system(MW, (), False) 

omega: angular velocity of the rotating frame. If an axisymmetric model is used, omega must be 0 kms/kpc
dt_bridge: bridge timestep. By now it is fixed
method: One of the composite methods. The default one is LEAPFROG but it can be used:

SPLIT_4TH_S_M6
SPLIT_4TH_S_M5
SPLIT_4TH_S_M4
SPLIT_6TH_SS_M11
SPLIT_6TH_SS_M13
SPLIT_8TH_SS_M21
SPLIT_10TH_SS_M35 

"""

from amuse.support.exceptions import AmuseException
import threading
from amuse.units import quantities
from amuse import datamodel
from amuse.ext.bridge import bridge
import numpy

from numpy import cos,sin

from amuse.datamodel import TransformedParticles

# same as below, retained for legacy
def inertial_to_rotating(t,omega,parts):
  x=parts.x
  y=parts.y
  vx=parts.vx
  vy=parts.vy
  rotating=parts.copy()
  rotating.x=x*cos(omega*t)+y*sin(omega*t)
  rotating.y=-x*sin(omega*t)+y*cos(omega*t)
  rotating.vx=(vx+y*omega)*cos(omega*t)+(vy-x*omega)*sin(omega*t)
  rotating.vy=-(vx+y*omega)*sin(omega*t)+(vy-x*omega)*cos(omega*t)
  return rotating
  
def rotating_to_inertial(t,omega,parts):     
  return inertial_to_rotating(t,-omega,parts)

class Rotating_Bridge(bridge):
    def __init__(self, omega, **kwargs):
        timestep=kwargs.pop('timestep', None)
        self.omega=omega        
        self.initial_angle=kwargs.pop('initial_angle', 0.)
        bridge.__init__(self, **kwargs)
        self.timestep=timestep

    def kick_system_rotational(self, system, partners, dt):
        parts=system.particles.copy()
        ax= quantities.zero
        ay= quantities.zero
        az= quantities.zero
        if(self.verbose):  
            print(system.__class__.__name__,"receives kick from", end=' ')
        for y in partners:
            if system is not y:
                if(self.verbose):  
                    print(y.__class__.__name__, end=' ')
                _ax,_ay,_az= y.get_gravity_at_point(parts.radius,parts.x,parts.y,parts.z)
                ax+=_ax
                ay+=_ay
                az+=_az

        if self.omega != quantities.zero:
            vx0=parts.vx.copy()
            vy0=parts.vy.copy()
            omega=2*self.omega
            a1_omega=(ax+self.omega**2*parts.x)/omega
            a2_omega=(ay+self.omega**2*parts.y)/omega
            parts.vx=(vx0-a2_omega)*numpy.cos(omega*dt)+(vy0+a1_omega)*numpy.sin(omega*dt)+a2_omega
            parts.vy=-(vx0-a2_omega)*numpy.sin(omega*dt)+(vy0+a1_omega)*numpy.cos(omega*dt)-a1_omega
            parts.vz=parts.vz+az*dt
        else:
            parts.vx=parts.vx+ax*dt
            parts.vy=parts.vy+ay*dt
            parts.vz=parts.vz+az*dt
       
        channel=parts.new_channel_to(system.particles)
        channel.copy_attributes(["vx","vy","vz"])   
        if(self.verbose):
            print(".. done")
            
    def kick_systems(self,dt):

        for x in self.systems:
            if self.do_sync[x]:
                if hasattr(x,"synchronize_model"):
                    if(self.verbose): print(x.__class__.__name__,"is synchronizing", end=' ')
                    x.synchronize_model()    
                    if(self.verbose):  print(".. done")
        for x in self.systems:
            if hasattr(x,"particles"):
                        self.kick_system_rotational(x, self.partners[x], dt)
        return 0

    
    @property
    def jacobi_potential_energy(self):
        parts=self.particles
        return -0.5*(parts.mass*self.omega**2*(parts.x**2+parts.y**2)).sum()

    def transform_inertial_to_rotating(self,x,y,vx,vy,inverse=False):
        angle = self.initial_angle + self.omega*self.model_time
        omega = self.omega
        if inverse:
            angle=-angle
            omega=-omega

        C1 = vx + omega*y
        C2 = vy - omega*x
        x_ =  x * numpy.cos(angle) + y * numpy.sin(angle)
        y_ = -x * numpy.sin(angle) + y * numpy.cos(angle)
        vx_ = C1*numpy.cos(angle) + C2*numpy.sin(angle)
        vy_ = C2*numpy.cos(angle) - C1*numpy.sin(angle)
        return x_,y_,vx_,vy_

    def transform_rotating_to_inertial(self,x,y,vx,vy):
        return self.transform_inertial_to_rotating(x,y,vx,vy, inverse=True)

    # this return a view on self.particles
    # which automatically updates 
    # (uses above transforms for this reason)
    @property
    def particles_inertial_frame(self):
        return TransformedParticles(self.particles, 
                                    ["x","y","vx","vy"], 
                                    self.transform_rotating_to_inertial,
                                    ["x","y","vx","vy"],
                                    self.transform_inertial_to_rotating,
                                    )

    @property
    def gas_particles_inertial_frame(self):
        return TransformedParticles(self.gas_particles, 
                                    ["x","y","vx","vy"], 
                                    self.transform_rotating_to_inertial,
                                    ["x","y","vx","vy"],
                                    self.transform_inertial_to_rotating,
                                    )
  
class RotatingBridgeInertialParticles(Rotating_Bridge):
    """
    same as above, except non-inertial frame is hidden. Note that:
    code.particles.get_subsets()[i]
    gets a view on the ith code particles
    
    """
    
    @property
    def particles(self):
        arr=[]
        for x in self.systems:
            if hasattr(x,"particles"):
                arr.append(x.particles)
        particles=datamodel.ParticlesSuperset(arr) 
        return TransformedParticles(particles, 
                                    ["x","y","vx","vy"], 
                                    self.transform_rotating_to_inertial,
                                    ["x","y","vx","vy"],
                                    self.transform_inertial_to_rotating,
                                    )  

    @property
    def gas_particles(self):
        arr=[]
        for x in self.systems:
            if hasattr(x,"gas_particles"):
                arr.append(x.gas_particles)
        particles=datamodel.ParticlesSuperset(arr) 

        return TransformedParticles(particles, 
                                    ["x","y","vx","vy"], 
                                    self.transform_rotating_to_inertial,
                                    ["x","y","vx","vy"],
                                    self.transform_inertial_to_rotating,
                                    )
