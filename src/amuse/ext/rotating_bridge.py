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
        timestep=kwargs.pop('timestep')
        self.omega= omega        
        bridge.__init__(self, **kwargs)
        self.timestep=timestep        

    def kick_system_rotational(self, system, partners, dt):
        parts=system.particles.copy()
        ax= quantities.zero
        ay= quantities.zero
        az= quantities.zero
        if(self.verbose):  
            print system.__class__.__name__,"receives kick from",
        for y in partners:
            if system is not y:
                if(self.verbose):  
                    print y.__class__.__name__,
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
            print ".. done"
            
    def kick_systems(self,dt):

        for x in self.systems:
            if self.do_sync[x]:
                if hasattr(x,"synchronize_model"):
                    if(self.verbose): print x.__class__.__name__,"is synchronizing",
                    x.synchronize_model()    
                    if(self.verbose):  print ".. done"
        for x in self.systems:
            if hasattr(x,"particles"):
                        self.kick_system_rotational(x, self.partners[x], dt)
        return 0

    
    @property
    def jacobi_potential_energy(self):
        parts=self.particles
        return -0.5*(parts.mass*self.omega**2*(parts.x**2+parts.y**2)).sum()
