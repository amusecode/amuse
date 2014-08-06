"""
Wrapper of the bridge-like integrator for amuse. To see how bridge works go to the file bridge.py in this folder.

This Rotating bridge makes the integration of eqs of motion in a counterclockwise rotating system.

Usage:
from amuse.ext.rotating_bridge import Rotating_Bridge


"""

from amuse.support.exceptions import AmuseException
import threading
from amuse.units import quantities
from amuse.units import units
from amuse.units import quantities
from amuse import datamodel
from amuse.ext.bridge import bridge
import numpy



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

        if self.omega.value_in(self.omega.unit) !=0.:
            vx0=parts.vx
            vy0=parts.vy
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
