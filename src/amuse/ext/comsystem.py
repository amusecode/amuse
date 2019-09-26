"""
function center_of_mass_system generates a thin wrapper class
to a grav dynamics class to calculate the dynamics in the center of mass frame
"""

from amuse.datamodel import Particles

def center_of_mass_system(baseclass):

  class comsystem(baseclass):
     def __init__(self,*args,**kwargs):
       self.particles_accessed=True
       self._particles=Particles(0)
       baseclass.__init__(self,*args,**kwargs)
     
     def evolve_model(self,*args,**kwargs):
       if self.particles_accessed:         
         self.com_position=self._particles.center_of_mass()
         self.com_velocity=self._particles.center_of_mass_velocity()
         com_time=self.model_time
         self._particles.synchronize_to(self.overridden().particles)
         self._particles.new_channel_to(self.overridden().particles).copy_attributes(["mass","radius"])
         self.overridden().particles.position=self._particles.position-self.com_position
         self.overridden().particles.velocity=self._particles.velocity-self.com_velocity
       self.overridden().evolve_model(*args,**kwargs)
       self.com_position+=self.com_velocity*(self.model_time-com_time)        
       self.particles_accessed=False
       
     @property
     def particles(self):
       if not self.particles_accessed:
         self._particles=self.overridden().particles.copy()
         self._particles.position+=self.com_position
         self._particles.velocity+=self.com_velocity
         self.particles_accessed=True
       return self._particles

  return comsystem


if __name__=="__main__":
  from amuse.community.huayno.interface import Huayno
  from amuse.ic.plummer import new_plummer_model
  from amuse.units import nbody_system

  comHuayno=center_of_mass_system(Huayno)
  grav=comHuayno()

  parts=new_plummer_model(100)
  
  grav.particles.add_particles(parts)
  
  grav.evolve_model(1| nbody_system.time)
  
  print(grav.particles)
  
  
  
  
  
  

