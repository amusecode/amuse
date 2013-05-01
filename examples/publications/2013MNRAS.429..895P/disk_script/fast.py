# issues:
# - for now, units in si 
# - a common coordinate system is used for all systems
# - sync of systems should be checked
# - timestepping: adaptive dt?
from amuse.units import units
from amuse.units.quantities import zero
import threading

from amuse import datamodel
def radius_or_hsmooth(parts):
  d=set(dir(parts))
  if "radius" in d:
    return parts.radius
  else:
    if "h_smooth" in d:
      return parts.h_smooth
    else:
      print d
      raise Exception    


def potential_energy(system, get_potential):
  parts=system.particles.copy()
  pot=get_potential(parts.radius,parts.x,parts.y,parts.z)
  return (pot*parts.mass).sum()/2

def kick_system(system, get_gravity, dt):
  parts=system.particles.copy()
  ax,ay,az=get_gravity(parts.radius,parts.x,parts.y,parts.z)
  parts.vx=parts.vx+dt*ax
  parts.vy=parts.vy+dt*ay
  parts.vz=parts.vz+dt*az
#  parts.copy_values_of_state_attributes_to(system.particles)
  channel = parts.new_channel_to(system.particles)
  channel.copy_attributes(["vx","vy","vz"])
  
class FAST(object):
  def __init__(self,verbose=False,timestep=None):
    """
    verbose indicates whether to output some run info
    """  
    self.systems=set()
    self.partners=dict()
    self.time_offsets=dict()
    self.time=0. | units.s
    self.do_sync=dict()
    self.verbose=verbose
    self.timestep=timestep
  
  def set_timestep(self,timestep):
    self.timestep=timestep
  
  def add_system(self, interface,  partners=set(),do_sync=True):
    """
    add a system to bridge integrator  
    """
    if hasattr(interface,"model_time"):
      self.time_offsets[interface]=(self.time-interface.model_time)
    else:
      self.time_offsets[interface]=0.        
    self.systems.add(interface)
    for p in partners:
      if not hasattr(interface,"get_gravity_at_point"):
          return -1
    self.partners[interface]=partners
    self.do_sync[interface]=do_sync  
    return 0
    
  def evolve_model(self,tend):
    """
    evolve combined system to tend, timestep fixes timestep
    """
    timestep=self.timestep
    if timestep is None:
      timestep=tend-self.time
    first=True
    while self.time < (tend-timestep/2.):    
      if first:      
        self.kick_systems(timestep/2.)
        first=False
      else:
        self.kick_systems(timestep)
      self.drift_systems(self.time+timestep)
      self.time=self.time+timestep
    if not first:
      self.kick_systems(timestep/2.)         
    return 0    
  
  def synchronize_model(self):
    """ 
    explicitly synchronize all components
    """
    for x in self.systems:
      if hasattr(x,"synchronize_model"):
        if(self.verbose): print x.__class__.__name__,"is synchronizing",
        x.synchronize_model()    
        if(self.verbose): print ".. done"
                          
  def get_potential_at_point(self,radius,x,y,z):
    err=0
    pot=0.*radius
    for x in self.systems:
      _pot,err=x.get_potential_at_point(radius,x,y,z)
      if err != 0: 
        break
      pot=pot+_pot
    return pot,err
      
  def get_gravity_at_point(self,radius,x,y,z):
    err=0
    ax=0.*radius
    ay=0.*radius
    az=0.*radius
    for x in self.systems:
      _ax,_ay,_az,err=x.get_gravity_at_point(radius,x,y,z)
      if err != 0: 
        break
      ax=ax+_ax
      ay=ay+_ay
      az=az+_az
    return ax,ay,az,err
    
  @property
  def potential_energy(self):
    Ep=zero
    for x in self.systems:
      Ep+=x.potential_energy
      if hasattr(x,"particles"):
        for y in self.partners[x]:
          Ep+=potential_energy(x,y.get_potential_at_point)
    return Ep
  
  @property
  def kinetic_energy(self):  
    Ek=zero
    for x in self.systems:
      Ek+=x.kinetic_energy
    return Ek

  @property
  def thermal_energy(self):  
    Eth=zero
    for x in self.systems:
      if hasattr(x,'thermal_energy'):
        Eth+=x.thermal_energy
    return Eth

  @property
  def model_time(self):  
    return self.time
        
  @property
  def particles(self):
    arr=[]
    for x in self.systems:
      if hasattr(x,"particles"):
        arr.append(x.particles)
    return datamodel.ParticlesSuperset(arr)          

  @property
  def gas_particles(self):
    arr=[]
    for x in self.systems:
      if hasattr(x,"gas_particles"):
        arr.append(x.gas_particles)
    return datamodel.ParticlesSuperset(arr)          

# 'private' functions

  def drift_systems(self,tend):
    threads=[]
    for x in self.systems:
      if hasattr(x,"evolve_model"):
        offset=self.time_offsets[x]
        if(self.verbose):
          print "evolving", x.__class__.__name__,
        threads.append(threading.Thread(target=x.evolve_model, args=(tend-offset,)) )
    for x in threads:
      x.start()
    for x in threads:
      x.join()
    if(self.verbose): 
      print ".. done"
    return 0

  def kick_systems(self,dt):
    for x in self.systems:
      if self.do_sync[x]:
        if hasattr(x,"synchronize_model"):
          if(self.verbose): print x.__class__.__name__,"is synchronizing",
          x.synchronize_model()    
          if(self.verbose):  print ".. done"
    for x in self.systems:
      if hasattr(x,"particles"):
        for y in self.partners[x]:
          if x is not y:
            if(self.verbose):  print x.__class__.__name__,"receives kick from",y.__class__.__name__,
            kick_system(x,y.get_gravity_at_point,dt)
            if(self.verbose):  print ".. done"
    return 0
