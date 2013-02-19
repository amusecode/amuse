"""
  bridge-like integrator for amuse
  
  the bridge class provides a bridge like coupling between different 
  gravitational integrators. In this way a system composed of multiple 
  components can be evolved taking account of the self gravity of the whole 
  system self consistently, while choosing the most appropiate integrator
  for the self-gravity of the component systems. This is mainly useful for
  systems  consist of two or more components that are either well separated
  spatially or have different scales (otherwise using a single integrator is
  more efficient) 

  The main idea is that systems experience each others gravity through 
  periodic velocty kicks with ordinary evolution in  between - the evolution
  is thus described by an alternation of drift (D) and  kick (K) operators,
  here chosen as:

       K(1/2 dt) D(dt) K(1/2 dt)    
  
  K(dt) denotes a kick of the velocities over a timestep dt, while D(dt)
  denotes  a drift, meaning secular evolution using self gravity of the
  system, over dt.

  implementation notes:
  
  In order to use bridge the component systems should be initialized as usual,
  then a bridge systems is initialized, after which one or more systems are
  added:
  
  from amuse.ext.bridge import bridge
   
  bridgesys=bridge(verbose=False)
  bridgesys.add_system(galaxy, (cluster,), False)
  bridgesys.add_system(cluster, (galaxy,), True )

  bridge builds on the full gravity interface, so unit handling etc is 
  guaranteed. Bridge itself is a (somewhat incomplete) gravity interface,
  so  the usual evolve, get_potential methods work (and bridge can be a
  component  in a bridge systems). Note that a single coordinate system should
  be used at the moment for all the components systems (different units are 
  allowed though). The call to add systems, for example:

  bridgesys.add_system(galaxy, False, (cluster,))
  
  has three arguments: the system, a flag to specify whether
  synchronization  is needed and a set with *interaction* partners. The
  interaction partners indicate which systems will kick the system. In the
  most simple case these  would be the set of other systems that are added,
  but usually this is not  what you want to get good performace. In some
  cases you want to ignore one  direction of interaction (eg. in a combined
  simulation of a galaxy and a  comet orbits around a star you may want the
  ignore the gravity of the comet), in other cases you want to use a
  different force calculator (eg integrating a cluster in  a galaxy where
  the galaxy is evolved with a tree code and the cluster with a direct sum
  code, one also would want to use a tree code to calculate the cluster
  gravity for the galaxy. In such a case one can derive a skeleton gravity
  interface from  the cluster system.  A module is provided with some
  examples of such *derived* systems, derived_grav_systems.py 

  Hints for good use:
  
  The bridgesys is flexible but care should be taken in order to obtain 
  valid results. For one thing, there is no restriction or check on the 
  validity of the assumption of well seperated dynamics: for example any 
  system could be split up and put together in bridge, but if the timestep
  is chosen to be larger than the timestep criterion of the code, the
  integration will show errors. 
  
  For good performance one should use derived systems to reduce the
  complexity where possible. 
  
  There is an issue with the synchronization: some codes do not end on the
  exact time of an evolve, or need an explicit sync call. In these cases it
  is up to the user to  determine whether bridge can be used (an explicit
  sync call may induce extra errors that degrade the order of the
  integrator).

"""  


# issues:
# - for now, units in si 
# - a common coordinate system is used for all systems
# - sync of systems should be checked
# - timestepping: adaptive dt?

import threading

from amuse.units import quantities
from amuse.units import units

from amuse import datamodel
def potential_energy(system, get_potential):
    parts=system.particles.copy()
    pot=get_potential(parts.radius,parts.x,parts.y,parts.z)
    return (pot*parts.mass).sum() / 2

def kick_system(system, get_gravity, dt):
    parts=system.particles.copy()
    ax,ay,az=get_gravity(parts.radius,parts.x,parts.y,parts.z)
    parts.vx=parts.vx+dt*ax
    parts.vy=parts.vy+dt*ay
    parts.vz=parts.vz+dt*az
    channel=parts.new_channel_to(system.particles)
    channel.copy_attributes(["vx","vy","vz"])   
#    parts.copy_values_of_all_attributes_to(system.particles)
  
class bridge(object):
    def __init__(self,verbose=False,method=None):
        """
        verbose indicates whether to output some run info
        """  
        self.systems=set()
        self.partners=dict()
        self.time_offsets=dict()
        self.time=quantities.zero
        self.do_sync=dict()
        self.verbose=verbose
        self.timestep=None
        self.method=method
    
    def add_system(self, interface,  partners=set(),do_sync=True):
        """
        add a system to bridge integrator  
        """
        if hasattr(interface,"model_time"):
            self.time_offsets[interface]=(self.time-interface.model_time)
        else:
            self.time_offsets[interface]=quantities.zero     
        self.systems.add(interface)
        for p in partners:
            if not hasattr(p,"get_gravity_at_point"):
                return -1
        self.partners[interface]=partners
        self.do_sync[interface]=do_sync  
        return 0
      
    def evolve_model(self,tend,timestep=None):
        """
        evolve combined system to tend, timestep fixes timestep
        """
        if timestep is None:
            if self.timestep is None:
                timestep=tend-self.time
            else:
                timestep = self.timestep
                
        if self.method==None:
          return self.evolve_joined_leapfrog(tend,timestep)
        else:
          return self.evolve_simple_steps(tend,timestep)          

    def evolve_simple_steps(self,tend,timestep):
        while self.time < (tend-timestep/2):
            self._drift_time=self.time
            self.method(self.kick_systems,self.drift_systems_dt, timestep)
            self.time=self.time+timestep
        return 0    

    def evolve_joined_leapfrog(self,tend,timestep):
        first=True
        while self.time < (tend-timestep/2):
             if first:      
                 self.kick_systems(timestep/2)
                 first=False
             else:
                 self.kick_systems(timestep)
             self.drift_systems(self.time+timestep)
             self.time=self.time+timestep
        if not first:
             self.kick_systems(timestep/2)         
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
    def model_time(self):  
         return self.time
      
    @property
    def potential_energy(self):
        Ep=quantities.zero
        for x in self.systems:
            Ep+=x.potential_energy
            if hasattr(x,"particles"):
                for y in self.partners[x]:
                    Ep += potential_energy(x,y.get_potential_at_point)
        return Ep
    
    @property
    def kinetic_energy(self):  
        Ek=quantities.zero
        for x in self.systems:
            Ek+=x.kinetic_energy
        return Ek
        
    @property
    def thermal_energy(self):  
        result=quantities.zero
        for x in self.systems:
            if hasattr(x,'thermal_energy'):
                result+=x.thermal_energy
        return result
          
    @property
    def particles(self):
        arr=[]
        for x in self.systems:
            if hasattr(x,"particles"):
                arr.append(x.particles)
        return datamodel.ParticlesSuperset(arr)                

# 'private' functions
    def drift_systems_dt(self,dt):
        self._drift_time+=dt
        self.drift_systems(self._drift_time)

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
