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

  bridgesys.add_system(galaxy, (cluster,), False)
  
  has three arguments: the system, a set with *interaction* partners and
  a flag to specify whether synchronization is needed . The
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
from amuse.units import units, constants, generic_unit_system, nbody_system
from amuse import datamodel
from amuse.support.exceptions import AmuseException

  

class AbstractCalculateFieldForCodes(object):
    """
    Calculated gravity and potential fields using the particles
    of other codes with the code provided.
    """
    
    def __init__(self, input_codes, verbose=False):
        """
        verbose indicates whether to output some run info
        """  
        self.codes_to_calculate_field_for = input_codes
        self.verbose=verbose
      
    def evolve_model(self,tend,timestep=None):
        """
        """
        
    def get_potential_at_point(self,radius,x,y,z):
        code = self._setup_code()
        try:
            for input_code in self.codes_to_calculate_field_for:
                code.particles.add_particles(input_code.particles)
            code.commit_particles()
            return code.get_potential_at_point(radius,x,y,z)
        finally:
            self._cleanup_code(code)
        
    def get_gravity_at_point(self,radius,x,y,z):
        code = self._setup_code()
        try:
            for input_code in self.codes_to_calculate_field_for:
                code.particles.add_particles(input_code.particles)
            code.commit_particles()
            return code.get_gravity_at_point(radius,x,y,z)
        finally:
            self._cleanup_code(code)
    
    def _setup_code(self):
        pass
        
    def _cleanup_code(self, code):
        pass
        

class CalculateFieldForCodes(AbstractCalculateFieldForCodes):
    """
    Calculated gravity and potential fields using the particles
    of other codes with the code provided. 
    The code is created for every calculation.
    """
    
    def __init__(self, code_factory_function, input_codes, verbose=False):
        AbstractCalculateFieldForCodes.__init__(self, input_codes, verbose)
        self.code_factory_function = code_factory_function
      
    def _setup_code(self):
        return self.code_factory_function()
        
    def _cleanup_code(self, code):
        code.stop()

class CalculateFieldForCodesUsingReinitialize(AbstractCalculateFieldForCodes):
    """
    Calculated gravity and potential fields using the particles
    of other codes with the code provided. 
    The code is created for every calculation.
    """
    
    def __init__(self, code, input_codes, verbose=False):
        AbstractCalculateFieldForCodes.__init__(self, input_codes, verbose)
        self.code = code
      
    def _setup_code(self):
        return self.code
        
    def _cleanup_code(self, code):
        code.reset()
        
class CalculateFieldForCodesUsingRemove(AbstractCalculateFieldForCodes):
    """
    Calculated gravity and potential fields using the particles
    of other codes with the code provided. 
    The code is created for every calculation.
    """
    
    def __init__(self, code, input_codes, verbose=False):
        AbstractCalculateFieldForCodes.__init__(self, input_codes, verbose)
        self.code = code
      
    def _setup_code(self):
        return self.code
        
    def _cleanup_code(self, code):
        code.particles.remove_particles(code.particles)


class CalculateFieldForParticles(object):
    """
    Calculates an field for a set of particles, the set
    of particles can be from another code.
    """
    
    def __init__(self, particles = None, gravity_constant = None,
            softening_mode="shared"):
        if particles is None:
            self.particles=datamodel.Particles()      
        else:
            self.particles = particles
        if gravity_constant is None:
            if len(particles) and hasattr(particles, 'mass'):
                try:
                    particles[0].mass.value_in(units.kg)
                    self.gravity_constant = constants.G
                except:
                    raise AmuseException("For generic units the gravity_constant must be specified")
            else:
                raise AmuseException("Particle data not yet available, so the gravity_constant must be specified")
        else:
            self.gravity_constant = gravity_constant
        
        if softening_mode == "individual" or softening_mode == "radius":
            self._softening_lengths_squared = self._softening_lengths_squared_individual
        elif softening_mode == "h_smooth":
            self._softening_lengths_squared = self._softening_lengths_squared_h_smooth
        else:
            self._softening_lengths_squared = self._softening_lengths_squared_shared
        self.smoothing_length_squared = quantities.zero
    
    def _softening_lengths_squared_individual(self):
        return self.particles.radius**2
    def _softening_lengths_squared_h_smooth(self):
        return self.particles.h_smooth**2
    def _softening_lengths_squared_shared(self):
        return self.smoothing_length_squared#.as_vector_with_length(len(self.particles))
    
    def cleanup_code(self):
        self.particles = datamodel.Particles()
        
    def evolve_model(self,tend,timestep=None):
        """
        """
        
    def get_potential_at_point(self,radius,x,y,z):
        positions = self.particles.position
        result = quantities.AdaptingVectorQuantity()
        
        for i in range(len(x)):
            dx = x[i] - positions.x
            dy = y[i] - positions.y
            dz = z[i] - positions.z
            dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
            dr = (dr_squared + self._softening_lengths_squared()).sqrt()
            energy_of_this_particle = (self.particles.mass / dr).sum()
            result.append(-self.gravity_constant * energy_of_this_particle)
        return result
    
    
    def get_gravity_at_point(self,radius,x,y,z):
        positions = self.particles.position
        m1 = self.particles.mass
        result_ax = quantities.AdaptingVectorQuantity()
        result_ay = quantities.AdaptingVectorQuantity()
        result_az = quantities.AdaptingVectorQuantity()
        for i in range(len(x)):
            dx = x[i] - positions.x
            dy = y[i] - positions.y
            dz = z[i] - positions.z
            dr_squared = ((dx * dx) + (dy * dy) + (dz * dz) + 
                self._softening_lengths_squared() + radius[i]**2)
            
            ax = -self.gravity_constant * (m1*dx/dr_squared**1.5).sum()
            ay = -self.gravity_constant * (m1*dy/dr_squared**1.5).sum()
            az = -self.gravity_constant * (m1*dz/dr_squared**1.5).sum()
            
            result_ax.append(ax)
            result_ay.append(ay)
            result_az.append(az)
        return result_ax, result_ay, result_az
        

class GravityCodeInField(object):
    
    
    def __init__(self, code, field_codes, do_sync=True, verbose=False, radius_is_eps=False, h_smooth_is_eps=False):
        """
        verbose indicates whether to output some run info
        """  
        self.code = code
        self.field_codes = field_codes
        
        if hasattr(self.code, 'model_time'):
            self.time = self.code.model_time
        else:
            self.time = quantities.zero
            
        self.do_sync=do_sync
        self.verbose=verbose
        self.timestep=None
        self.radius_is_eps = radius_is_eps
        self.h_smooth_is_eps = h_smooth_is_eps
        if not hasattr(self.code,"parameters"):
            self.zero_smoothing=True
        elif not hasattr(self.code.parameters,"epsilon_squared"):
            self.zero_smoothing=True
        else:
            self.zero_smoothing=False    
    
      
    def evolve_model(self,tend,timestep=None):
        """
        evolve combined system to tend, timestep fixes timestep
        """
 
        if timestep is None:
            timestep = self.timestep
 
        first=True
        while self.time < (tend-timestep/2.):    
            if first:      
                self.kick(timestep/2.)
                first=False
            else:
                self.kick(timestep)
             
            self.drift(self.time+timestep)

            self.time+=timestep
             
        if not first:
             self.kick(timestep/2.)         
        
 
    
    def synchronize_model(self):
        """ 
        explicitly synchronize all components
        """
        if hasattr(self.code,"synchronize_model"):
            if(self.verbose):
                print self.code.__class__.__name__,"is synchronizing",
            
            self.code.synchronize_model()    
            
            if(self.verbose):
                print ".. done"
                            
    def get_potential_at_point(self,radius,x,y,z):
        return self.code.get_potential_at_point(radius,x,y,z)
        
    def get_gravity_at_point(self,radius,x,y,z):
        return self.code.get_gravity_at_point(radius,x,y,z)

    @property
    def model_time(self):  
         return self.time
      
    @property
    def potential_energy(self):
        if not hasattr(self.code, 'particles'):
            return quantities.zero
            
        result = self.code.potential_energy
        particles=self.code.particles.copy()
        
        for y in self.field_codes:
            energy = self.get_potential_energy_in_field_code(particles, y)
            result += energy
        return result
    
    @property
    def kinetic_energy(self):  
        return self.code.kinetic_energy
        
    @property
    def thermal_energy(self):  
        if hasattr(self.code,'thermal_energy'):
            return self.code.thermal_energy
        else:
            return quantities.zero
          
    @property
    def particles(self):
        return self.code.particles         
    
    @property
    def gas_particles(self):
        if hasattr(self.code, "gas_particles"):
            return self.code.gas_particles
        else:
            raise AttributeError
    
    @property
    def dm_particles(self):
        if hasattr(self.code, "dm_particles"):
            return self.code.dm_particles
        else:
            raise AttributeError
    
    def drift(self, tend):
        if not hasattr(self.code,"evolve_model"):
            return
        if (self.verbose):
            print self.code.__class__.__name__, "is evolving to", tend
        
        self.code.evolve_model(tend)
        
        
        if(self.verbose): 
            print ".. done"

    def kick(self, dt):
        
        if not hasattr(self.code, 'particles') or len(self.code.particles)==0:
            return quantities.zero
        
        particles=self.code.particles.copy()
        kinetic_energy_before = particles.kinetic_energy()
        
        for field_code in self.field_codes:
            if(self.verbose):
                print self.code.__class__.__name__,"receives kick from",field_code.__class__.__name__,
            
            self.kick_with_field_code(
                particles,
                field_code,
                dt
            )
            
            if(self.verbose):
                print ".. done"
        
        channel=particles.new_channel_to(self.code.particles)
        channel.copy_attributes(["vx","vy","vz"])   
        
        kinetic_energy_after = particles.kinetic_energy()
        return kinetic_energy_after - kinetic_energy_before
        
    
    def _softening_lengths(self, particles):
        if self.radius_is_eps:
            return particles.radius
        elif self.h_smooth_is_eps:
            return particles.h_smooth
        elif self.zero_smoothing:
            return 0.*particles.x
        else:
            return (self.code.parameters.epsilon_squared**0.5).as_vector_with_length(len(particles))
    
    def get_potential_energy_in_field_code(self, particles, field_code):
        pot=field_code.get_potential_at_point(
            self._softening_lengths(particles),
            particles.x,
            particles.y,
            particles.z
        )
        return (pot*particles.mass).sum() / 2

    def kick_with_field_code(self, particles, field_code, dt):
        ax,ay,az=field_code.get_gravity_at_point(
            self._softening_lengths(particles),
            particles.x,
            particles.y,
            particles.z
        )
        self.update_velocities(particles, dt, ax, ay, az)
        
    def update_velocities(self,particles, dt,  ax, ay, az):
        particles.vx += dt * ax
        particles.vy += dt * ay
        particles.vz += dt * az 
  
class Bridge(object):
    def __init__(self, timestep = None, verbose=False, use_threading=True,method=None):
        """
        verbose indicates whether to output some run info
        """  
        self.codes=[]
        self.time=quantities.zero
        self.verbose=verbose
        self.timestep=timestep
        self.kick_energy = quantities.zero
        self.use_threading = use_threading
        self.time_offsets = dict()
        self.method=method
    
    def add_system(self, interface, partners=set(),do_sync=True, 
            radius_is_eps=False, h_smooth_is_eps=False):
        """
        add a system to bridge integrator  
        """
        
        if hasattr(interface, "particles"):
            code = GravityCodeInField(interface, partners, do_sync, self.verbose, 
                radius_is_eps, h_smooth_is_eps)
            self.add_code(code)
        else:
            if len(partners):
                raise Exception("You added a code without particles, but with partners, this is not supported!")
            self.add_code(interface)
    
    def add_code(self, code):
        self.codes.append(code)
        if hasattr(code,"model_time"):
            self.time_offsets[code]=(self.time-code.model_time)
        else:
            self.time_offsets[code]=quantities.zero   
         
      
    def evolve_model(self, tend, timestep=None):
        """
        evolve combined system to tend, timestep fixes timestep
        """
        if timestep is None:
            if self.timestep is None:
                timestep=tend-self.time
            else:
                timestep = self.timestep

        if self.method is None:
          return self.evolve_joined_leapfrog(tend,timestep)
        else:
          return self.evolve_simple_steps(tend,timestep)          

    def evolve_simple_steps(self,tend,timestep):
        while self.time < (tend-timestep/2):
            self._drift_time=self.time
            self.method(self.kick_codes,self.drift_codes_dt, timestep)
            self.time=self.time+timestep

    def evolve_joined_leapfrog(self,tend,timestep):
        first=True            
        while self.time < (tend-timestep/2.):    
             if first:      
                 self.kick_codes(timestep/2.)
                 first=False
             else:
                 self.kick_codes(timestep)
             
             self.drift_codes(self.time+timestep)
             
             self.time += timestep
             
        if not first:
            self.kick_codes(timestep/2.)
            
            
    def synchronize_model(self):
        """ 
        explicitly synchronize all components
        """
        for x in self.codes:
            if hasattr(x,"synchronize_model"):
                if(self.verbose): print x.__class__.__name__,"is synchronizing",
                x.synchronize_model()    
                if(self.verbose): print ".. done"
    
    def stop(self):
        for one_code in self.codes:
            if hasattr(one_code, "stop"):
                one_code.stop()
    
    def get_potential_at_point(self,radius,x,y,z):
        pot=quantities.zero
        for x in self.codes:
            _pot=x.get_potential_at_point(radius,x,y,z)
            pot=pot+_pot
        return pot
        
    def get_gravity_at_point(self,radius,x,y,z):
        ax=quantities.zero
        ay=quantities.zero
        az=quantities.zero
        for x in self.codes:
            _ax,_ay,_az=x.get_gravity_at_point(radius,x,y,z)
            ax=ax+_ax
            ay=ay+_ay
            az=az+_az
        return ax,ay,az

    @property
    def model_time(self):  
         return self.time
      
    @property
    def potential_energy(self):
        result=quantities.zero
        for x in self.codes:
            result+=x.potential_energy
        return result
    
    @property
    def kinetic_energy(self):  
        result=quantities.zero
        for x in self.codes:
            result+=x.kinetic_energy
        return result #- self.kick_energy
        
    @property
    def thermal_energy(self):  
        result=quantities.zero
        for x in self.codes:
            if hasattr(x,'thermal_energy'):
                result+=x.thermal_energy
        return result
          
    @property
    def particles(self):
        array=[]
        for x in self.codes:
            if hasattr(x,"particles"):
                array.append(x.particles)
        if len(array) == 0:
            raise AttributeError
        elif len(array) == 1:
            return array[0]
        return datamodel.ParticlesSuperset(array)                
    
    @property
    def gas_particles(self):
        array=[]
        for x in self.codes:
            if hasattr(x,"gas_particles"):
                array.append(x.gas_particles)
        if len(array) == 0:
            raise AttributeError
        elif len(array) == 1:
            return array[0]
        return datamodel.ParticlesSuperset(array)                
    
    @property
    def dm_particles(self):
        array=[]
        for x in self.codes:
            if hasattr(x,"dm_particles"):
                array.append(x.dm_particles)
            elif hasattr(x,"particles"):
                array.append(x.particles)
        if len(array) == 0:
            raise AttributeError
        elif len(array) == 1:
            return array[0]
        return datamodel.ParticlesSuperset(array)                

# 'private' functions
    def drift_codes_dt(self,dt):
        self._drift_time+=dt
        self.drift_codes(self._drift_time)

    def drift_codes(self,tend):
        threads=[]
                
        for x in self.codes:
            offset=self.time_offsets[x]
            if hasattr(x,"drift"):
                threads.append(threading.Thread(target=x.drift, args=(tend-offset,)) )
            elif hasattr(x,"evolve_model"):
                threads.append(threading.Thread(target=x.evolve_model, args=(tend-offset,)) )
        
        if self.use_threading:
            for x in threads:
                x.start()
            
            for x in threads:
                x.join()
        else:
            for x in threads:
                x.run()

    def kick_codes(self,dt):
       
        de = quantities.zero
        for x in self.codes:
            if hasattr(x,"kick"):
                de += x.kick(dt)
        
        self.kick_energy += de
