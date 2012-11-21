import os
import sys
import math
import traceback
import numpy

from amuse.units import units
from amuse.units import constants
from amuse.community.sse.interface import SSE
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA

from amuse.community.cachedse.interface import CachedStellarEvolution, ParticlesTimeseries

from amuse import datamodel

def _fb_search_endpoint_reached(star):
    if (0 <= star.stellar_type.value_in(units.stellar_type) and \
            star.stellar_type.value_in(units.stellar_type) <= 9) or  \
            star.stellar_type.value_in(units.stellar_type) == 16:
        return False
    return True

def is_remnant_stellar_type(stellar_type):
  return stellar_type.value_in(units.stellar_type) in range(10,17)

class FallbackStellarEvolution(object):

    """ 
    FallbackStellarEvolution (formerly EVtwin2SSE) started as a modification 
    of the EVtwin stellar evolution code that hands over execution to 
    SSE when EVtwin crashes. It now is set up to work with general codes.
    
    The handover algorithm performs a RMS search on the relative differences of 
    mass, radius and luminosity between the last point known EVtwin state and 
    the entire stellar history of SSE until the first state of a stellar remnant type.

    :argument enforce_monotonic_mass_evolution: (False) flag to enforce that the mass can only go down
    :argument verbose: (True) be verbose?
    :argument rms_weights: ([1.,1.,1.]) RMS weighting for [mass,radius,luminosity] 
    
    changelog:
    2012-11-21 FIP code can handle N>1 stars, end_time, rms_weights, option to enforce monotonic mass evolution. 
    2012-11-20 FIP changed to FallbackStellarEvolution. 
    """
    
    def __init__(self, main_code_factory = EVtwin, fallback_code_factory = SSE,
                     enforce_monotonic_mass_evolution=False,verbose=True,rms_weights=[1.,1.,1.]):

        self._main_se = main_code_factory()
        self._fallback_se = fallback_code_factory()
        self.particles = datamodel.Particles()
        self.EVtwinAgeAtSwitch = dict()
        self.EVtwinException = dict()
        self.ActiveModel = dict()
        self._FBTimeseries=dict()
        self.model_time=0| units.Myr
        self.enforce_monotonic_mass_evolution=enforce_monotonic_mass_evolution
        self.verbose=verbose
        if len(rms_weights)!=3:
          raise Exception("weights should have len 3")
        self.rms_weights=numpy.array(rms_weights,'d')
        self.rms_weights/=self.rms_weights.sum()
        if self.verbose:
          print "started FallbackStellarEvolution with:"
          print "main SE code:", self._main_se.__class__.__name__
          print "fallback SE code:", self._fallback_se.__class__.__name__
          if self.enforce_monotonic_mass_evolution:
            print "enforcing monotonic mass evolution"
          print "normalized rms weights are %5.3f (mass), %5.3f (radius), %5.3f (luminosity)"% \
            (self.rms_weights[0],self.rms_weights[1],self.rms_weights[2])
        
    def cache_underlying_models(self, cacheDir):
        self._main_se = CachedStellarEvolution(self._main_se, cacheDir)
        self._fallback_se = CachedStellarEvolution(self._fallback_se, cacheDir)

# note: commit and recommit of parameters (and also of particles) should be checked at some point
# it works, but could lead to unexpected results (because sse evolution is delayed)
#    def commit_parameters(self):
#        self._main_se.commit_parameters()
#        self._fallback_se.commit_parameters()

    def commit_particles(self):
        new=self.particles.difference(self._main_se.particles)
        removed=self._main_se.particles.difference(self.particles)

        # remove all particles from underlying models
        if len(removed)>0:
            self._main_se.particles.remove_particles(removed)
        if len(removed)>0:
            self._fallback_se.particles.remove_particles(removed)

        # initialize EVtwin, transfer state
        self._main_se.particles.add_particles(new)

        for part in self.particles:
          self.ActiveModel[part]=self._main_se
          self._transfer_state(part)

        # initialize SSE
        sse_part=self._fallback_se.particles.add_particles(new)
        for part in self.particles:
          self._FBTimeseries[part]=ParticlesTimeseries(part.as_particle_in_set(self._fallback_se.particles))
          self._FBTimeseries[part].add_timepoint()
          
        self.model_time=self.particles.age.min()  

    # copy current state from underlying <active model>.particles to self.particles
    def _transfer_state(self,particle,age_offset=None):
        ActiveModelParticle=particle.as_particle_in_set(self.ActiveModel[particle].particles)
        particle.mass = ActiveModelParticle.mass
        particle.age = (ActiveModelParticle.age+age_offset) if age_offset else ActiveModelParticle.age 
        particle.luminosity = ActiveModelParticle.luminosity
        particle.temperature = ActiveModelParticle.temperature
        particle.stellar_type = ActiveModelParticle.stellar_type
        particle.radius = ActiveModelParticle.radius

    def evolve_model(self,tend=None):

        if tend is not None:
          while self.model_time<tend:
            self.evolve_model()
          return

        for particle in self.particles:
            if particle.age>self.model_time: continue

            if self.ActiveModel[particle] == self._main_se:
                evtwin_part=particle.as_particle_in_set(self._main_se.particles)
                try:
                    prev_age = evtwin_part.age
                    evtwin_part.evolve_one_step()
                    if (prev_age == evtwin_part.age):
                        raise Exception("Evtwin model timestep is zero.")
                    self._transfer_state(particle)
    
                # EVtwin crashed; switch to SSE
                except Exception as ex:
    
                    self.EVtwinAgeAtSwitch[particle] = evtwin_part.age
                    self.EVtwinException[particle] = ex
                    self.ActiveModel[particle] = self._fallback_se
                    if self.verbose:
                      print "FallbackStellarEvolution switching models, %s (age = %s) threw exception: %s" % \
                        (self._main_se.__class__.__name__,self.EVtwinAgeAtSwitch[particle],self.EVtwinException[particle])
    
                    # run SSE for just long enough to get data for the RMS search
                    sse_part=particle.as_particle_in_set(self._fallback_se.particles)
                    while not _fb_search_endpoint_reached(sse_part):
                        sse_part.evolve_one_step()
                        self._FBTimeseries[particle].add_timepoint()
                    if self.verbose:
                      print "FallbackStellarEvolution switch: evolved SSE to: %s " % (sse_part.age,)
    
                    sse_track=self._FBTimeseries[particle].particles[0]
                    self._FB_rms_search(evtwin_part, sse_track)
                        # TODO: Add ModelSwitchFailed exception when RMS statistics is above some threshold?
                    if self.verbose:
                      print ("FallbackStellarEvolution switch parameters: %s %s %s %s" %  
                                (sse_track.SSEIndexAtSwitch, sse_track.SSENextStateIndex, sse_track.SSEAgeAtSwitch, sse_track.RMSErrorAtSwitch))
    
                    self._evolve_model_FB(particle)
            
            # model has been switched to SSE
            else:
                self._evolve_model_FB(particle)
        self.model_time=self.particles.age.min()

#	def _plausible_stellar_type_transition(evt_state, sse_state):
#
#		return ( \
#			# no change / advancement in stellar type
#			(evt_state <= sse_state) or \
#
#			# no differentiation between MS star and convective low-mass star
#			(evt_state <= 1 and sse_state <= 1) 
#		)

    # returns the optimal index for SSE rms
    def _FB_rms_search(self, evtwin_star, sse_track):

        sse_track.SSEIndexAtSwitch = float("nan")
        sse_track.SSEAgeAtSwitch = float("nan")
        sse_track.RMSErrorAtSwitch = float("inf")

        # TODO heuristic for fixing non-physical stellar type transitions
        #evtwin_final_known_state = -1
        #for i in range(len(evt_raw['stellar_types'])):
        #	if evt_raw['stellar_types'][i] != 16:
        #		evtwin_final_known_state = evt_raw['stellar_types'][i]

        for i in range(len(sse_track.age)):

            # TODO
            #if not plausible_stellar_type_transition(evtwin_final_known_state, sse_raw['stellar_types'][i]):
            #	continue

            rel_diff_mass = (sse_track.mass[i] - evtwin_star.mass) / evtwin_star.mass
            if self.enforce_monotonic_mass_evolution and rel_diff_mass>0:
              continue
            rel_diff_radius = (sse_track.radius[i] - evtwin_star.radius) / evtwin_star.radius
            rel_diff_luminosity = (sse_track.luminosity[i] - evtwin_star.luminosity) / evtwin_star.luminosity

            rms = (   self.rms_weights[0]*(rel_diff_mass)**2  \
                    + self.rms_weights[1]*(rel_diff_radius)**2  \
                    + self.rms_weights[2]*(rel_diff_luminosity)**2)

            if (rms < sse_track.RMSErrorAtSwitch):
                sse_track.SSEIndexAtSwitch = i
                sse_track.SSENextStateIndex = i
                sse_track.SSEAgeAtSwitch = sse_track.age[i] #- (10E-3 | units.Myr) # ugly way to "cheat the convergence check"
                sse_track.RMSErrorAtSwitch = rms
                            # TODO calculate fudge factors for m, r, L, T?

    #self._transfer_state_FB()
    def _evolve_model_FB(self,star):
        
        sse_part=star.as_particle_in_set(self._fallback_se.particles)
        sse_track=self._FBTimeseries[star].particles[0]

        # advance SSE if necessary
        while (sse_track.SSENextStateIndex >= len(sse_track.age)):
            sse_part.evolve_one_step()
            self._FBTimeseries[star].add_timepoint()

        # update state state
        star.age = sse_track.age[ sse_track.SSENextStateIndex ] - sse_track.SSEAgeAtSwitch + self.EVtwinAgeAtSwitch[star]
        star.mass = sse_track.mass[ sse_track.SSENextStateIndex ]
        star.radius = sse_track.radius[ sse_track.SSENextStateIndex ]
        star.luminosity = sse_track.luminosity[ sse_track.SSENextStateIndex ]
        star.temperature = sse_track.temperature[ sse_track.SSENextStateIndex ]
        star.stellar_type = sse_track.stellar_type[ sse_track.SSENextStateIndex ]

        # advance index
        sse_track.SSENextStateIndex = sse_track.SSENextStateIndex + 1

    def stop(self):
        self._main_se.stop()
        self._fallback_se.stop()
  
if __name__ == '__main__':

    stellar_evolution = FallbackStellarEvolution()

    stars = datamodel.Particles(4)
    stars.mass = [0.5,1.0,40.0,100.0] | units.MSun

    stars = stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()

    print stellar_evolution.model_time
    for star in stars:
      print star.stellar_type, 
      print stellar_evolution.ActiveModel[star].__class__.__name__
    print
    while stellar_evolution.model_time < 13.2 | units.Gyr:

        stellar_evolution.evolve_model()
        print stellar_evolution.model_time
        for star in stars:
          print star.stellar_type,',', 
          print stellar_evolution.ActiveModel[star].__class__.__name__,'|',
        print

    stellar_evolution.stop()
