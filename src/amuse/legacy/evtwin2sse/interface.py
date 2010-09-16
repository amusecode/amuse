import os
import sys
import math
import traceback

from amuse.support.units import units
from amuse.support.units import constants

from amuse.support.data import core

from amuse.legacy.sse.interface import SSE
from amuse.legacy.evtwin.interface import EVtwin

from amuse.legacy.cachedse.interface import CachedStellarEvolution, ParticlesTimeseries


class EVtwin2SSE:

    """
    EVtwin2SSE is a modification of the EVtwin stellar evolution code that hands 
    over execution to SSE when EVtwin crashes.
    
    The handover algorithm performs a RMS search on the relative differences of 
    mass, radius and luminosity between the last point known EVtwin state and 
    the entire stellar history of SSE until the first state of a stellar remnant type.

    Note that although the code allows running several stars simultaneously the 
    timepoint where EVtwin crashes is strongly dependant on the inital mass of 
    the stars. It is advisable to run the EVtwin2SSE code for single initial 
    masses at a time.
    
    """

    class ParticleCache:
        pass
    
    def __init__(self, pEVtwin = None, pSSE = None):

        # specifying the stellar evolution objects as parameters in the constructor 
        # allows setting up caching in a convenient way
        if pEVtwin is None:
            self._EVtwin = EVtwin()
        else:
            self._EVtwin = pEVtwin

        if pSSE is None:
            self._SSE = SSE()
        else:
            self._SSE = pSSE

        # initialize member variables
        self.particles = core.Particles()
        self.EVtwinAgeAtSwitch = float("nan") | units.Myr
        self.EVtwinException = None
        self.ActiveModel = self._EVtwin # self.ActiveModel.__class__.__name__ contains name of active model
        
        self._EVtwin_particlesh = None
        self._SSE_particlesh = None

    def cache_underlying_models(self, cacheDir):
        self._EVtwin = CachedStellarEvolution(self._EVtwin, cacheDir)
        self._SSE = CachedStellarEvolution(self._SSE, cacheDir)

    def initialize_module_with_current_parameters(self):
        self._EVtwin.initialize_module_with_current_parameters()
        self._SSE.initialize_module_with_current_parameters()

    def initialize_stars(self):
        self.ActiveModel = self._EVtwin # self.ActiveModel.__class__.__name__ contains name of active model

        # remove all particles from underlying models
        if not (self._EVtwin_particlesh is None):
            self._EVtwin.particles.remove_particles(self._EVtwin_particlesh)
        if not (self._SSE_particlesh is None):
            self._SSE.particles.remove_particles(self._SSE_particlesh)

        # initialize EVtwin, transfer state
        self._EVtwin_particlesh = self._EVtwin.particles.add_particles(self.particles)
        self._EVtwin.initialize_stars()
        self._transfer_state_EVtwin()

        # initialize SSE
        self._SSE_particlesh = self._SSE.particles.add_particles(self.particles)
        self._SSE.initialize_stars()
        self._SSETimeseries = ParticlesTimeseries(self._SSE.particles)
        self._SSETimeseries.add_timepoint()

    # copy current state from underlying <active model>.particles to self.particles
    def _transfer_state_EVtwin(self):

        for particle, ActiveModelParticle in zip(self.particles, self.ActiveModel.particles):        
            particle.mass                        = ActiveModelParticle.mass
            particle.age                        = ActiveModelParticle.age
            particle.luminosity        = ActiveModelParticle.luminosity
            particle.temperature        = ActiveModelParticle.temperature
            particle.stellar_type        = ActiveModelParticle.stellar_type
            particle.radius                        = ActiveModelParticle.radius

    def evolve_model(self): # TODO add evolve_until argument

        # EVtwin is working so far
        #if math.isnan(self.EVtwinAgeAtSwitch.value_in(units.Myr)):
        if self.ActiveModel == self._EVtwin:
            #print "evolve_model() EVtwin"	

            # try advancing a timestep with EVtwin
            try:
                prev_age = self._EVtwin.particles.age
                self._EVtwin.evolve_model()
                if (prev_age == self._EVtwin.particles.age):
                    raise Exception("Evtwin model timestep is zero.")
                self._transfer_state_EVtwin()

            # EVtwin crashed; switch to SSE
            except Exception as ex:

                self.EVtwinAgeAtSwitch = self._EVtwin.particles.age
                self.EVtwinException = ex
                self.ActiveModel = self._SSE
                print "Evtwin2SSE switching models, EVtwin (age = %s) threw exception: %s" % (self._EVtwin.particles.age, self.EVtwinException)

                # run SSE for just long enough to get data for the RMS search
                while not _sse_search_endpoint_reached(self._SSE.particles):
                    self._SSE.evolve_model()
                    self._SSETimeseries.add_timepoint()
                print "Evtwin2SSE switch: evolved SSE to: %s " % (self._SSE.particles.age,)

                for evtwin_star, sse_track in zip(self._EVtwin.particles, self._SSETimeseries.particles):
                    self._SSE_rms_search(evtwin_star, sse_track)
                    # TODO: Add ModelSwitchFailed exception when RMS statistics is above some threshold?
                    print "Evtwin2SSE switch parameters: %s %s %s %s" %  
                            (sse_track.SSEIndexAtSwitch, sse_track.SSENextStateIndex, sse_track.SSEAgeAtSwitch, sse_track.RMSErrorAtSwitch)

                self._evolve_model_SSE()

        # model has been switched to SSE
        else:
            #print "evolve_model() SSE"	
            self._evolve_model_SSE()

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
    def _SSE_rms_search(self, evtwin_star, sse_track):

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
            rel_diff_radius = (sse_track.radius[i] - evtwin_star.radius) / evtwin_star.radius
            rel_diff_luminosity = (sse_track.luminosity[i] - evtwin_star.luminosity) / evtwin_star.luminosity

            rms = (rel_diff_mass.value_in(units.none))**2  
                    + (rel_diff_radius.value_in(units.none))**2  
                    + (rel_diff_luminosity.value_in(units.none))**2 

            if (rms < sse_track.RMSErrorAtSwitch):
                sse_track.SSEIndexAtSwitch = i
                sse_track.SSENextStateIndex = i
                sse_track.SSEAgeAtSwitch = sse_track.age[i] - (10E-3 | units.Myr) # ugly way to "cheat the convergence check"
                sse_track.RMSErrorAtSwitch = rms
                            # TODO calculate fudge factors for m, r, L, T?

    #self._transfer_state_SSE()
    def _evolve_model_SSE(self):
        #print "evolve_model_SSE()"
        for evtwin_star, star, sse_track in zip(self._EVtwin.particles, self.particles, self._SSETimeseries.particles):

            # advance SSE if necessary
            while (sse_track.SSENextStateIndex >= len(sse_track.age)):
                self._SSE.evolve_model()
                self._SSETimeseries.add_timepoint()

            # update state state
            star.age = sse_track.age[ sse_track.SSENextStateIndex ] - sse_track.SSEAgeAtSwitch + self.EVtwinAgeAtSwitch
            star.mass = sse_track.mass[ sse_track.SSENextStateIndex ]
            star.radius = sse_track.radius[ sse_track.SSENextStateIndex ]
            star.luminosity = sse_track.luminosity[ sse_track.SSENextStateIndex ]
            star.temperature = sse_track.temperature[ sse_track.SSENextStateIndex ]
            star.stellar_type = sse_track.stellar_type[ sse_track.SSENextStateIndex ]

            # advance index
            sse_track.SSENextStateIndex = sse_track.SSENextStateIndex + 1

    def stop(self):
        self._EVtwin.stop()
        self._SSE.stop()


def _sse_search_endpoint_reached(stars):
    for star in stars:
        if (0 <= star.stellar_type.value_in(units.stellar_type) and  
                star.stellar_type.value_in(units.stellar_type) <= 9) or  
                star.stellar_type.value_in(units.stellar_type) == 16:
            return False
    return True


if __name__ == '__main__':

    # simple test for the EVtwinContSSE script
    stellar_evolution = EVtwinContSSE()
    stellar_evolution.initialize_module_with_current_parameters() 

    star = core.Particle()
    star.mass = 2.0 | units.MSun

    star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.initialize_stars()

    stopped_evolving = False

    stellar_remnant_counter = 10
    print "%s\t%s\t%s\t%s\t%s\t%s" % (star.age, star.mass, star.radius, star.luminosity, star.stellar_type, stellar_evolution.activeModel)
    while stellar_remnant_counter > 0 and not stopped_evolving:

        if (is_remnant_stellar_type( star.stellar_type )):
            stellar_remnant_counter -= 1
        
        previous_age = star.age
        try:
            stellar_evolution.evolve_model()
            stopped_evolving = (star.age == previous_age) # Check whether the age has stopped increasing

        except Exception as ex:
            traceback.print_exc(file=sys.stdout)
#			print str(ex)
#			stopped_evolving = True

        print "%s\t%s\t%s\t%s\t%s\t%s" % (star.age, star.mass, star.radius, star.luminosity, star.stellar_type, stellar_evolution.activeModel)

    stellar_evolution.particles.remove_particle(star)
    stellar_evolution.stop()
