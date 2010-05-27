import os
import pickle
import fnmatch

#from scipy.interpolate import bisplrep, bisplev

import numpy as np

from amuse.support.units import units
from amuse.support.units import constants

from amuse.support.data import core

from amuse.legacy.sse.interface import SSE
from amuse.legacy.evtwin.interface import EVtwin
from amuse.legacy.mesa.interface import MESA


class CachedStellarEvolution:

	"""
	The CachedStellarEvolution wraps around any other stellar evolution class and then caches the evolution. In particular,

		a) if the class has not been called before with this particular combination of {model, initial_mass}, it runs the actual underlying stellar evolution code while storing the stellar state parameters after each call to .evolve_model() into a .pkl file in the cache directory;

		b) if the class finds this particular {model, initial_mass} combination to exist in the cache, it uses the results stored in the .pkl files and doesn't run the actual underlying stellar evolution object, thus making the stellar evolution code "run" significantly faster.

	CachedStellarEvolution does not save the internal state of the underlying stellar evolution code and thus does not protect against "cache underruns". In other words, if a star with mass 1.0MSun has been cached until 100Myr, it will throw an exception upon a call to .evolve_model(end_time = 200 | units.Myr). 
	
	CachedStellarEvolution can use partial cache. In other words, when evolving a set of particles with initial masses m={1.0, 2.0, 3.0}MSun with cached results only available for m={1.0, 3.0}MSun, it will use cached results for these and run the underlying stellar evolution model only for m=2.0MSun.

	Examples of the usage of the class can be found in the scripts:
		test/examples/test_HRdiagram_cluster.py
		test/examples/test_HRdiagram_tracks.py
		
	Both scripts implement the argument '--cache', which can be used to determine the location of the cache directory.

	The caching is then implemented by first creating a regular stellar evolution code and then the feeding this to the constructor for CachedStellarEvolution (along with the location of the cache on disk, stored in options.cacheDir for our example).

	if not (options.cacheDir is None):
		print "Using cache directory: %s" % (options.cacheDir)
		code = CachedStellarEvolution(code, options.cacheDir)

	The 'code' object can now be used like any stellar evolution code, e.g. reading the state with code.particles and evolving with code.evolve_model().
	
	In test_HRdiagram_cluster.py, the initial masses are determined randomly, which means that we need to specify a seed in order for the cache to have a meaningful effect. This can be done through adding the additional argument, e.g. '--seed 6' on the command line.

	"""

	# extra fields that should logically reside in the particles-list but are non-quantity values
	class ParticleCache:
		pass

	def __init__(self, baseStellarEvolution, cacheDirectory):
		self.particles = core.Particles()
		self.particlesCache = []

		self.baseStellarEvolution = baseStellarEvolution
		self.cacheDirectory = cacheDirectory
		
		# create cache directory if it has not been created already
		if not os.path.exists(self.cacheDirectory):
			os.makedirs(self.cacheDirectory)
			print "__init__: created cache directory"
		else:
			print "__init__: cache directory exists"
	
	def initialize_module_with_current_parameters(self):
		self.baseStellarEvolution.initialize_module_with_current_parameters()

	def _particle_update(self, index):
	
		particle = self.particles[index]
		particleCache = self.particlesCache[index]

		# if particle is cached: load next state, update
		if particleCache.isCached:
			#print "_particle_update(): using cache for m=%s" % (particle.mass)
			state = pickle.load(particleCache.cachefh)
			particle.mass 			= state['mass']         | units.MSun
			particle.age		   	= state['age']          | units.Myr
			particle.luminosity 	= state['luminosity']   | units.LSun
			particle.temperature	= state['temperature']  | units.K
			particle.stellar_type	= state['stellar_type'] | units.stellar_type
			particle.radius			= state['radius']       | units.RSun

		# particle is not cached: update state, save state
		else:
			#print "_particle_update(): evolving for m=%s" % (particle.mass)
			particle.mass 			= particleCache.baseParticle.mass
			particle.age		   	= particleCache.baseParticle.age
			particle.luminosity 	= particleCache.baseParticle.luminosity
			particle.temperature	= particleCache.baseParticle.temperature
			particle.stellar_type	= particleCache.baseParticle.stellar_type
			particle.radius			= particleCache.baseParticle.radius

			pickle.dump(dict( \
				mass	   		= particleCache.baseParticle.mass.value_in(units.MSun), \
				age		   		= particleCache.baseParticle.age.value_in(units.Myr), \
				luminosity 		= particleCache.baseParticle.luminosity.value_in(units.LSun), \
				temperature		= particleCache.baseParticle.temperature.value_in(units.K), \
				stellar_type	= particleCache.baseParticle.stellar_type.value_in(units.stellar_type), \
				radius			= particleCache.baseParticle.radius.value_in(units.RSun), \
			), particleCache.cachefh)
	
	def initialize_stars(self):

		# loop over all particles
		self.particleCache = []

		for index in range(len(self.particles)):
		
			self.particlesCache.append(self.ParticleCache())
			particle = self.particles[index]
			particleCache = self.particlesCache[index]
			
			particleCache.cachefn = "%s/model-%s-mass-%s.pkl" % (
				self.cacheDirectory,
				self.baseStellarEvolution.__class__.__name__, 
				("%0.8f" % particle.mass.value_in(units.MSun)).zfill(12)
			)

			# open cache file
			if (os.path.isfile(particleCache.cachefn)):
				print "Using cache: %s" % (particleCache.cachefn, )
				particleCache.isCached = True
				particleCache.cachefh = open(particleCache.cachefn, 'rb')

			else:
				print "Creating cache: %s" % (particleCache.cachefn, )
				particleCache.isCached = False
				particleCache.cachefh = open(particleCache.cachefn, 'wb')
				particleCache.baseParticle = self.baseStellarEvolution.particles.add_particle(particle)			

		# initialize uncached stars
		self.baseStellarEvolution.initialize_stars()

		# initialize states for all particles
		for index in range(len(self.particles)):
			self._particle_update(index)


	def evolve_model(self, end_time = None):

		# evolve non-cached models
		if len(self.baseStellarEvolution.particles) > 0:
			self.baseStellarEvolution.evolve_model(end_time)

		# update state on every star once (enough for default timestep advancement)
		for index in range(len(self.particles)):
			self._particle_update(index)
	
		# further advance states for cached stars and end_time != None
		if not (end_time is None):
			for index in range(len(self.particles)):
				particle = self.particles[index]
				particleCache = self.particlesCache[index]

				while particle.age < end_time and particleCache.isCached:
					self._particle_update(index)

	def stop(self):

		self.baseStellarEvolution.stop()

		for index in range(len(self.particles)):
			particle = self.particles[index]
			particleCache = self.particlesCache[index]
			particleCache.cachefh.close()

	
