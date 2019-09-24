import os
import pickle
import fnmatch

#import numpy as np

from amuse.units import units
from amuse.units import constants

from amuse import datamodel
class CachedStellarEvolution:

    """
    The CachedStellarEvolution wraps around any other stellar evolution class 
    and then caches the evolution. In particular,

            a) if the class has not been called before with this particular 
            combination of {model, initial_mass}, it runs the actual underlying 
            stellar evolution code while storing the stellar state parameters after 
            each call to .evolve_model() into a .pkl file in the cache directory;
            
            b) if the class finds this particular {model, initial_mass} combination 
            to exist in the cache, it uses the results stored in the .pkl files and 
            doesn't run the actual underlying stellar evolution object, thus making 
            the stellar evolution code "run" significantly faster.

    CachedStellarEvolution does not save the internal state of the underlying 
    stellar evolution code and thus does not protect against "cache underruns". 
    In other words, if a star with mass 1.0MSun has been cached until 100Myr, it 
    will throw an exception upon a call to .evolve_model(end_time = 200 | 
    units.Myr). 
    
    CachedStellarEvolution can use partial cache. In other words, when evolving 
    a set of particles with initial masses m={1.0, 2.0, 3.0}MSun with cached 
    results only available for m={1.0, 3.0}MSun, it will use cached results and 
    run the underlying stellar evolution model only for m=2.0MSun.

    Examples of the usage of the class can be found in the scripts:
            test/examples/test_HRdiagram_cluster.py
            test/examples/test_HRdiagram_tracks.py
            
    Both scripts implement the argument '--cache', which can be used to 
    determine the location of the cache directory.

    The caching is then implemented by first creating a regular stellar 
    evolution code and then the feeding this to the constructor for 
    CachedStellarEvolution (along with the location of the cache on disk, stored 
    in options.cacheDir in our example).

    if not (options.cacheDir is None):
            print "Using cache directory: %s" % (options.cacheDir)
            code = CachedStellarEvolution(code, options.cacheDir)

    The 'code' object can now be used like any stellar evolution code, e.g. 
    reading the state with code.particles and evolving with code.evolve_model().
    
    In test_HRdiagram_cluster.py, the initial masses are determined randomly, 
    which means that we need to specify a seed in order for the cache to have a 
    meaningful effect. This can be done through adding the additional argument, 
    e.g. '--seed 6' on the command line.

    """

    # extra fields that should logically reside in the particles-list but are non-quantity values
    class ParticleCache:
        pass

    def __init__(self, baseStellarEvolution, cacheDirectory):
        self.particles = datamodel.Particles()
        self.particlesCache = []

        self.baseStellarEvolution = baseStellarEvolution
        self.cacheDirectory = cacheDirectory
        
        # create cache directory if it has not been created already
        if not os.path.exists(self.cacheDirectory):
            os.makedirs(self.cacheDirectory)
                #print "__init__: created cache directory"
        else:
            #print "__init__: cache directory exists"
            pass
    
    def commit_parameters(self):
        self.baseStellarEvolution.commit_parameters()
        
    def initialize_module_with_current_parameters(self):
        self.commit_parameters()

    def _particle_update(self, index):

        particle = self.particles[index]
        particleCache = self.particlesCache[index]

        # if particle is cached: load next state, update
        if particleCache.isCached:
            #print "_particle_update(): using cache for m=%s" % (particle.mass)
            state = pickle.load(particleCache.cachefh)
            if isinstance(state, dict):
                particle.mass                        = state['mass']         | units.MSun
                particle.age                        = state['age']          | units.Myr
                particle.luminosity        = state['luminosity']   | units.LSun
                particle.temperature        = state['temperature']  | units.K
                particle.stellar_type        = state['stellar_type'] | units.stellar_type
                particle.radius                        = state['radius']       | units.RSun
            
            # cache contains an exception
            else:
                #print "_particle_update(): loaded exception: %s" % (state,)
                raise(state)

        # particle is not cached: update state, save state
        else:
            #print "_particle_update(): evolving for m=%s" % (particle.mass)
            particle.mass                        = particleCache.baseParticle.mass
            particle.age                        = particleCache.baseParticle.age
            particle.luminosity        = particleCache.baseParticle.luminosity
            particle.temperature        = particleCache.baseParticle.temperature
            particle.stellar_type        = particleCache.baseParticle.stellar_type
            particle.radius                        = particleCache.baseParticle.radius

            pickle.dump(dict(  
                    mass                        = particleCache.baseParticle.mass.value_in(units.MSun),  
                    age                                = particleCache.baseParticle.age.value_in(units.Myr),  
                    luminosity                = particleCache.baseParticle.luminosity.value_in(units.LSun),  
                    temperature                = particleCache.baseParticle.temperature.value_in(units.K),  
                    stellar_type        = particleCache.baseParticle.stellar_type.value_in(units.stellar_type),  
                    radius                        = particleCache.baseParticle.radius.value_in(units.RSun),  
            ), particleCache.cachefh)
    
    def commit_particles(self):

        # clean up non-cached particles
        for particleCache in self.particlesCache:
            if not particleCache.isCached:
                self.baseStellarEvolution.particles.remove_particle(particleCache.baseParticle)

        # loop over all particles
        self.particlesCache = []

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
                print("Using cache: %s" % (particleCache.cachefn, ))
                particleCache.isCached = True
                particleCache.cachefh = open(particleCache.cachefn, 'rb')

            else:
                print("Creating cache: %s" % (particleCache.cachefn, ))
                particleCache.isCached = False
                particleCache.cachefh = open(particleCache.cachefn, 'wb')
                particleCache.baseParticle = self.baseStellarEvolution.particles.add_particle(particle)                        

        # initialize uncached stars
        self.baseStellarEvolution.commit_particles()

        # initialize states for all particles
        for index in range(len(self.particles)):
            self._particle_update(index)

    def evolve_model(self, end_time = None):

        # evolve non-cached models
        if len(self.baseStellarEvolution.particles) > 0:

            try:
                self.baseStellarEvolution.evolve_model(end_time)

            # store the exception object for all cached particles
            except Exception as ex:
                #print "CachedStellarEvolution.evolve_model(): caching exception:%s" % (ex,)
                for (particle, particleCache) in zip(self.particles, self.particlesCache):
                    if not particleCache.isCached:
                        pickle.dump(ex, particleCache.cachefh)

                raise(ex)
                
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

        # clean up non-cached particles
        for particleCache in self.particlesCache:
            if not particleCache.isCached:
                self.baseStellarEvolution.particles.remove_particle(particleCache.baseParticle)

    # return a list of available cached masses
    def cached_mass(self):
        model_files = fnmatch.filter(os.listdir(self.cacheDirectory), 'model-%s-*.pkl' % (self.baseStellarEvolution.__class__.__name__))
        model_files.sort()
        cached_masses = []

        for model_file in model_files:
            model_fileh = open("%s/%s" % (self.cacheDirectory, model_file), 'rb')
            state = pickle.load(model_fileh)
            cached_masses.append(state['mass'] | units.MSun)
            model_fileh.close()

        return cached_masses

class ParticlesTimeseries:

    """
    Helper object for collecting/keeping stellar evolution states as timeseries. Example use case
    looks something like the following.
    
            stellar_evolution = ...
            Particles = stellar_evolution.add_particles(...)
            ParticlesSeries = ParticlesTimeseries(Particles)
            
            while(...)
                    stellar_evolution.evolve_model()
                    ParticlesSeries.update()
            
            plot(ParticlesSeries.age_Myr(), ParticlesSeries.mass_MSun())

    """

    class ParticleTimeseries:
        pass

    def __init__(self, pParticlesBase):
        self.particlesBase = pParticlesBase.as_set()
        self.particles = []
        
        for i in range(len(self.particlesBase)):
            self.particles.append(self.ParticleTimeseries())

        for particle in self.particles:
            particle.mass                        = [] 
            particle.age                        = [] 
            particle.luminosity        = [] 
            particle.temperature        = [] 
            particle.stellar_type        = [] 
            particle.radius                        = [] 

    def add_timepoint(self):

        for particle, baseParticle in zip(self.particles, self.particlesBase):

            particle.mass.append(baseParticle.mass) 
            particle.age.append(baseParticle.age)
            particle.luminosity.append(baseParticle.luminosity)
            particle.temperature.append(baseParticle.temperature)
            particle.stellar_type.append(baseParticle.stellar_type)
            particle.radius.append(baseParticle.radius)

    def list_value_in(self, list, unit):
        result = []
        for element in list:
            result.append(element.value_in(unit))
        return result

    def mass_MSun(self, index=0):
        return self.list_value_in(self.particles[index].mass, units.MSun)

    def age_Myr(self, index=0):
        return self.list_value_in(self.particles[index].age, units.Myr)

    def radius_RSun(self, index=0):
        return self.list_value_in(self.particles[index].radius, units.RSun)

    def luminosity_LSun(self, index=0):
        return self.list_value_in(self.particles[index].luminosity, units.LSun)

    def temperature_K(self, index=0):
        return self.list_value_in(self.particles[index].temperature, units.K)

    def stellar_type(self, index=0):
        return self.list_value_in(self.particles[index].stellar_type, units.stellar_type)
