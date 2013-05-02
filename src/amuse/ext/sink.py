"""
Sinks

This module contains functions to create new sink particles. These can be used 
to model accretion, for example unto protostars or compact objects.
"""

import numpy
from amuse.units import units
from amuse.units.quantities import zero, AdaptingVectorQuantity
from amuse.datamodel import Particle, ParticlesOverlay, ParticlesSubset
from amuse.datamodel import Particles, ParticlesSuperset

from amuse.support.exceptions import AmuseException

__all__ = ["new_sink_particles"]


def angular_momentum(mass,position,velocity):
    """
    Returns the angular momentum of the particles.
    """
    try:
      return mass.reshape((-1,1)) * position.cross(velocity)
    except:
      return mass * position.cross(velocity)

class SinkParticles(ParticlesOverlay):
    
    def __init__(self, original_particles, sink_radius=None, mass=None, position=None, 
                   velocity=None, angular_momentum=None):
        ParticlesOverlay.__init__(self, original_particles)
                
        self.sink_radius = sink_radius or original_particles.radius
        
        if not hasattr(original_particles, "mass"):
            self.mass = mass or (0 | units.kg)
        
        if not hasattr(original_particles, "x"):
            for attribute, value in zip(["x","y","z"], position or ([0, 0, 0] | units.m)):
                setattr(self, attribute, value)
        
        if not hasattr(original_particles, "vx"):
            for attribute, value in zip(["vx","vy","vz"], velocity or ([0, 0, 0] | units.m / units.s)):
                setattr(self, attribute, value)

        if not hasattr(original_particles, "lx"):
            for attribute, value in zip(["lx","ly","lz"], angular_momentum or ([0, 0, 0] | units.kg * units.m**2 / units.s)):
                setattr(self, attribute, value)
  
    def add_particles_to_store(self, keys, attributes = [], values = []):
        (
            (attributes_inbase, values_inbase), 
            (attributes_inoverlay, values_inoverlay)
        ) = self._split_attributes_and_values(attributes, values)
        
        self._private.overlay_set.add_particles_to_store(keys, attributes_inoverlay, values_inoverlay)
        
        #
        # The sink particles have a little different concept of "overlay particles"
        # apparently the sink particles are positioned on all particles (the complete superset gas + sink)
        # and adding sink_particles will not work, subsets must be summed
        #
        particles = self._private.base_set._original_set()._subset(keys)
        self._private.base_set = self._private.base_set + particles
        
    
    def add_sinks(self, original_particles, sink_radius=None, mass=None, position=None, 
                    velocity=None, angular_momentum=None):
        
        new_sinks = self.add_particles(original_particles)
        new_sinks.sink_radius = sink_radius or original_particles.radius
        if not hasattr(original_particles, "mass"):
            new_sinks.mass = mass or (0 | units.kg)
        
        if not hasattr(original_particles, "x"):
            for attribute, value in zip(["x","y","z"], position or ([0, 0, 0] | units.m)):
                setattr(new_sinks, attribute, value)
        
        if not hasattr(original_particles, "vx"):
            for attribute, value in zip(["vx","vy","vz"], velocity or ([0, 0, 0] | units.m / units.s)):
                setattr(new_sinks, attribute, value)

        if not hasattr(original_particles, "lx"):
            for attribute, value in zip(["lx","ly","lz"], angular_momentum or ([0, 0, 0] | units.kg * units.m**2 / units.s)):
                setattr(self, attribute, value)

    def add_sink(self, particle):
        self.add_sinks(particle.as_set())
    
    def accrete(self, orgparticles):
        particles=orgparticles.copy()
        others = (particles - self.get_intersecting_subset_in(particles))
        too_close = []
        for pos, r_squared in zip(self.position, self.sink_radius**2):
            subset = others[(others.position-pos).lengths_squared() < r_squared]
            too_close.append(subset)
        try:
            all_too_close = sum(too_close, particles[0:0])
        except AmuseException as ex:
            too_close = self.resolve_duplicates(too_close, particles)
            all_too_close = sum(too_close, particles[0:0])
        if len(all_too_close):
            corrected_masses = AdaptingVectorQuantity()
            corrected_positions = AdaptingVectorQuantity()
            corrected_velocities = AdaptingVectorQuantity()
            corrected_angular_momenta = AdaptingVectorQuantity()
            for subset, m, pos, vel, Lin in zip(too_close, self.mass, self.position, self.velocity, self.angular_momentum):
                if len(subset):
                    total_mass = subset.total_mass() + m
                    cmpos=(m*pos + subset.total_mass()*subset.center_of_mass())/total_mass
                    cmvel=(m*vel + subset.total_mass()*subset.center_of_mass_velocity())/total_mass
                    L=Lin+angular_momentum(m,pos-cmpos,vel-cmvel)+angular_momentum(subset.mass,subset.position-cmpos,subset.velocity-cmvel).sum(axis=0)
                    corrected_masses.append(total_mass)                    
                    corrected_positions.append(cmpos)
                    corrected_velocities.append(cmvel)
                    corrected_angular_momenta.append(L)
                else:
                    corrected_masses.append(m)
                    corrected_positions.append(pos)
                    corrected_velocities.append(vel)
                    corrected_angular_momenta.append(Lin)                    
            self.mass = corrected_masses
            self.position = corrected_positions
            self.velocity = corrected_velocities
            self.angular_momentum = corrected_angular_momenta
            
            orgparticles.remove_particles(all_too_close)
    
    def resolve_duplicates(self, too_close, particles):
        # Find the particles that are within more than one sink's radius
        duplicates = particles[0:0]
        keys = set()
        for subset in too_close:
            for particle in subset:
                if (particle.key in keys) and (particle.key not in duplicates.key):
                    duplicates += particle
                else:
                    keys.add(particle.key)
        
        # Determine which sink's attraction is strongest
        strongest_sinks = []
        for duplicate in duplicates:
            candidate_sinks = []
            for index, subset in enumerate(too_close):
                if duplicate in subset:
                    candidate_sinks.append(index)
            attraction = self[candidate_sinks].mass/(self[candidate_sinks].position-duplicate.position).lengths_squared()
            strongest_sinks.append(candidate_sinks[numpy.where(attraction==attraction.amax())[0][0]])
        # Define a new list with particles to be accreted, without the duplicates
        result = []
        for index, subset in enumerate(too_close):
            for duplicate, strongest_sink in zip(duplicates, strongest_sinks):
                if duplicate in subset and not index == strongest_sink:
                    subset -= duplicate
            result.append(subset)
        return result
    

def new_sink_particles(original_particles, *list_arguments, **keyword_arguments):
    """
    Returns new sink particles. These are bound to the 'original_particles' in 
    the sense that they share their attributes. However, the sink particles differ 
    from the original ones in two ways, that is to say, they have:
    (1) an additional attribute 'sink_radius'
    (2) a function 'accrete(particles)' to accrete from particles those that lie 
        within this radius.
    
    :argument original_particles: the particles to be modeled as sinks (required)
    :argument sink_radius: the radii of the sinks (default: original_particles.radius)
    :argument mass: masses of the sinks if not supplied by the original_particles (default: zero)
    :argument position: positions of the sinks if not supplied by the original_particles (default: the origin)
    """
    return SinkParticles(original_particles, *list_arguments, **keyword_arguments)

