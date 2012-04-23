"""
Sinks

This module contains functions to create new sink particles. These can be used 
to model accretion, for example unto protostars or compact objects.
"""

import numpy
from amuse.units import units
from amuse.units.quantities import zero, AdaptingVectorQuantity
from amuse.datamodel import Particle, Particles, ParticlesSubset
from amuse.support.exceptions import AmuseException

__all__ = ["new_sink_particles"]

class SinkParticles(Particles):
    
    def __init__(self, original_particles, sink_radius=None, mass=None, position=None, velocity=None):
        Particles.__init__(self)
        self.add_particles_to_store(original_particles.get_all_keys_in_store())
        object.__setattr__(self, "original_particles", original_particles)
        object.__setattr__(self, "_additional_attributes", [])
        
        self.sink_radius = sink_radius or original_particles.radius
        
        if not hasattr(original_particles, "mass"):
            self.mass = mass or (0 | units.kg)
        
        if not hasattr(original_particles, "x"):
            for attribute, value in zip(["x","y","z"], position or ([0, 0, 0] | units.m)):
                setattr(self, attribute, value)
        
        if not hasattr(original_particles, "vx"):
            for attribute, value in zip(["vx","vy","vz"], velocity or ([0, 0, 0] | units.m / units.s)):
                setattr(self, attribute, value)
    
    def __setattr__(self, attribute_name, value):
        if hasattr(self.original_particles, attribute_name):
            setattr(self.original_particles, attribute_name, value)
        else:
            Particles.__setattr__(self, attribute_name, value)
            if not attribute_name in self._additional_attributes:
                self._additional_attributes.append(attribute_name)
    
    def get_attribute_names_defined_in_store(self):
        if hasattr(self, "original_particles"):
            attribute_names = self.original_particles.get_attribute_names_defined_in_store()
        else:
            attribute_names = []
        attribute_names.extend(self._additional_attributes)
        return attribute_names
    
    def get_all_values_of_attribute_in_store(self, name_of_the_attribute):
        if name_of_the_attribute in self._additional_attributes:
            return Particles.get_all_values_of_attribute_in_store(self, name_of_the_attribute)
        else:
            return self.original_particles.get_all_values_of_attribute_in_store(name_of_the_attribute)
    
    def get_values_in_store(self, keys, attributes, by_key = True):
        if set(attributes).isdisjoint(self._additional_attributes):
            if not by_key:
                keys = self.key
            return self.original_particles.get_values_in_store(keys, attributes, by_key = True)
        if set(attributes).issubset(self._additional_attributes):
            return Particles.get_values_in_store(self, keys, attributes, by_key = by_key)
        
        additional_attributes = list(set(attributes).intersection(self._additional_attributes))
        original_attributes = list(set(attributes).difference(additional_attributes))
        results_additional = Particles.get_values_in_store(self, keys, additional_attributes, by_key = by_key)
        if not by_key:
            keys = self.key
        results_original = self.original_particles.get_values_in_store(keys, original_attributes, by_key = True)
        
        result = []
        for attribute in attributes:
            try:
                result.append(results_additional[additional_attributes.index(attribute)])
            except ValueError:
                result.append(results_original[original_attributes.index(attribute)])
        return result
    
    def add_sinks(self, original_particles, sink_radius=None, mass=None, position=None):
        new_sinks = self.add_particles(original_particles)
        new_sinks.sink_radius = sink_radius or original_particles.radius
        
        if not hasattr(original_particles, "mass"):
            new_sinks.mass = mass or zero
        
#~        if not hasattr(original_particles, "position"):
#~            new_sinks.position = position or [0, 0, 0] | new_sinks.sink_radius.unit
        if not hasattr(original_particles, "x"):
            if not position:
                position = [0, 0, 0] | self.sink_radius.unit
            self.x = position.x
            self.y = position.y
            self.z = position.z
        
        if not hasattr(original_particles, "vx"):
            if not velocity:
                velocity = [zero, zero, zero] | self.sink_radius.unit/self.sink_radius.unit # Will this work for nbody and SI?
            self.vx = velocity.x
            self.vy = velocity.y
            self.vz = velocity.z
        
        object.__setattr__(self, "original_particles", self.original_particles + original_particles)
    
    def add_sink(self, particle):
        self.add_sinks(particle.as_set())
    
    def accrete(self, particles):
        others = (particles - self.get_intersecting_subset_in(particles))
        too_close = []
        for pos, r_squared in zip(self.position, self.sink_radius**2):
            too_close.append(others.select_array(
                lambda p_pos : (p_pos-pos).lengths_squared() < r_squared, ["position"]))
        try:
            all_too_close = sum(too_close, particles[0:0])
        except AmuseException:
            too_close = self.resolve_duplicates(too_close, particles)
            all_too_close = sum(too_close, particles[0:0])
        if len(all_too_close):
            corrected_masses = AdaptingVectorQuantity()
            corrected_positions = AdaptingVectorQuantity()
            corrected_velocities = AdaptingVectorQuantity()
            for subset, m, pos, vel in zip(too_close, self.mass, self.position, self.velocity):
                if len(subset):
                    total_mass = subset.total_mass() + m
                    corrected_masses.append(total_mass)
                    corrected_positions.append((m*pos + subset.total_mass()*subset.center_of_mass())/total_mass)
                    corrected_velocities.append((m*vel + subset.total_mass()*subset.center_of_mass_velocity())/total_mass)
                else:
                    corrected_masses.append(m)
                    corrected_positions.append(pos)
                    corrected_velocities.append(vel)
            self.mass = corrected_masses
            self.position = corrected_positions
            self.velocity = corrected_velocities
            
            particles.remove_particles(all_too_close)
    
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
            strongest_sinks.append(candidate_sinks[numpy.where(attraction==attraction.amax())[0]])
        
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

