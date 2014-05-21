"""
Sinks

This module contains functions to create new sink particles. These can be used
to model accretion, for example unto protostars or compact objects.
"""

import numpy
from amuse.units import units, quantities
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
                   velocity=None, angular_momentum=None, looping_over="sinks"):
        ParticlesOverlay.__init__(self, original_particles)

        self._private.looping_over=looping_over

        self.sink_radius = sink_radius or original_particles.radius

        if not hasattr(original_particles, "mass"):
            self.mass = mass or (([0.]*len(self)) | units.kg)

        if not hasattr(original_particles, "x"):
            self.position=position  or (([[0.,0.,0.]]*len(self)) | units.m)

        if not hasattr(original_particles, "vx"):
            self.velocity=velocity  or (([[0.,0.,0.]]*len(self)) | units.m/units.s)

        if not hasattr(original_particles, "lx"):
            self.angular_momentum=angular_momentum  or (([[0.,0.,0.]]*len(self)) | units.g*units.m**2/units.s)

    def accrete(self,orgparticles):
        if self._private.looping_over=="sinks":
          return self.accrete_looping_over_sinks(orgparticles)
        else:
          return self.accrete_looping_over_sources(orgparticles)

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
            new_sinks.mass = mass or (([0.]*len(new_sinks)) | units.kg)

        if not hasattr(original_particles, "x"):
            new_sinks.position=position or (([[0.,0.,0.]]*len(new_sinks)) | units.m)

        if not hasattr(original_particles, "vx"):
            new_sinks.velocity=velocity or (([[0.,0.,0.]]*len(new_sinks)) | units.m/units.s)

        if not hasattr(original_particles, "lx"):
            new_sinks.angular_momentum=angular_momentum or (([[0.,0.,0.]]*len(new_sinks)) | units.g*units.m**2/units.s)

    def add_sink(self, particle):
        self.add_sinks(particle.as_set())

    def select_too_close(self, others):
        too_close = []
        for pos, r_squared in zip(self.position, self.sink_radius**2):
            subset = others[(others.position-pos).lengths_squared() < r_squared]
            too_close.append(subset)
        return too_close

    def accrete_looping_over_sinks(self, orgparticles):
        particles=orgparticles.copy()
        others = (particles - self.get_intersecting_subset_in(particles))
        too_close = self.select_too_close(others)
        try:
            all_too_close = sum(too_close, particles[0:0])
        except AmuseException as ex:
            too_close = self.resolve_duplicates(too_close, particles)
            all_too_close = sum(too_close, particles[0:0])
        if len(all_too_close):
            self.aggregate_mass(too_close)
            orgparticles.remove_particles(all_too_close)
        return all_too_close

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

    def accrete_looping_over_sources(self, orgparticles):
        if len(self) == 0:
            return
        particles=orgparticles.copy()
        others = (particles - self.get_intersecting_subset_in(particles))
        too_close = [particles[0:0] for p in self]
        all_too_close=particles[0:0]
        positions=self.position
        masses=self.mass
        sink_radii2=self.sink_radius**2
        for p in others:
            d2=(positions-p.position).lengths_squared()
            a=numpy.where(d2<sink_radii2)[0]
            if len(a) > 0:
                amin=(d2[a]/masses[a]).argmin()
                too_close[a[amin]]+=p
                all_too_close+=p
        if len(all_too_close):
            self.aggregate_mass(too_close)
            orgparticles.remove_particles(all_too_close)
        return all_too_close

    def aggregate_mass(self,too_close):
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


class AbstractShape(object):
    """
    Abstract superclass of all shapes.
    This class defines common code for all shapes.
    """
    def __or__(self, other):
        return CompoundShape(self, other)

    def within_shape(self, position, sources):
        return sources[self.select(position, sources)]

class CompoundShape(AbstractShape):
    def __init__(self, *sub_shapes):
        self.sub_shapes = sub_shapes

    def __or__(self, other_shape):
        if hasattr(other_shape, "sub_shapes"):
            self.sub_shapes += other_shape.sub_shapes
        else:
            self.sub_shapes += (other_shape,)
        return self

    def select(self, position, sources):
        selection = numpy.zeros(len(sources), dtype=bool)
        for sub_shape in self.sub_shapes:
            selection = numpy.logical_or(selection, sub_shape.select(position, sources))
        return selection

class Sphere(AbstractShape):
    def __init__(self, radius):
        self.radius = radius

    def select(self, position, sources):
        return (sources.position-position).lengths_squared() < self.radius**2

class Spheroid(AbstractShape):
    def __init__(self, dimensions, orientation=None):
        self.dimensions = dimensions
        self.orientation=orientation

    def select(self, position, sources):
        normalized_shape = self.dimensions / self.dimensions.max()
        max_r_squared = self.dimensions.max()**2

        rel_position = (sources.position-position)
        if self.orientation is not None:
            #TODO rotate relative_position using orientation
            raise AmuseException("spheroid orientation not implemented yet")
        return (rel_position/normalized_shape).lengths_squared() < max_r_squared

class Disc(Spheroid):
    def __init__(self, radius, height, **kwargs):
        dimensions = quantities.as_vector_quantity([radius, radius, height])
        Spheroid.__init__(self, dimensions, **kwargs)

class NonSphericalSinkParticles(SinkParticles):
    def __init__(self, original_particles, shapes, *args, **kwargs):
        SinkParticles.__init__(self, original_particles, *args, **kwargs)

        if isinstance(shapes, AbstractShape):
            shapes = [shapes] * len(original_particles)

        self.shapes = shapes

    def select_too_close(self, sources):
        too_close = []
        for pos, shape in zip(self.position, self.shapes):
            subset = shape.within_shape(pos, sources)
            too_close.append(subset)
        return too_close

    def accrete_looping_over_sources(self, orgparticles):
        raise AmuseException("Looping over sources not supported for non spherical sink particles")

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
    :argument shapes: the sink particles can be made non spherical by adding a
        shape object for each particle or a single shape for all particles.
        Note that this is slower then spherical accretion without a shape added.
    """
    shapes = keyword_arguments.pop('shapes', None)
    if shapes is None:
        return SinkParticles(original_particles, *list_arguments, **keyword_arguments)
    else:
        return NonSphericalSinkParticles(original_particles, shapes, *list_arguments, **keyword_arguments)

