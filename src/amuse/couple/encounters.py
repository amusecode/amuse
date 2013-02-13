"""
This module defines the classe to handle handle close
encounters between particles. 

It is used by the multiples module.
"""

from amuse.datamodel import Particles
from amuse.units import constants

class AbstractHandleEncounter(object):
    """Abstract base class for all strategies to
    handle encounters.
    
    We have different important scales in the encounter:
    
    1. Small scale of the interaction. This is the smallest distance
       between any two particles in the encounter. Only binaries
       with an aphelion smaller than this distance (times a factor)
       will be handled as a hard binary.
    2. Large scale of the interaction. This is the total diameter of
       the containing sphere of the interaction. For two body interactions
       this is the same as 1. This scale (times a factor)
       will be used to find neighbour particles to include in the handling
       of the encounter.
    3. Initial sphere radius. This is the radius of the sphere containing
       all the particles in the interaction. The sphere itself is centered
       on the center of mass of the particles. This radius is used to
       scale back (or forward) all particles after the interaction 
       calculation is done. 
       
    After an interaction the following should be True of the root 
    particles:
    
    1. The particles are moving apart
    2. The particles on the outside are just insisde the initial sphere
       radius
    3. The distances between all pairs of particles is larger than the 
       small scale of interaction.    
    """
    
    # look for neighbours of the interaction
    # neighbours_factor * large scale of the interaction
    NEIGHBOURS_FACTOR=1
    
    def __init__(self, 
        particles_in_encounter, 
        particles_in_field = None, 
        existing_multiples = None, 
        existing_binaries = None,
        G = constants.G
    ):
        
        if existing_multiples is None:
            self.existing_multiples = Particles()
        else:
            self.existing_multiples = existing_multiples
            
        if existing_binaries is None:
            self.existing_binaries = Particles()
        else:
            self.existing_binaries = existing_binaries
            
        if particles_in_field is None:
            self.particles_in_field = Particles()
        else:
            self.particles_in_field = particles_in_field
            
        self.new_binaries = Particles()
        self.new_multiples = Particles()
        
        self.dissolved_binaries = Particles()
        self.dissolved_multiples = Particles()
        
        self.particles_close_to_encounter = Particles()
        self.particles_in_encounter = particles_in_encounter
        self.particles_in_field = particles_in_field
        self.particles_close_to_encounter = Particles()
        
        self.all_component_particles_in_encounter = Particles()
        self.resolved_component_particles_in_encounter = Particles()
    
    
    def start(self):
        
        self.determine_scale_of_particles_in_the_encounter()
        
        self.select_neighbours_from_field()
        
        self.fill_set_with_all_component_particles_in_encounter()
        
        self.determine_initial_sphere_of_component_particles_in_encounter()
        
        self.move_all_component_particles_to_initial_sphere_frame_of_reference()
        
        self.evolve_particles_in_encounter_until_end_state()
        
        self.determine_structure_of_the_evolved_state()
        
        self.remove_soft_binaries_from_resolve_component_particles()
        
        self.scale_resolved_component_particles_to_initial_sphere()
        
        self.move_resolved_component_particles_from_initial_sphere_frame_of_reference_to_original()
        
    def fill_set_with_all_component_particles_in_encounter(self):
        for x in self.particles_in_encounter:
            components = self.get_component_particles_of_a_particle(x)
            self.all_component_particles_in_encounter.add_particles(components)
            
        for x in self.particles_close_to_encounter:
            components = self.get_component_particles_of_a_particle(x)
            self.all_component_particles_in_encounter.add_particles(components)

    def get_component_particles_of_a_particle(self, particle):
        if particle in self.existing_multiples:
            multiple = particle.as_particle_in_set(self.existing_multiples)
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            result = Particles()
            for node in tree.iter_descendant_leafs():
                result.add_particle(node.particle)
            result.position += particle.position
            result.velocity += particle.velocity
            return result
        else:
            return particle.as_set()

    def determine_scale_of_particles_in_the_encounter(self):
        # determine large scale from the distance of the farthest particle to the center of mass
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_encounter.position-center_of_mass).lengths()
        max_distance = distances.max()
        self.large_scale_of_particles_in_the_encounter = max_distance * 2
        
        # determine small scale from the smallest distance between all pairs in the encounter
        # for two body interaction this scale is the same as the large scale
        positions = self.particles_in_encounter.position
        transpose_positions = positions.reshape((len(self.particles_in_encounter), 1, 3))
        distances_between_all_particles = ((transpose_positions - positions)**2).sum(axis=2).sqrt()
        distances_between_different_particles  = distances_between_all_particles[distances_between_all_particles > 0*max_distance]
        min_distance = distances_between_different_particles.min()
        self.small_scale_of_particles_in_the_encounter = min_distance    
           
    def determine_initial_sphere_of_component_particles_in_encounter(self):
        self.initial_sphere_position = self.all_component_particles_in_encounter.center_of_mass()
        self.initial_sphere_velocity = self.all_component_particles_in_encounter.center_of_mass_velocity()
        
        distances = (self.all_component_particles_in_encounter.position-self.initial_sphere_position).lengths()
        self.initial_sphere_radius = distances.max()
    
    def move_all_component_particles_to_initial_sphere_frame_of_reference(self):
        self.all_component_particles_in_encounter.position -= self.initial_sphere_position
        self.all_component_particles_in_encounter.velocity -= self.initial_sphere_velocity
        
    def move_resolved_component_particles_from_initial_sphere_frame_of_reference_to_original(self):
        self.resolved_component_particles_in_encounter.position += self.initial_sphere_position
        self.resolved_component_particles_in_encounter.velocity += self.initial_sphere_velocity
        
    def select_neighbours_from_field(self): 
        if len(self.particles_in_field) == 0:
            return
            
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_field.position-center_of_mass).lengths()
        
        indices_for_sort = distances.argsort()
        
        sorted_particles = self.particles_in_field[indices_for_sort]
        sorted_distances = distances[indices_for_sort]
        
        for particle, sorted_distance in zip(sorted_particles, sorted_distances):
            if sorted_distance <= self.large_scale_of_particles_in_the_encounter * self.NEIGHBOURS_FACTOR:
                self.particles_close_to_encounter.add_particle(particle)
                
    def get_potential_energy_of_particles_in_field(self, particles, field):
        """
        Returns the potential energy of a set of particles in
        a field given by another set of particles. Implemented in python,
        subclasses should reimplement this function to use a code. 
        """
        return particles.potential_energy_in_field(field)
        

    def evolve_particles_in_encounter_until_end_state(self):
        """
        Resolves the system of the component particles in the encounter
        (the components of the particles in the encounter and
        the componentns of the neighbouring particles). Fills a new
        set with the resolved particles. Implementation on the abstract
        class is a no-op, need to re-implement this on a subclass
        """
        self.resolved_component_particles_in_encounter.add_particles(self.all_component_particles_in_encounter)
        
    
    def determine_structure_of_the_evolved_state(self):
        """
        Based on the evolved solution determine the hierarchical structure
        of the particles (i.e. binary, triples etc). 
        Implementation on the abstract class is a no-op, need to re-implement this on a subclass
        """
        self.resolved_component_particles_in_encounter.child1 = None
        self.resolved_component_particles_in_encounter.child2 = None

    
    def remove_soft_binaries_from_resolve_component_particles(self):
        """
        Remove binaries with a aphelion (largest separation between the
        parts) larger that the small scale of the encounter from the
        resolved component list.
        """
        pass
        
    def scale_resolved_component_particles_to_initial_sphere(self):
        """
        Scale the system so that all particles are just inside the initial sphere.
        Particles should be moving apart.
        Implementation should be equvalent to moving the system back in time (or forward
        if the system is smaller than the initial scale).
        
        Implementation on the abstract class is a no-op, need to re-implement this on a subclass
        """
        pass
    
    def determine_multiples_in_the_resolved_component_particles(self):
        """
        
        """
