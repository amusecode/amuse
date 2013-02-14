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
        
        self.updated_binaries = Particles()
        
        self.dissolved_binaries = Particles()
        self.dissolved_multiples = Particles()
        
        self.particles_close_to_encounter = Particles()
        self.particles_in_encounter = particles_in_encounter
        self.particles_in_field = particles_in_field
        self.particles_close_to_encounter = Particles()
        
        self.all_singles_in_encounter = Particles()
        self.singles_and_multiples_after_evolve = Particles()
    
    
    def start(self):
        
        self.determine_scale_of_particles_in_the_encounter()
        
        self.select_neighbours_from_field()
        
        self.determine_singles_from_particles_in_encounter()
        
        self.determine_initial_sphere_of_singles_in_encounter()
        
        self.move_all_singles_to_initial_sphere_frame_of_reference()
        
        self.evolve_singles_in_encounter_until_end_state()
        
        self.determine_structure_of_the_evolved_state()
        
        self.remove_soft_binaries_from_evolved_state()
        
        self.scale_evolved_state_to_initial_sphere()
        
        self.move_evolved_state_to_original_frame_of_reference()
        
        self.determine_multiples_in_the_evolved_state()
        
    def determine_singles_from_particles_in_encounter(self):
        for x in self.particles_in_encounter:
            components = self.get_singles_of_a_particle(x)
            self.all_singles_in_encounter.add_particles(components)
            
        for x in self.particles_close_to_encounter:
            components = self.get_singles_of_a_particle(x)
            self.all_singles_in_encounter.add_particles(components)

    def get_singles_of_a_particle(self, particle):
        if particle in self.existing_multiples:
            multiple = particle.as_particle_in_set(self.existing_multiples)
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            result = Particles()
            for node in tree.iter_descendant_leafs():
                result.add_particle(node.particle)
            result.position += particle.position
            result.velocity += particle.velocity
            
            self.dissolved_multiples.add_particle(multiple)
            self.existing_multiples.remove_particle(multiple)
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
           
    def determine_initial_sphere_of_singles_in_encounter(self):
        self.initial_sphere_position = self.all_singles_in_encounter.center_of_mass()
        self.initial_sphere_velocity = self.all_singles_in_encounter.center_of_mass_velocity()
        
        distances = (self.all_singles_in_encounter.position-self.initial_sphere_position).lengths()
        self.initial_sphere_radius = distances.max()
    
    def move_all_singles_to_initial_sphere_frame_of_reference(self):
        self.all_singles_in_encounter.position -= self.initial_sphere_position
        self.all_singles_in_encounter.velocity -= self.initial_sphere_velocity
        
    def move_evolved_state_to_original_frame_of_reference(self):
        self.singles_and_multiples_after_evolve.position += self.initial_sphere_position
        self.singles_and_multiples_after_evolve.velocity += self.initial_sphere_velocity
        
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
        

    def evolve_singles_in_encounter_until_end_state(self):
        """
        Resolves the system of the component particles in the encounter
        (the components of the particles in the encounter and
        the componentns of the neighbouring particles). Fills a new
        set with the resolved particles. Implementation on the abstract
        class is a no-op, need to re-implement this on a subclass
        """
        self.singles_and_multiples_after_evolve.add_particles(self.all_singles_in_encounter)
        
    
    def determine_structure_of_the_evolved_state(self):
        """
        Based on the evolved solution determine the hierarchical structure
        of the particles (i.e. binary, triples etc). 
        Implementation on the abstract class is a no-op, need to re-implement this on a subclass
        """
        self.singles_and_multiples_after_evolve.child1 = None
        self.singles_and_multiples_after_evolve.child2 = None

    
    def remove_soft_binaries_from_evolved_state(self):
        """
        Remove binaries with a aphelion (largest separation between the
        parts) larger that the small scale of the encounter from the
        resolved component list.
        """
        pass
        
    def scale_evolved_state_to_initial_sphere(self):
        """
        Scale the system so that all particles are just inside the initial sphere.
        Particles should be moving apart.
        Implementation should be equvalent to moving the system back in time (or forward
        if the system is smaller than the initial scale).
        
        Implementation on the abstract class is a no-op, need to re-implement this on a subclass
        """
        pass
    
    def determine_multiples_in_the_evolved_state(self):
        """
        Called after culling and scaling the evolved state. What is left
        are multiples (binaries, triples etc) that need to be handled
        as a single particle or singles
        """
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        
        # TODO it would be nice to update a multiple
        # instead of dissolving and creating it
        # as we do now, we will assume that a multiple with the same
        # singles in it is the same multiple
        #    for multiple in self.dissolved_multiples:
        #       look_up_table[keys of all singles] = multiple
        
        binary_lookup_table = {}
        for binary in self.existing_binaries:
            key = (binary.child1.key,binary.child2.key)
            binary_lookup_table[key] = binary
        
        # a branch in the tree is node with two children
        for root_node in tree.iter_branches():
            root_particle = root_node.particle
            
            multiple_components = Particles()
            # descendants are all children and grandchildren etc. etc.
            for child in root_node.iter_descendants():
                component_particle = multiple_components.add_particle(child.particle)
            
            self.update_binaries(root_node, binary_lookup_table)
            
            # put all components in frame of reference of the root particle
            multiple_components.position -= root_particle.position
            multiple_components.velocity -= root_particle.velocity
            
            # create multiple partice and store it
            multiple_particle = root_particle.copy()
            multiple_particle.components = multiple_components
            self.new_multiples.add_particle(multiple_particle)
            
        
    def update_binaries(self, root_node, binary_lookup_table):
        # a binary tree node is a node with two children
        # the children are leafs (have no children of their own)
        if root_node.is_binary():
            binary_found_in_encounter = root_node.particle
                
            key = tuple([x.particle.key for x in root_node.iter_children()])
            if key in binary_lookup_table:
                binary_known_in_system = binary_lookup_table[key]
                self.update_binary(binary_found_in_encounter, binary_known_in_system)
            else:
                 self.new_binaries.add_particle(binary_found_in_encounter)
        else:
            for branch in root_node.iter_descendant_branches():
                if branch.is_binary():
                    binary_found_in_encounter = branch.particle
                    key = tuple([x.particle.key for x in branch.iter_children()])
                    if key in binary_lookup_table:
                        binary_known_in_system = binary_lookup_table[key]
                        self.update_binary(binary_found_in_encounter, binary_known_in_system)
                    else:
                        self.new_binaries.add_particle(binary_found_in_encounter)
                
    def update_binary(self, binary_found_in_encounter, binary_known_in_system):
        binary_copy = self.updated_binaries.add_particle(binary_known_in_system)
        binary_copy.child1 = binary_found_in_encounter.child1.copy()
        binary_copy.child2 = binary_found_in_encounter.child2.copy()
        binary_copy.position = binary_found_in_encounter.position
        binary_copy.velocity = binary_found_in_encounter.velocity
