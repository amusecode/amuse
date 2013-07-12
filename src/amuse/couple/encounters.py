"""
This module defines the classe to handle handle close
encounters between particles. 

It is used by the multiples module.
"""

from amuse.datamodel import Particle
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.units.quantities import as_vector_quantity

from amuse.support import options
from amuse.support import code

import logging
import numpy

class AbstractHandleEncounter(object):
    """Abstract base class for all strategies to handle encounters.
    
    We have different scales in the encounter:
    
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
       all the particles in the interaction (encounter + neighbours).
       The sphere itself is centered on the center of mass of the particles. 
       This radius is used to scale back (or forward) all particles 
       after the interaction calculation is done. 
       
    After an interaction the following should be True of the particles:
    
    1. The particles are moving apart.
    
    2. The particles on the outside are just insisde the initial sphere
       radius.
       
    3. The distances between all pairs of particles is larger than the 
       small scale of interaction.    
    """
    
    # look for neighbours of the interaction
    # neighbours_factor * large scale of the interaction
    NEIGHBOURS_FACTOR=1
    
    # a hard binary is defined as the small scale of the
    # interaction times this factor
    HARD_BINARY_FACTOR=3
    
    # Initial separation for the scattering experiment, relative
    # to the small scale of the interaction.
    SCATTER_FACTOR=10
    
    
    def __init__(self,
        kepler_code = None,
        G = constants.G
    ):
        self.G = G
        self.kepler_orbits = KeplerOrbits(kepler_code)
        self.reset()
        
    def reset(self):
        self.particles_in_field = Particles()
        
        self.particles_in_encounter = Particles()
        self.particles_close_to_encounter = Particles()
        
        self.all_particles_in_encounter = ParticlesSuperset([self.particles_in_encounter, self.particles_close_to_encounter])
        
        self.existing_multiples = Particles()
        self.existing_binaries = Particles()
        
        self.new_binaries = Particles()
        self.new_multiples = Particles()
        
        self.updated_binaries = Particles()
        
        self.dissolved_binaries = Particles()
        self.dissolved_multiples = Particles()
        
        self.captured_singles = Particles()
        self.released_singles = Particles()
        
        self.all_singles_in_encounter = Particles()
        self.all_singles_close_to_encounter = Particles()
        self.all_singles_in_evolve = ParticlesSuperset([self.all_singles_in_encounter, self.all_singles_close_to_encounter])
        
        self.singles_and_multiples_after_evolve = Particles()
    
        self.kepler_orbits.reset()
            
    def execute(self):
        
        self.determine_scale_of_particles_in_the_encounter()
        
        self.select_neighbours_from_field()
        
        self.determine_initial_sphere_of_particles_in_encounter()
        
        self.scale_up_system_if_two_body_scattering()
        
        self.determine_singles_from_particles_and_neighbours_in_encounter()
        
        self.move_all_singles_to_initial_sphere_frame_of_reference()
        
        self.evolve_singles_in_encounter_until_end_state()
        
        self.determine_structure_of_the_evolved_state()
        
        self.remove_soft_binaries_from_evolved_state()
        
        self.determine_multiples_in_the_evolved_state()
        
        self.determine_captured_singles_from_the_multiples()
        
        self.determine_released_singles_from_the_multiples()
        
        self.determine_particles_after_encounter()
        
        self.scale_evolved_state_to_initial_sphere()
        
        self.move_evolved_state_to_original_frame_of_reference()
        
        self.update_positions_of_subsets()
    
    
    def determine_scale_of_particles_in_the_encounter(self):
        # limit all scaling to the sum of radii of two particles
        # the min distance may not be smaller
        # the max distance may not be smaller
        max_sum_radii = None
        radii = self.particles_in_encounter.radius
        for i, radius in enumerate(radii[:-1]):
            max_sum_radii_i = (radius + radii[i+1:]).max()
            if max_sum_radii is None or max_sum_radii_i > max_sum_radii:
                max_sum_radii = max_sum_radii_i
                
        
        # determine large scale from the distance of the farthest particle to the center of mass
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_encounter.position-center_of_mass).lengths()
        max_distance = distances.max() * 2 # times 2 as we are relative to the center of mass
        max_distance = max(max_sum_radii,  max_distance)
        self.large_scale_of_particles_in_the_encounter = max_distance 
        
        # determine small scale from the smallest distance between all pairs in the encounter
        # for two body interaction this scale is the same as the large scale
        positions = self.particles_in_encounter.position
        transpose_positions = positions.reshape((len(self.particles_in_encounter), 1, 3))
        distances_between_all_particles = ((transpose_positions - positions)**2).sum(axis=2).sqrt()
        distances_between_different_particles  = distances_between_all_particles[distances_between_all_particles > 0*max_distance]
        min_distance = distances_between_different_particles.min()
        min_distance = max(max_sum_radii,  min_distance)
        self.small_scale_of_particles_in_the_encounter = min_distance    
        
    def determine_singles_from_particles_and_neighbours_in_encounter(self):
        for x in self.particles_in_encounter:
            components = self.get_singles_of_a_particle(x)
            self.all_singles_in_encounter.add_particles(components)
        for x in self.particles_close_to_encounter:
            components = self.get_singles_of_a_particle(x)
            self.all_singles_close_to_encounter.add_particles(components)
        
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
            return result
        else:
            return particle.as_set()

    def determine_initial_sphere_of_particles_in_encounter(self):
        self.initial_sphere_position = self.all_particles_in_encounter.center_of_mass()
        self.initial_sphere_velocity = self.all_particles_in_encounter.center_of_mass_velocity()
        
        distances = (self.all_particles_in_encounter.position-self.initial_sphere_position).lengths()
        self.initial_sphere_radius = max(distances.max(), self.small_scale_of_particles_in_the_encounter / 2.0)
    
    def move_all_singles_to_initial_sphere_frame_of_reference(self):
        self.all_singles_in_evolve.position -= self.initial_sphere_position
        self.all_singles_in_evolve.velocity -= self.initial_sphere_velocity
        
    def move_evolved_state_to_original_frame_of_reference(self):
        self.particles_after_encounter.position += self.initial_sphere_position
        self.particles_after_encounter.velocity += self.initial_sphere_velocity
        
    def select_neighbours_from_field(self): 
        if len(self.particles_in_field) == 0:
            return
            
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_field.position-center_of_mass).lengths()
        
        near_distance = self.large_scale_of_particles_in_the_encounter * self.NEIGHBOURS_FACTOR
        near_particles = self.particles_in_field[distances <= near_distance]
        
        self.particles_close_to_encounter.add_particles(near_particles)
    
    def scale_up_system_if_two_body_scattering(self):
        if len(self.particles_close_to_encounter) > 0:
            return # no two body scattering if close perturbers
        if not (len(self.particles_in_encounter) == 2):
            return
        
        delta_position, delta_velocity = self.kepler_orbits.expand_binary(
            self.particles_in_encounter,
            self.SCATTER_FACTOR * self.small_scale_of_particles_in_the_encounter
        )
        self.particles_in_encounter.position += delta_position
        self.particles_in_encounter.velocity += delta_velocity
        
        
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
        self.singles_and_multiples_after_evolve.add_particles(self.all_singles_in_evolve)
        
    
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
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        
        nodes_to_break_up = []
        
        hard_binary_radius = self.small_scale_of_particles_in_the_encounter * self.HARD_BINARY_FACTOR
        
        # a branch in the tree is a node with two children
        # the iter_branches will return only the branches under this node
        roots_to_check = list(tree.iter_branches())
        while len(roots_to_check)>0:
            root_node = roots_to_check.pop()
            children = root_node.get_children_particles()
            semimajor_axis, _ = self.kepler_orbits.get_semimajor_axis_and_eccentricity_for_binary_components(
                children[0],
                children[1]
            )
            print "semimajor_axis", semimajor_axis, _, hard_binary_radius,  semimajor_axis < hard_binary_radius
            if semimajor_axis < hard_binary_radius:
                continue
            nodes_to_break_up.append(root_node.particle)
            
            # if we will break up a level in a triple/multiple, we 
            # also will check the binaries under that level.
            roots_to_check.extend(root_node.iter_branches())
            
        # as this is a binary tree with no pointer up the tree
        # we can break up a binary by removing the parent particle
        for root_node in nodes_to_break_up:
            self.singles_and_multiples_after_evolve.remove_particle(root_node)
                
    def scale_evolved_state_to_initial_sphere(self):
        """
        Scale the system so that all particles are just inside the initial sphere.
        Particles should be moving apart.
        Implementation should be equivalent to moving the system back in time (or forward
        if the system is smaller than the initial scale).
        """
        pass
        
        
    
    def determine_multiples_in_the_evolved_state(self):
        """
        Called after culling and scaling the evolved state. What is left
        are:
            1. multiples (binaries, triples etc) that need to be handled
            as a single particle 
            2. singles
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
            binary_lookup_table[binary.child1.key] = binary
            binary_lookup_table[binary.child2.key] = binary
            
        
        #multiples_lookup_table = {}
        #for multiples in self.existing_multiples:
        #    subtree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        #    tree = components.new_binary_tree_wrapper()
        #    leaves = Particles()
        #    for node in tree.iter_descendant_leafs():
        #        leaves.add_particle(node.particle)
        #    key = tuples(leaves.key.sorted())
        #    multiples_lookup_table[key] = binary
        
        # a branch in the tree is a node with two children
        for root_node in tree.iter_branches():
                
            root_particle = root_node.particle
            
            multiple_components = Particles()
            # descendant_leafs are all children and grandchildren and ... without children
            for child in root_node.iter_descendant_leafs():
                component_particle = multiple_components.add_particle(child.particle)
            
            self.update_binaries(root_node, binary_lookup_table)
            
            # put all components in frame of reference of the root particle
            multiple_components.position -= root_particle.position
            multiple_components.velocity -= root_particle.velocity
            
            # create multiple partice and store it
            multiple_particle = root_particle.copy()
            multiple_particle.child1 = None
            multiple_particle.child2 = None
            multiple_particle.components = multiple_components
            multiple_particle.radius = multiple_components.position.lengths().max() * 2
            self.new_multiples.add_particle(multiple_particle)
            
    def determine_captured_singles_from_the_multiples(self):
        for particle in self.particles_in_encounter:
            if particle in self.existing_multiples:
                continue
                
            for multiple in self.new_multiples:
                if particle in multiple.components:
                    self.captured_singles.add_particle(particle)
            
        for particle in self.particles_close_to_encounter:
            if particle in self.existing_multiples:
                continue
                
            for multiple in self.new_multiples:
                if particle in multiple.components:
                    self.captured_singles.add_particle(particle)
        
    
    def determine_released_singles_from_the_multiples(self):
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        for root_node in tree.iter_leafs():
            
            particle = root_node.particle
            if particle in self.particles_in_encounter or particle in self.particles_close_to_encounter:
                continue
            
            found = False
            for multiple in self.new_multiples:
                if particle in multiple.components:
                    found = True
                    break
            
            if not found:
                self.released_singles.add_particle(particle)
            
        
    def determine_particles_after_encounter(self):
        particles_after_encounter = Particles()
        particles_after_encounter.add_particles(self.particles_in_encounter)
        particles_after_encounter.add_particles(self.particles_close_to_encounter)
        particles_after_encounter.remove_particles(self.dissolved_multiples)
        particles_after_encounter.remove_particles(self.captured_singles)
        particles_after_encounter.add_particles(self.released_singles)
        particles_after_encounter.add_particles(self.new_multiples)
        
        channel = self.singles_and_multiples_after_evolve.new_channel_to(particles_after_encounter)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
        self.particles_after_encounter = particles_after_encounter
    
    def update_positions_of_subsets(self):
        channel = self.particles_after_encounter.new_channel_to(self.new_binaries)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
        channel = self.particles_after_encounter.new_channel_to(self.new_multiples)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
        channel = self.particles_after_encounter.new_channel_to(self.released_singles)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
            
    def update_binaries(self, root_node, binary_lookup_table):
        # a binary tree node is a node with two children
        # the children are leafs (have no children of their own)
        if root_node.is_binary():
            binary_found_in_encounter = root_node.particle
            
            children = list(root_node.iter_children())
            key0 = children[0].particle.key
            key1 = children[1].particle.key
            
            if key0 in binary_lookup_table:
                if key1 in binary_lookup_table:
                    binary0 = binary_lookup_table[key0]
                    binary1 = binary_lookup_table[key1]
                    if binary0 is binary1:
                        binary_known_in_system = binary0
                        self.update_binary(binary_found_in_encounter, binary_known_in_system)
                    
                    else:
                        x = self.new_binaries.add_particle(binary_found_in_encounter)
                        self.dissolved_binaries.add_particle(binary0)
                        self.dissolved_binaries.add_particle(binary1)
                        del binary_lookup_table[binary0.child1.key]
                        del binary_lookup_table[binary0.child2.key]
                        del binary_lookup_table[binary1.child1.key]
                        del binary_lookup_table[binary1.child2.key]
                        binary_lookup_table[key0] = x
                        binary_lookup_table[key1] = x
                else:
                    x = self.new_binaries.add_particle(binary_found_in_encounter)
                    binary0 = binary_lookup_table[key0]
                    self.dissolved_binaries.add_particle(binary0)
                    del binary_lookup_table[binary0.child1.key]
                    del binary_lookup_table[binary0.child2.key]
                    binary_lookup_table[key0] = x
                    binary_lookup_table[key1] = x
            elif key1 in binary_lookup_table:
                x = self.new_binaries.add_particle(binary_found_in_encounter)
                binary1 = binary_lookup_table[key1]
                self.dissolved_binaries.add_particle(binary1)
                del binary_lookup_table[binary1.child1.key]
                del binary_lookup_table[binary1.child2.key]
                binary_lookup_table[key0] = x
                binary_lookup_table[key1] = x
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
        
class HandleEncounter(AbstractHandleEncounter):
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code,
        G = nbody_system.G
    ):
        self.resolve_collision_code = resolve_collision_code
        self.interaction_over_code = interaction_over_code
        AbstractHandleEncounter.__init__(
            self,
            kepler_code,
            G
        )
    
    def reset(self):
        AbstractHandleEncounter.reset(self)
        
        self.resolve_collision_code.reset()
        if not self.interaction_over_code is None:
            self.interaction_over_code.reset()
    
    def evolve_singles_in_encounter_until_end_state(self):
        code = self.resolve_collision_code
        code.reset()
        
        print self.all_singles_in_evolve.key
        code.particles.add_particles(self.all_singles_in_evolve)
        
        interaction_over = code.stopping_conditions.interaction_over_detection
        interaction_over.enable()
        
        end_time = 10000 | nbody_system.time
        if len(self.all_singles_in_evolve) == 2:
            end_time = 100 | nbody_system.time
        code.evolve_model(end_time)
        print "evolve done", code.model_time,  interaction_over.is_set()
            
        if interaction_over.is_set():
            # Create a tree in the module representing the binary structure.
            code.update_particle_tree()

            # Return the tree structure.
            code.update_particle_set()
            self.singles_and_multiples_after_evolve.add_particles(code.particles)
        else:        
            raise Exception(
                "Did not finish the small-N simulation before end time {0}".
                format(end_time)
            )
    
    
    def determine_structure_of_the_evolved_state(self):
        """
        Based on the evolved solution determine the hierarchical structure
        of the particles (i.e. binary, triples etc). 
        Structure is determined during evolve singles....
        """
        pass
        
     
    def scale_evolved_state_to_initial_sphere(self):
        """
        Scale the system so that all particles are just inside the initial sphere.
        Particles should be moving apart.
        Implementation should be equivalent to moving the system back in time (or forward
        if the system is smaller than the initial scale).
        """
        
        self.scale_code = ScaleSystem(self.kepler_orbits, self.G)
        self.scale_code.scale_particles_to_sphere(self.particles_after_encounter, self.initial_sphere_radius)
        
        
class KeplerOrbits(object):
    
    def __init__(self, kepler_code):
        self.kepler_code = kepler_code
    
    def reset(self):
        pass
        
    def get_semimajor_axis_and_eccentricity_for_binary_components(self, particle1, particle2):
        
        particles = Particles()
        particles.add_particle(particle1)
        particles.add_particle(particle2)
        
        self.kepler_code.initialize_from_particles(particles)
        
        return self.kepler_code.get_elements()


    def move_binary(self, scale, true_anomaly, receeding):
        
        if receeding:
            # Note: Always end up on an outgoing orbit.  If
            # periastron > scale, we will be just past periapsis.
            if true_anomaly < 0:
                self.kepler_code.advance_to_periastron()
                self.kepler_code.advance_to_radius(scale)
            else:
                if self.kepler_code.get_separation() < scale:
                    self.kepler_code.advance_to_radius(scale)
                else:
                    self.kepler_code.return_to_radius(scale)
        else:
            if true_anomaly > 0:
                self.kepler_code.return_to_periastron()
                self.kepler_code.return_to_radius(scale)
            else:
                if self.kepler_code.get_separation() < scale:
                    self.kepler_code.return_to_radius(scale)
                else:
                    self.kepler_code.advance_to_radius(scale)
                    

            
    def compress_binary(self, particles, scale, receeding = True):
        """
        Returns the change in positions and velocities for 
        the two-body system consisting of 'particle1' and 'particle2'.
        After applying the change the particles will lie
        inside distance 'scale' of one another.  
        The final orbit will be receding (moving away from each other).
        """    
        separation = (particles[1].position - particles[0].position).length()
        if separation <= scale: 
            # particles are already close together, no scaling done
            # AVE is this correct, will the particle(s) be receding?
            #     or should some movement always happen
            return particles.position * 0, particles.velocity * 0
        
        
        self.kepler_code.initialize_from_particles(particles)

        true_anomaly, mean_anomaly = self.kepler_code.get_angles()
        semimajor_axis, eccentricity = self.kepler_code.get_elements()
        
        periapsis, apoapsis = self.get_periapsis_and_apoapsis(semimajor_axis, eccentricity)

        # closest distance plus 1% of the distance between peri and apo
        limit = periapsis + 0.01*(apoapsis-periapsis)
        if periapsis < scale and limit > scale:
            limit = scale
            
        # we cannot scale to smaller than the periapsis distance
        if scale < limit:
            scale = limit
            
        self.move_binary(scale, true_anomaly, receeding)
        
        rel_position = as_vector_quantity(self.kepler_code.get_separation_vector())
        rel_velocity = as_vector_quantity(self.kepler_code.get_velocity_vector())
    
        return self.deltas_to_update_binary(particles, rel_position, rel_velocity)
    
    def expand_binary(self, particles, scale, receeding = False):
        """
        Returns the change in positions and velocities for 
        the two-body system consisting of 'particle1' and 'particle2'.
        After applying the change the particles will lie
        close to distance 'scale' of one another.  
        The particles will be moving towards each other.
        """
        
        separation = (particles[1].position - particles[0].position).length()     
        
        if separation > scale:
            return particles.position * 0, particles.velocity * 0
            
        
        self.kepler_code.initialize_from_particles(particles)

        true_anomaly,   mean_anomaly = self.kepler_code.get_angles()
        semimajor_axis, eccentricity = self.kepler_code.get_elements()
        
        periapsis, apoapsis = self.get_periapsis_and_apoapsis(semimajor_axis, eccentricity)
        # largest distance minus 1% of the distance between peri and apo
        limit = apoapsis - 0.01*(apoapsis-periapsis)

        # we cannot scale to larger than the apoapsis distance
        # but if eccentricity > 1 we can!
        if eccentricity <= 1 and scale > limit:
            scale = limit
            
        self.move_binary(scale, true_anomaly, receeding)
    
        rel_position = as_vector_quantity(self.kepler_code.get_separation_vector())
        rel_velocity = as_vector_quantity(self.kepler_code.get_velocity_vector())
        
        return self.deltas_to_update_binary(particles, rel_position, rel_velocity)
        
    def deltas_to_update_binary(self, particles, relative_position, relative_velocity):
        total_mass = particles.mass.sum()
        center_of_mass_position = particles.center_of_mass()
        center_of_mass_velocity = particles.center_of_mass_velocity()
        
        positions_to_center_of_mass = particles.position - center_of_mass_position
        velocities_to_center_of_mass = particles.velocity - center_of_mass_velocity
        
        f = particles[1].mass / total_mass 
        fractions = numpy.asarray([f, -(1-f)]).reshape(2,1)
        
        delta_positions  = (relative_position * fractions) - positions_to_center_of_mass
        delta_velocities = (relative_velocity * fractions) - velocities_to_center_of_mass
        
        return delta_positions, delta_velocities
    
    def get_periapsis_and_apoapsis(self, semimajor_axis, eccentricity):
        
        # periapsis == smallest distance
        # apoapsis == largest distance
        if eccentricity < 1: 
            # we hava an ellipsis
            periapsis = semimajor_axis * (1-eccentricity)
            apoapsis  = semimajor_axis * (1+eccentricity)
        else: 
            # we have a parabola or a hyperbola
            periapsis = semimajor_axis * (eccentricity-1)
            apoapsis  = semimajor_axis + periapsis
            # apoapsis is infinity, but this is better
            # as we only use it for the limit
        
        return periapsis, apoapsis
        
    def sample_binary(self, particles, end_time, number_of_points):
        
        self.kepler_code.initialize_from_particles(particles)
        
        sample_points = quantities.linspace(0.0 * end_time, end_time, number_of_points)
        positions = [] | particles.x.unit
        velocities = [] | particles.vx.unit
        for time in sample_points:
            self.kepler_code.transform_to_time(time)
            
            rel_position = as_vector_quantity(self.kepler_code.get_separation_vector())
            rel_velocity = as_vector_quantity(self.kepler_code.get_velocity_vector())
            
            dpos, dvel = self.deltas_to_update_binary(particles, rel_position, rel_velocity)
            positions.append(particles.position + dpos)
            velocities.append(particles.velocity + dvel)
        return positions, velocities
    
class ScaleSystem(object):
    
    def __init__(self, kepler_orbits, G = nbody_system.G):
        self.kepler_orbits = kepler_orbits
        self.G = G
    
    def move_particle(self, particle, delta_position, delta_velocity):
        """"
        Move a particle and all of its descendants by delta position
        and velocity
        """
        
        particle.position += delta_position
        particle.velocity += delta_velocity
        
        tree = particle.as_set().new_binary_tree_wrapper()
        
        descendants = tree.get_descendants_subset()
        descendants.position += delta_position
        descendants.velocity += delta_velocity
    
    def get_particles_with_minimum_separation(self, particles):
        positions = particles.position
        radii = particles.radius
        
        minimum_separation = None
        for i in range(len(particles) - 1):
            i_position = positions[i]
            j_positions = positions[i+1:]
            
            i_radius = radii[i+1:]
            j_radii = radii[i+1:]
            
            delta_positions = i_position - j_positions
            dr = delta_positions.lengths()
            sum_radii = i_radius + j_radii
            
            delta = dr - sum_radii
            index = delta.argmin()
            min_delta =  delta[index]
            if (
                    minimum_separation is None 
                or 
                    min_delta < minimum_separation
            ):
                minimum_separation = min_delta
                particle_i = particles[i]
                particle_j = particles[i+1+index]
        return particle_i, particle_j
            


    
    def scale_particles_to_sphere(self, particles, radius):
        """
        Rescale the system of particles to lie within a sphere
        of the given radius.
        System is moved to the center of mass.
        System may be compressed or expanded.
        
        note:: 
            radius can be zero -> the system will be scaled to minimum seperation (using the radii), 
            the radii of the particles can be zero -> the system will be scaled to radius
            note this is not implemented for 2 body yet!!!
        """
        

        center_of_mass_position = particles.center_of_mass()
        center_of_mass_velocity = particles.center_of_mass_velocity()

        particles.position -= center_of_mass_position
        particles.velocity -= center_of_mass_velocity
        
        # special case, 1 body
        if len(particles) == 1:
            "The position and velocity of this particle must be zero"
            return
        
        kinetic_energy = particles.kinetic_energy()
        potential_energy = particles.potential_energy(G = self.G)
        particle0, particle1 = self.get_particles_with_minimum_separation(particles)
        sphere_radius = particles.position.lengths().max()
        
        distance = (particle0.position - particle1.position).length()
        sum_of_radii = particle0.radius + particle1.radius
        separation = distance - sum_of_radii
        
            
        # special case, 2 bodies, we can use kepler to 
        # do the scaling in a consistent, energy perserving way
        if len(particles) == 2:
            print distance, sum_of_radii, radius
            if distance < sum_of_radii:
                scale = max(2*radius, sum_of_radii)
                delta_p, delta_v = self.kepler_orbits.expand_binary(particles, scale, receeding = True)
            elif separation > radius:
                scale = max(2 * radius, sum_of_radii)
                delta_p, delta_v = self.kepler_orbits.compress_binary(particles, scale, receeding = True)
            else:
                delta_p, delta_v = self.kepler_orbits.expand_binary(particles, 2 * radius, receeding = True)
            particles.position += delta_p
            particles.velocity += delta_v
            particles.velocity += center_of_mass_velocity
            return
        
        
        # for all other situations, we revert to scaling
        # where we perserve energy by scaling
        # the velocities
    
        # we need to scale up, as the separation between particles is less than zero
        if distance < sum_of_radii:
            # use the largest scaling factor
            factor_position = max(sum_of_radii / distance, (2 * radius) / distance)
            
        # we need to scale up, as the minimum distance is less than the sphere diameter
        elif distance < 2 * radius:
            factor_position = (2 * radius) / distance
            
        # we need to scale down, the minimum distance is larger than the radius
        else:
            # we have room to scale down
            if distance > sum_of_radii:
                if (2 * radius) < sum_of_radii:
                    factor_position = sum_of_radii / distance
                else:
                    factor_position = (2 * radius) / distance
            # we have no room available for any scaling
            else:
                factor_position = 1.0
        factor_velocity_squared = 1.0 - (1.0/factor_position-1.0) * potential_energy/kinetic_energy
        if factor_velocity_squared < 0.0:
            from amuse.units import units
            print particles.position
            print particles.velocity
            print particles.radius
            print "radius", radius
            print "distance", distance.as_quantity_in(units.AU)
            print "sum_of_radii", sum_of_radii
            print "factor_position", factor_position
            raise Exception("cannot scale the velocities")
            
        factor_velocity = numpy.sqrt(factor_velocity_squared)
        
        particles.position = center_of_mass_position + factor_position*(particles.position-center_of_mass_position)
        particles.velocity = center_of_mass_velocity + factor_velocity*(particles.velocity-center_of_mass_velocity)


class Binaries(Particles):
    
    def __init__(self, singles):
        Particles.__init__(self)
        self._private.singles = singles
        self.add_particle_function_attribute('components', self.get_children_subset)
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        if len(keys) == 0:
            return
            
        given_attributes = set(attributes)
        
        if not "child1" in given_attributes:
            raise Exception("a binary must always have a child1 attribute")
            
        if not "child2" in given_attributes:
            raise Exception("a binary must always have a child2 attribute")
            
        all_attributes = []
        all_values = []
        for attribute, value in zip(attributes, values):
            all_attributes.append(attribute)
            if attribute == 'child1' or attribute == 'child2':
                value = value.copy_with_link_transfer(None, self._private.singles)
                all_values.append(value)
            else:
                all_values.append(value)
        
        return super(Binaries, self).add_particles_to_store(keys, all_attributes, all_values)
    
    def remove_particles_from_store(self, keys):
        if len(keys) == 0:
            return
        
        pass    
        
        return super(Binaries, self).remove_particles_from_store(keys)
    def get_children_subset(self, binaries, particle):
        return self._private.singles._subset(keys = (particle.child1.key, particle.child2.key))
        
class MultiplesStoppingConditions(object):
    
    def __init__(self):
        self.multiples_change_detection = code.StoppingCondition('multiples_change_detection')
        self.binaries_change_detection = code.StoppingCondition('binaries_change_detection')
        self.encounter_detection = code.StoppingCondition('encounter_detection')
    
    def unset(self):
        self.multiples_change_detection.unset()
        self.binaries_change_detection.unset()
        self.encounter_detection.unset()
        
    def disable(self):
        self.multiples_change_detection.disable()
        self.binaries_change_detection.disable()
        self.encounter_detection.disable()
        
    def is_set(self):
        return (
            self.multiples_change_detection.is_set() or 
            self.binaries_change_detection.is_set() or 
            self.encounter_detection.is_set()
        )
        
class Multiples(options.OptionalAttributes):

    def __init__(self, 
            gravity_code = None,
            handle_encounter_code = None,
            G = nbody_system.G,
            **opts
        ):
            
        options.OptionalAttributes.__init__(self, **opts)
    
       
        self.gravity_code = gravity_code
        self.handle_encounter_code = handle_encounter_code
        self.G = G
        
        self.stopping_conditions = MultiplesStoppingConditions()
        
        self.reset()
        
        
    
    def reset(self):
        self.particles = Particles()
        
        self.multiples = Particles()
        self.singles   = Particles()
        self.singles_in_binaries = Particles()
        self.binaries  = Binaries(self.singles_in_binaries)
        
        self.gravity_code.reset()
        self.stopping_condition = self.gravity_code.stopping_conditions.collision_detection
        self.stopping_condition.enable()
        
        self.channel_from_code_to_model = self.gravity_code.particles.new_channel_to(self.particles)
        self.channel_from_model_to_code = self.particles.new_channel_to(self.gravity_code.particles)
        
        self.stopping_conditions.disable()
    
    def commit_particles(self):
        if len(self.multiples) == 0:
            if not len(self.binaries) == 0:
                for binary in self.binaries:
                    multiple = self.multiples.add_particle(binary)
                    components = binary.components().copy()
                    components.child1 = None
                    components.child2 = None
                    multiple.components = components
                    multiple.mass = components.mass.sum()
                    # TODO radius!
                    multiple.radius = (binary.child1.position - binary.child2.position).length() * 2
                    multiple.position = components.center_of_mass()
                    multiple.velocity = components.center_of_mass_velocity()
                    components.position -= multiple.position 
                    components.velocity -= multiple.velocity 
                    #self.multiples.add_particle(multiple)
        
        if len(self.singles) == 0:
            self.singles.add_particles(self.particles)
            
        if len(self.particles) == 0:
            self.particles.add_particles(self.singles)
            self.particles.add_particles(self.multiples)
            
        self.gravity_code.particles.add_particles(self.particles)
        
        
    def evolve_model(self, time):
        self.stopping_conditions.unset()
        
        self.model_time = self.gravity_code.model_time
        
        previous_time = None
        while self.model_time < time:
            self.gravity_code.evolve_model(time)
            self.model_time = self.gravity_code.model_time
            self.channel_from_code_to_model.copy()
            
            self.particles.new_channel_to(self.multiples).copy_attributes(["x","y","z","vx","vy","vz"])
            if self.stopping_condition.is_set():
                self.handle_stopping_condition()
                self.particles.synchronize_to(self.gravity_code.particles)
                self.channel_from_model_to_code.copy()
            
            if not previous_time is None and previous_time == self.model_time:
                break
            
            if self.stopping_conditions.is_set():
                break
            
            previous_time = self.model_time
            
            if len(self.particles) == 1:
                break
        
    def synchronize_model(self):
        """
        updates the singles to the right position
        """
        
        
    @property
    def all_singles(self):
        result = self.particles.copy()
        for multiple in self.multiples:
            result.remove_particle(multiple)
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            subset = Particles()
            for node in tree.iter_descendant_leafs():
                subset.add_particle(node.particle)
            subset.position += multiple.position
            subset.velocity += multiple.velocity
            result.add_particles(subset)
        return result
        
    def handle_stopping_condition(self):
        encounters = self.determine_encounters()
        if len(encounters[0]) > 1:
            if hasattr(self, 'plot_func'):
                self.plot_func(encounters[0])
        for particles_in_encounter in encounters:
            self.handle_encounter(particles_in_encounter)
    
    def handle_encounter(self, particles_in_encounter):
        code = self.handle_encounter_code
    
        code.reset()
        
        code.particles_in_encounter.add_particles(particles_in_encounter)
        code.particles_in_field.add_particles(self.particles - particles_in_encounter)
        code.existing_binaries.add_particles(self.binaries)
        code.existing_multiples.add_particles(self.multiples)
        
        print "handling encounter"
        print code.particles_in_encounter.key
        #print code.particles_in_encounter.position
        #print code.particles_in_encounter.velocity
        #print code.particles_in_encounter.radius
        #print code.particles_in_encounter.mass
        code.execute()
        
        print "handling encounter done"
        print "number of multiples: ", len(code.new_multiples)
        # update particles (will have singles and multiples)
        self.particles.remove_particles(code.dissolved_multiples)
        self.particles.remove_particles(code.captured_singles)
        self.particles.add_particles(code.new_multiples)
        self.particles.add_particles(code.released_singles)
        
        self.singles.remove_particles(code.captured_singles)
        self.singles.add_particles(code.released_singles)
        
        if self.stopping_conditions.encounter_detection.is_enabled():
            model = Particles()
            print len(code.new_multiples)
            particles_before_encounter = Particles()
            particles_before_encounter.add_particles(code.all_particles_in_encounter)
            particles_after_encounter = Particles()
            particles_after_encounter.add_particles(code.particles_after_encounter)
            
            particle = Particle()
            particle.particles_before_encounter = particles_before_encounter
            particle.particles_after_encounter = particles_after_encounter
            model.add_particle(particle)
            self.stopping_conditions.encounter_detection.set(model)
            
        if self.stopping_conditions.multiples_change_detection.is_enabled():
            if len(code.new_multiples) > 0 or len(code.dissolved_multiples) > 0:
                self.stopping_conditions.multiples_change_detection.set(
                    code.new_multiples,
                    code.dissolved_multiples
                )
                
        # update multiples
        self.multiples.remove_particles(code.dissolved_multiples)
        self.multiples.add_particles(code.new_multiples)
        
        # update binaries
        self.binaries.remove_particles(code.dissolved_binaries)
        self.binaries.add_particles(code.new_binaries)
        
        if self.stopping_conditions.binaries_change_detection.is_enabled():
            if len(code.new_binaries) > 0 or len(code.dissolved_binaries) > 0:
                self.stopping_conditions.binaries_change_detection.set(
                    code.new_binaries,
                    code.dissolved_binaries
                )
        
        channel = code.particles_after_encounter.new_channel_to(self.particles)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
    
    def determine_encounters(self):
        particles0 = self.stopping_condition.particles(0)
        particles1 = self.stopping_condition.particles(1)
        print (particles0.position - particles1.position).lengths()
        encounters = []
        
        from_key_to_encounter = {}
        for particle0, particle1 in zip(particles0, particles1):
            key0 = particle0.key
            key1 = particle1.key
            if key0 in from_key_to_encounter:
                if key1 in from_key_to_encounter:
                    encounter0 = from_key_to_encounter[key0]
                    encounter1 = from_key_to_encounter[key1]
                    if not encounter0 is encounter1:
                        encounter0.add_particles(encounter1)
                        encounter1.remove_particles(encounter1)
                        for x in encounter0:
                            from_key_to_encounter[x.key] = encounter0
                else:
                    encounter = from_key_to_encounter[key0]
                    encounter.add_particle(particle1)
                    from_key_to_encounter[key1] = encounter
            elif key1 in from_key_to_encounter:
                encounter = from_key_to_encounter[key1]
                encounter.add_particle(particle0)
                from_key_to_encounter[key0] = encounter
            else:
                encounter = Particles()
                encounter.add_particle(particle0)
                encounter.add_particle(particle1)
                encounters.append(encounter)
                from_key_to_encounter[key0] = encounter
                from_key_to_encounter[key1] = encounter
        
        return [x for x in encounters if len(x) > 0]
