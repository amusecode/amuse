"""
This module defines the classe to handle handle close
encounters between particles. 

It is used by the multiples module.
"""

from amuse.datamodel import Particle
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.datamodel import trees
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.units.quantities import as_vector_quantity
from amuse.units.quantities import zero

from amuse.support import options
from amuse.support import code
from amuse.support import interface
from amuse import io

import logging
import numpy
import logging
import sys

LOG_ENERGY = logging.getLogger('energy')
LOG_ENCOUNTER = logging.getLogger('encounter')

# todo to be more compatible with multiples module:
# - not handling the encounter in case of the neighbours
# - splitting a resulting binary in case of a perturber

class AbstractSelectNeighboursMixin(object):

    def define_neighbours_selection_parameters(self, handler):
        pass
        
    def select_neighbours_from_field(self): 
        raise NotImplementedError
        
class EmptySelectNeighboursMixin(AbstractSelectNeighboursMixin):

    def select_neighbours_from_field(self): 
        pass # do not select any neighbours
        

class SelectNeighboursByDistanceMixin(AbstractSelectNeighboursMixin):

    def __init__(self):
        self.neighbours_factor = 1.0
        
        
    def get_neighbours_factor(self):
        return self.neighbours_factor
        
    def set_neighbours_factor(self, value):
        self.neighbours_factor = value
        
    def define_neighbours_selection_parameters(self, handler):
        handler.add_method_parameter(
            "get_neighbours_factor",
            "set_neighbours_factor",
            "neighbours_factor",
            "look for neighbours of the interaction, neighbours_factor * large scale of the interaction",
            default_value = 1.0
        )
        
    def select_neighbours_from_field(self): 
        if len(self.particles_in_field) == 0:
            return
            
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_field.position-center_of_mass).lengths()
        
        near_distance = self.large_scale_of_particles_in_the_encounter * self.neighbours_factor
        near_particles = self.particles_in_field[distances <= near_distance]
        
        self.particles_close_to_encounter.add_particles(near_particles)
        
        LOG_ENCOUNTER.info("neighbor particles (mutliples or singles): {0}".format(self.particles_close_to_encounter.key))
        
        
class SelectNeighboursByPerturbationMixin(AbstractSelectNeighboursMixin):

    def __init__(self):
        self.neighbor_perturbation_limit = 0.1
        self.wide_perturbation_limit = 0.01 
        self.retain_binary_apocenter = True
        
    def get_neighbor_perturbation_limit(self):
        return self.neighbor_perturbation_limit
        
    def set_neighbor_perturbation_limit(self, value):
        self.neighbor_perturbation_limit = value
        
    def get_wide_perturbation_limit(self):
        return self.wide_perturbation_limit
        
    def set_wide_perturbation_limit(self, value):
        self.wide_perturbation_limit = value
        
    def define_neighbours_selection_parameters(self, handler):
        handler.add_method_parameter(
            "get_neighbor_perturbation_limit",
            "set_neighbor_perturbation_limit",
            "neighbor_perturbation_limit",
            "look for neighbours of the interaction if these neighbours might perturb the collission elements",
            default_value = 0.1
        )
        handler.add_method_parameter(
            "get_wide_perturbation_limit",
            "set_wide_perturbation_limit",
            "wide_perturbation_limit",
            "split resulting binary in case of a possible perturber",
            default_value = 0.01
        )
        
    def select_neighbours_from_field(self): 
        if len(self.particles_in_field) == 0:
            return
            
        center_of_mass = self.particles_in_encounter.center_of_mass()
        distances = (self.particles_in_field.position-center_of_mass).lengths()
        perturbation = self.particles_in_field.mass / distances**3 
        max_perturber_index = perturbation.argmax()
        
        self.perturber_in_field = self.particles_in_field[max_perturber_index]
        self.perturber_distance = distances[max_perturber_index]
        
        factor =  0.5*(self.particles_in_encounter.mass.sum())/self.large_scale_of_particles_in_the_encounter**3 
        minimum_perturbation = self.neighbor_perturbation_limit*factor
        print("ENCOUNTERS:", "minimum_perturbation", minimum_perturbation, "radius", self.large_scale_of_particles_in_the_encounter, "factor", factor)
        near_particles = self.particles_in_field[numpy.logical_or(perturbation > minimum_perturbation , distances < self.large_scale_of_particles_in_the_encounter)]
        print("NP:", near_particles)
        LOG_ENCOUNTER.info("perturbations({0}): {1}".format(minimum_perturbation, perturbation[perturbation > minimum_perturbation]))
        
        self.particles_close_to_encounter.add_particles(near_particles)
        
        LOG_ENCOUNTER.info("neighbor particles (mutliples or singles): {0}".format(self.particles_close_to_encounter.key))
        
        
    def remove_soft_binaries_from_evolved_state(self):
        """
        Remove binaries with a aphelion (largest separation between the
        parts) larger that the small scale of the encounter from the
        resolved component list.
        """
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        
        nodes_to_break_up = []
        
        # a branch in the tree is a node with two children
        # the iter_branches will return only the branches under this node
        roots = list(tree.iter_branches())
        roots_to_check = list(roots) 
        
        # a leaf in the tree is a node with no children (a node can have no children or two children)
        # the iter_leafs will return only the leafs under this node
        singles = Particles()
        for x in tree.iter_leafs():
            singles.add_particle(x.particle)
        
        while len(roots_to_check)>0:
            root_node = roots_to_check.pop()
            children = root_node.get_children_particles()
            semimajor_axis, eccentricity = self.kepler_orbits.get_semimajor_axis_and_eccentricity_for_binary_components(
                children[0],
                children[1]
            )
            periapsis, apoapsis = self.kepler_orbits.get_periapsis_and_apoapsis(
                semimajor_axis,
                eccentricity
            )
            binary_scale = apoapsis if self.retain_binary_apocenter else semimajor_axis * 2
            
            others = Particles()
            for x in roots:
                if not x is root_node:
                    others.add_particle(x.particle)
            others.add_particles(singles)
            others.add_particle(self.perturber_in_field)
            distances = (others.position-root_node.particle.position).lengths()
            perturbation = others.mass / distances**3 
            max_perturber_index = perturbation.argmax()
            
            distance = distances[max_perturber_index]
            max_perturbation = perturbation[max_perturber_index]
            max_perturbation = 2*max_perturbation*binary_scale**3/root_node.particle.mass
            print("max_perturbation:", max_perturbation, self.wide_perturbation_limit)
            if max_perturbation > self.wide_perturbation_limit:
                "break it up!"
                nodes_to_break_up.append(root_node.particle)
            print(roots, roots.index(root_node))
            del roots[roots.index(root_node)]
            print(roots)
            # if we will break up a level in a triple/multiple, we 
            # also will check the binaries under that level.
            roots_to_check.extend(root_node.iter_branches())
            
        # as this is a binary tree with no pointer up the tree
        # we can break up a binary by removing the parent particle
        for root_node in nodes_to_break_up:
            self.singles_and_multiples_after_evolve.remove_particle(root_node)
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
       
    After an interaction the following should be true of the particles:
    
    1. The particles are moving apart.
    
    2. The particles on the outside are just insisde the initial sphere
       radius.
       
    3. The distances between all pairs of particles is larger than the 
       small scale of interaction.    
    """
    
    
    def __init__(self,
        kepler_code = None,
        G = constants.G
    ):
        self.G = G
        self.kepler_orbits = KeplerOrbits(kepler_code)
        
        handler = interface.HandleParameters(self)
        self.define_parameters(handler)
        self.parameters = handler.get_attribute('parameters', None)
        
        self.hard_binary_factor = 3.0
        self.scatter_factor = 10.0
        self.small_scale_factor = 3.0
        
        self.reset()
    
    def before_set_parameter(self):
        pass
        
    def before_get_parameter(self):
        pass
        
    def define_parameters(self, handler):
        self.define_neighbours_selection_parameters(handler)
        
        handler.add_method_parameter(
            "get_hard_binary_factor",
            "set_hard_binary_factor",
            "hard_binary_factor",
            "a hard binary is defined as the small scale of the interaction times this factor",
            default_value = 3.0
        )
        handler.add_method_parameter(
            "get_scatter_factor",
            "set_scatter_factor",
            "scatter_factor",
            "Initial separation for the scattering experiment, relative to the small scale of the interaction.",
            default_value = 3.0
        )
        
    def get_hard_binary_factor(self):
        return self.hard_binary_factor
        
    def set_hard_binary_factor(self, value):
        self.hard_binary_factor = value
        
    def get_scatter_factor(self):
        return self.scatter_factor
        
    def set_scatter_factor(self, value):
        self.scatter_factor = value
        
    def reset(self):
        self.particles_in_field = Particles()
        
        self.particles_in_encounter = Particles()
        self.particles_close_to_encounter = Particles()
        self.multiples_in_encounter = Particles()
        
        self.all_particles_in_encounter = ParticlesSuperset([self.particles_in_encounter, self.particles_close_to_encounter])
        
        self.existing_multiples = Particles()
        self.existing_binaries = Particles()
        
        self.new_binaries = Particles()
        self.new_multiples = Particles()
        
        self.updated_binaries = Particles()
        self.updated_multiples = Particles()
        
        self.dissolved_binaries = Particles()
        self.dissolved_multiples = Particles()
        
        self.captured_singles = Particles()
        self.released_singles = Particles()
        
        self.all_singles_in_encounter = Particles()
        self.all_singles_close_to_encounter = Particles()
        self.all_singles_in_evolve = ParticlesSuperset([self.all_singles_in_encounter, self.all_singles_close_to_encounter])
        
        self.singles_and_multiples_after_evolve = Particles()
    
        self.kepler_orbits.reset()
        self.scatter_energy_error = zero
            
    def execute(self):
        
        self.determine_scale_of_particles_in_the_encounter()
        
        self.select_neighbours_from_field()
        
        self.determine_initial_energies()
        
        self.determine_initial_sphere_of_particles_in_encounter()
        print("scale_up system")
        self.scale_up_system_if_two_body_scattering()
        print("determine initial")
        
        self.determine_initial_multiple_energy()
        print("determin singles and energies")
        
        self.determine_singles_and_energies_from_particles_and_neighbours_in_encounter()
        print("determine initial singles")
        
        self.determine_initial_singles_energies()
        print("move all singles")
        self.move_all_singles_to_initial_sphere_frame_of_reference()
        print("evolve singles") 
        self.evolve_singles_in_encounter_until_end_state()
        
        #self.particles_before_scaling = self.singles_and_multiples_after_evolve .copy()
        print("determine structure") 
        self.determine_structure_of_the_evolved_state()
        
        self.scale_evolved_state_to_initial_sphere()
        
        self.remove_soft_binaries_from_evolved_state()
        
        self.determine_multiples_in_the_evolved_state()
        
        self.determine_captured_singles_from_the_multiples()
        
        self.determine_released_singles_from_the_multiples()
        
        self.determine_particles_after_encounter()
        
        
        self.move_evolved_state_to_original_frame_of_reference()
        
        self.update_positions_of_subsets()
        
        self.determine_final_multiple_energy()
        
        self.determine_final_energies()
        
        self.create_energy_report()
        
        self.error_on_fly_away()
    
    def error_on_fly_away(self):
        center_of_mass = self.particles_after_encounter.center_of_mass()
        distances = (self.particles_after_encounter.position - center_of_mass).lengths()
        particles_to_far = self.particles_after_encounter[(distances > self.large_scale_of_particles_in_the_encounter)]
        if len(particles_to_far) > 0:
            print("distances:", distances)
            print("lsi:", self.large_scale_of_particles_in_the_encounter)
            print("scaled:", distances / self.large_scale_of_particles_in_the_encounter)
            #raise Exception('a particle is too far!');
        
    def determine_initial_energies(self):
        
        self.initial_energy =  self.all_particles_in_encounter.kinetic_energy()
        self.initial_energy += self.all_particles_in_encounter.potential_energy(G=self.G)
        
        # particles_close_to_encounter come from the field, 
        # so these have to be removed for the potential calculation
        self.initial_potential_in_field = self.all_particles_in_encounter.potential_energy_in_field(
            self.particles_in_field - self.particles_close_to_encounter,
            G=self.G
        )
        
        LOG_ENERGY.info("E0={0}, PHI0={1}".format(self.initial_energy, self.initial_potential_in_field))
        
    
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
        #max_distance = max(max_sum_radii,  max_distance)
        self.large_scale_of_particles_in_the_encounter = max_distance 
        
        # determine small scale from the smallest distance between all pairs in the encounter
        # for two body interaction this scale is the same as the large scale
        positions = self.particles_in_encounter.position
        transpose_positions = positions.reshape((len(self.particles_in_encounter), 1, 3))
        distances_between_all_particles = ((transpose_positions - positions)**2).sum(axis=2).sqrt()
        distances_between_different_particles  = distances_between_all_particles[distances_between_all_particles > 0*max_distance]
        min_distance = distances_between_different_particles.min()
       
        min_distance = max(max_sum_radii,  min_distance)
        self.small_scale_of_particles_in_the_encounter = self.small_scale_factor * min_distance    
        
    def determine_initial_multiple_energy(self):
        self.initial_multiple_energy = zero
        print(self.particles_in_encounter)
        for x in self.particles_in_encounter:
            energy = self.get_energy_of_a_multiple(x)
            self.initial_multiple_energy += energy
             
        for x in self.particles_close_to_encounter:
            energy = self.get_energy_of_a_multiple(x)
            self.initial_multiple_energy += energy
    
    def determine_final_multiple_energy(self):
        self.final_multiple_energy = zero
        for x in self.particles_after_encounter:
            energy = self.get_final_energy_of_a_multiple(x)
            self.final_multiple_energy += energy
             
        
    def determine_singles_and_energies_from_particles_and_neighbours_in_encounter(self):
        for x in self.particles_in_encounter:
            components = self.break_up_multiple_and_return_singles_of_a_particle(x)
            self.all_singles_in_encounter.add_particles(components)
             
        for x in self.particles_close_to_encounter:
            components = self.break_up_multiple_and_return_singles_of_a_particle(x)
            self.all_singles_close_to_encounter.add_particles(components)
            
            
    def determine_initial_singles_energies(self):
        
        self.initial_singles_energy =  self.all_singles_in_evolve.kinetic_energy()
        self.initial_singles_energy += self.all_singles_in_evolve.potential_energy(G=self.G)
        
        self.delta_phi_1 = self.initial_singles_energy - self.initial_energy - self.initial_multiple_energy
             
        LOG_ENERGY.info("E1={0}, EMul0={1}, DPHI1={2}".format(
                self.initial_singles_energy, 
                self.initial_multiple_energy,
                self.delta_phi_1
            )
        )
      
    
    def determine_final_energies(self):
        
        self.final_kinetic_energy =  self.particles_after_encounter.kinetic_energy()
        self.final_energy = self.final_kinetic_energy + self.particles_after_encounter.potential_energy(G=self.G)
        
        self.final_potential_in_field = self.particles_after_encounter.potential_energy_in_field(
            self.particles_in_field - self.particles_close_to_encounter,
            G=self.G
        )
        
        self.delta_phi_2 = (
            self.initial_singles_energy + self.scatter_energy_error - 
            self.final_multiple_energy - 
            self.final_energy
        )
             
        LOG_ENERGY.info("E3={0}, PHI2={1}, EMul1={2}, DPHI2={3}".format(
                self.final_energy, 
                self.final_potential_in_field,
                self.final_multiple_energy,
                self.delta_phi_2
            )
        )
        
    def create_energy_report(self):
        self.delta_energy = self.final_energy - self.initial_energy
        self.delta_potential_in_field = self.final_potential_in_field - self.initial_potential_in_field
        self.delta_multiple_energy = self.final_multiple_energy - self.initial_multiple_energy
        self.delta_internal_potential = self.delta_phi_2 - self.delta_phi_1
        
        if self.final_kinetic_energy == zero:
            return
            
        tidal_factor = self.delta_energy / self.final_kinetic_energy
        if  abs(tidal_factor) >  1e-2:
            LOG_ENERGY.warn("Tidal correction is needed, {0}".format(tidal_factor))
        
      
    def break_up_multiple_and_return_singles_of_a_particle(self, particle):
        if particle in self.existing_multiples:
            multiple = particle.as_particle_in_set(self.existing_multiples)
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            result = Particles()
            for node in tree.iter_descendant_leafs():
                result.add_particle(node.particle)
            
            result.position += particle.position
            result.velocity += particle.velocity
            
            self.multiples_in_encounter.add_particle(multiple)
            return result
        else:
            return particle.as_set()

    def get_energy_of_a_multiple(self, particle):
        if particle in self.existing_multiples:
            multiple = particle.as_particle_in_set(self.existing_multiples)
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            singles = Particles()
            for node in tree.iter_descendant_leafs():
                singles.add_particle(node.particle)
            
            # tree is stored in rest state,
            # no energy of central particle
            energy  = singles.kinetic_energy()
            energy += singles.potential_energy(G = self.G)
            
            return  energy
        else:
            return zero
            
    def get_final_energy_of_a_multiple(self, particle):
        if particle in self.new_multiples:
            multiple = particle.as_particle_in_set(self.new_multiples)
        elif particle in self.updated_multiples:
            multiple = particle.as_particle_in_set(self.updated_multiples)
        else:
            return zero
        components = multiple.components
        tree = components.new_binary_tree_wrapper()
        singles = Particles()
        for node in tree.iter_descendant_leafs():
            singles.add_particle(node.particle)
        
        # tree is stored in rest state,
        # no energy of central particle
        energy  = singles.kinetic_energy()
        energy += singles.potential_energy(G = self.G)
        
        return  energy
        

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
        
    
    def scale_up_system_if_two_body_scattering(self):
        if len(self.particles_close_to_encounter) > 0:
            return # no two body scattering if close perturbers
        if not (len(self.particles_in_encounter) == 2):
            return
        
        print("initial_scatter_scale:", self.scatter_factor * self.small_scale_of_particles_in_the_encounter, self.small_scale_of_particles_in_the_encounter, self.scatter_factor)
        delta_position, delta_velocity = self.kepler_orbits.expand_binary(
            self.particles_in_encounter,
            self.scatter_factor * self.small_scale_of_particles_in_the_encounter
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
        
        self.scatter_energy_error = zero
    
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
        
        hard_binary_radius = self.small_scale_of_particles_in_the_encounter * self.hard_binary_factor
        
        # a branch in the tree is a node with two children
        # the iter_branches will return only the branches under this node
        roots_to_check = list(tree.iter_branches())
        while len(roots_to_check)>0:
            root_node = roots_to_check.pop()
            children = root_node.get_children_particles()
            semimajor_axis, eccentricity = self.kepler_orbits.get_semimajor_axis_and_eccentricity_for_binary_components(
                children[0],
                children[1]
            )
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
        
        multiple_lookup_table = {}
        for multiple in self.existing_multiples:
            for particle in self.singles_of_a_multiple(multiple):
                multiple_lookup_table[particle.key] = multiple
                
        binary_lookup_table = {}
        for binary in self.existing_binaries:
            binary_lookup_table[binary.child1.key] = binary
            binary_lookup_table[binary.child2.key] = binary
            
        
        # a branch in the tree is a child node with two children
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
            
            existing_multiple = self.lookup_existing_multiple_with_components(multiple_components, multiple_lookup_table)
            
            # create or copy multiple particle and store it
            if existing_multiple is None:
                multiple_particle = root_particle.copy()
                multiples_set = self.new_multiples
            else:
                multiple_particle = existing_multiple.copy()
                multiple_particle.position = root_particle.position
                multiple_particle.velocity = root_particle.velocity
                multiple_particle.mass = root_particle.mass
                multiples_set = self.updated_multiples
                
            multiple_particle.child1 = None
            multiple_particle.child2 = None
            multiple_particle.components = multiple_components
            multiple_particle.radius = self.determine_radius_from_components(multiple_components)
            
            multiples_set.add_particle(multiple_particle)
        
        for multiple in self.multiples_in_encounter:
            if not multiple in self.updated_multiples:
                self.dissolved_multiples.add_particle(multiple)
                
        # a leaft in the tree is a child node with no children
        for root_node in tree.iter_leafs():
            self.update_binaries_from_single(root_node.particle, binary_lookup_table)
    
    def determine_radius_from_components(self, components):
        return components.position.lengths().max() * 2
        
    def singles_of_a_multiple(self, multiple):
        components = multiple.components
        tree = components.new_binary_tree_wrapper()
        singles = Particles()
        for node in tree.iter_descendant_leafs():
            singles.add_particle(node.particle)
        return singles
        
    def lookup_existing_multiple_with_components(self, components, multiple_lookup_table):
        found_multiple = None
        if components[0].key in multiple_lookup_table:
            found_multiple = multiple_lookup_table[components[0].key]
        else:
            return None
        for x in components[1:]:
            if x.key in multiple_lookup_table:
                if not found_multiple == multiple_lookup_table[x.key]:
                    return None
            else:
                return None
        return found_multiple
        
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
        
        channel = self.particles_after_encounter.new_channel_to(self.updated_multiples)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
        channel = self.particles_after_encounter.new_channel_to(self.released_singles)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
        
    def update_binaries_from_single(self, single, binary_lookup_table):
        
        key = single.key   
        if key in binary_lookup_table:
            binary = binary_lookup_table[key]
            self.dissolved_binaries.add_particle(binary)
            del binary_lookup_table[binary.child1.key]
            del binary_lookup_table[binary.child2.key]
            
    def update_binaries(self, root_node, binary_lookup_table):
        # a binary tree node is a node with two children
        # the children are leafs (have no children of their own)
        if root_node.is_binary():
            self.lookup_and_update_binary(root_node, binary_lookup_table)
        else:
            for single in root_node.iter_leafs():
                self.update_binaries_from_single(single.particle, binary_lookup_table)
            for branch in root_node.iter_descendant_branches():
                if branch.is_binary():
                    self.lookup_and_update_binary(branch, binary_lookup_table)
                else:
                    for single in branch.iter_leafs():
                        self.update_binaries_from_single(single.particle, binary_lookup_table)
                        
    def lookup_and_update_binary(self, root_node, binary_lookup_table):
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
    
    def update_binary(self, binary_found_in_encounter, binary_known_in_system):
        binary_copy = self.updated_binaries.add_particle(binary_known_in_system)
        binary_copy.child1 = binary_found_in_encounter.child1.copy()
        binary_copy.child2 = binary_found_in_encounter.child2.copy()
        binary_copy.position = binary_found_in_encounter.position
        binary_copy.velocity = binary_found_in_encounter.velocity

        
        

class HandleEncounterWithCollisionCode(AbstractHandleEncounter):
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code = None,
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
        for x in self.all_singles_in_evolve:
            print(x)
        code = self.resolve_collision_code
        code.reset()
        code.particles.add_particles(self.all_singles_in_evolve)
        
        initial_scatter_energy = code.get_total_energy()
        
        end_time = 10000 | nbody_system.time
        if len(self.all_singles_in_evolve) == 2:
            end_time = 100 | nbody_system.time

        code.evolve_model(0.0001 * end_time)
        interaction_over = code.stopping_conditions.interaction_over_detection
        interaction_over.enable()
        
        LOG_ENCOUNTER.info("evolving singles in encounter")
        print(self.all_singles_in_evolve)
        code.evolve_model(end_time)
        LOG_ENCOUNTER.info("evolving singles in encounter finished model_time = {0}".format(code.model_time))
        
        print("i over:", interaction_over.is_set())
        if interaction_over.is_set():
            # Create a tree in the module representing the binary structure.
            code.update_particle_tree()

            # Return the tree structure.
            code.update_particle_set()

            print(code.particles)
        
            final_scatter_energy = code.get_total_energy()
        
            self.scatter_energy_error = final_scatter_energy - initial_scatter_energy 
        
            self.singles_and_multiples_after_evolve.add_particles(code.particles)
            self.particles_before_scaling = code.particles.copy()
            LOG_ENERGY.info('scatter_energy_error={0}'.format(self.scatter_energy_error))
            return
            
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
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        roots_and_singles = tree.get_children_subset()
        self.scale_code.scale_particles_to_sphere(roots_and_singles, 1.01 * self.initial_sphere_radius)
        

class HandleEncounterWithSmallN(AbstractHandleEncounter):
    debug_encounters = False
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code = None,
        G = nbody_system.G
    ):
        self.resolve_collision_code = resolve_collision_code
        self.interaction_over_code = interaction_over_code
        AbstractHandleEncounter.__init__(
            self,
            kepler_code,
            G
        )
    def evolve_singles_in_encounter_until_end_state(self):

        pre = 'encounter:'	# identifier for all output here

                          
        # Take the system described by particles and evolve it forward
        # in time until it is over.  Don't update global quantities,
        # don't interpret the outcome.  Return the energy error due to
        # the smallN integration.

        # Temporarily avoid "is_over" problems.  If we allow
        # collisions to stop early -- when they become too large or
        # last too long -- then we need will logic to manage the
        # intermediate state that results.  TODO
        final_scatter_scale = 1.e30 | nbody_system.length
        
        
        timescale = self.get_timescale()
        print("time_scale =", timescale)
        if self.resolve_collision_code.unit_converter is None:
            end_time = 1.e4 * abs(timescale) # nbody_system.time
            delta_t = min(10*abs(timescale), 1.0 | nbody_system.time)
        else:
            end_time = 1.e4 * abs(timescale)
            delta_t = 10*abs(timescale)
        print("end_time =", end_time)
        print("delta_t =", delta_t)

        resolve_collision_code = self.resolve_collision_code
        resolve_collision_code.reset()

        time = 0 * end_time
        resolve_collision_code.set_time(time)
        resolve_collision_code.particles.add_particles(self.all_singles_in_evolve)
        resolve_collision_code.commit_particles()
    #self.particles_before_scaling = self.all_singles_in_evolve.copy()

        # Channel to copy values from the code to the set in memory.
        # channel = resolve_collision_code.particles.new_channel_to(particles)

        initial_scatter_energy = self.get_total_energy(resolve_collision_code)

        print(pre, 'number_of_stars =', len(self.all_singles_in_evolve), ' ', self.all_singles_in_evolve.key)
        print(pre, 'initial energy =', initial_scatter_energy)
        #print particles

        delta_t_max = 64*delta_t
        while delta_t_max < end_time/10: 
            delta_t_max *= 2

        if self.debug_encounters:
            delta_t *= 0.1

        initial_delta_t = delta_t
        print(pre, 'evolving to time', end_time)
        print(pre, 'initial step =', initial_delta_t)
        # if self.debug_encounters:
        #     print(pre, '### START ENCOUNTER ###')
        #     print(pre, '### snapshot at time %f' % 0.0)
        #     for p in particles:
        #         print(pre, '### id=%d, x=%f, y=%f, z=%f,'\
        #               'vx=%f, vy=%f, vz=%f' % \
        #                 (p.id, p.x.number, p.y.number, p.z.number,
        #                  p.vx.number, p.vy.number, p.vz.number))

        resolve_collision_code.set_break_scale(final_scatter_scale)

        while time < end_time:

            tt = time
            time += delta_t
            print(pre, '...to time', time)

            # Work with internal steps of initial_delta_t to allow
            # checks for quasi-stable motion.

            while tt < time:

                tt += initial_delta_t
                if tt > time: 
                    tt = time
                print(pre, '    ...', time, tt, \
                      'model_time =', resolve_collision_code.model_time)
                resolve_collision_code.evolve_model(tt)
                print(pre, '    ...back:', \
                      ': model_time =', resolve_collision_code.model_time)
                tt = resolve_collision_code.model_time

                # Note: Return with tt != time means we have exceeded
                # the size limit and don't need to check is_over().

                # DEBUGGING:
                # if self.debug_encounters:
                #     print(pre, '### snapshot at time %f' % time.number)
                #     #resolve_collision_code.update_particle_tree()
                #     #resolve_collision_code.update_particle_set()
                #     resolve_collision_code.particles.synchronize_to(particles)
                #     channel.copy()
                #     for p in particles:
                #             print(pre, '### id=%d, x=%f, y=%f, z=%f,'\
                #               'vx=%f, vy=%f, vz=%f' % \
                #                 (p.id, p.x.number, p.y.number, p.z.number,
                #                  p.vx.number, p.vy.number, p.vz.number))

                # The argument final_scatter_scale is used to limit
                # the size of the system.  It has to be supplied again
                # because the code that determines if the scattering
                # is over isn't necessarily the same as
                # resolve_collision_code.  Currently, only smallN has
                # an "is_over()" function.  TODO
                #
                # Return values:	0 - not over
                #			1 - over
                #			2 - quasi-stable system
                #			3 - not over, but size exceeded limit
                #
                # Note that this is really a stopping condition, and
                # should eventually be handled that way.  TODO

                # We are currently ignoring any possibility of a
                # physical collision during the multiples encounter.
                # TODO

                over = resolve_collision_code.is_over(final_scatter_scale,
                                                      0)    # verbose = 0

                if over:
                    final_scatter_energy = self.get_total_energy(resolve_collision_code)
                    scatter_energy_error = final_scatter_energy - initial_scatter_energy
                    print(pre, 'over =', over, 'at time', tt)
                    #print pre, 'initial energy =', initial_scatter_energy
                    #print pre, 'final energy =', final_scatter_energy
                    #print pre, 'energy error =', scatter_energy_error
                    print(pre, 'fractional energy error =', scatter_energy_error/initial_scatter_energy)
                    if self.debug_encounters:
                            print(pre, '### END ENCOUNTER ###')

                    # Create a tree in the module representing the binary structure.
                    resolve_collision_code.update_particle_tree()

                    # TODO: what happens if we reach over = 2 or 3?

                    # Note that center of mass particles are now part
                    # of the particle set...

                    # Return the tree structure to AMUSE.  Children
                    # are identified by get_children_of_particle in
                    # interface.??, and the information is returned in
                    # the copy operation.

                    resolve_collision_code.update_particle_set()
                    self.singles_and_multiples_after_evolve.add_particles(resolve_collision_code.particles)
                    self.particles_before_scaling = resolve_collision_code.particles.copy()

                    return scatter_energy_error

                if tt >= 0.9999999*time:
                    break

            time = resolve_collision_code.model_time
            if not self.debug_encounters:
                if delta_t < delta_t_max and time > 0.999999*4*delta_t:
                    delta_t *= 2
                    print(pre, 'setting delta_t =', delta_t)

        raise Exception(
            pre + "Did not finish the small-N simulation before end time {0}".
            format(end_time)
        )
    def get_total_energy(self, code):
        # ??? from Steve: what is get_binary_energy()?
        try:
            binaries_energy = code.get_binary_energy()  # include binaries
        except:                                         # if code understands
            binaries_energy = zero
        total_energy = code.potential_energy + code.kinetic_energy \
                             + binaries_energy

        return total_energy

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
        tree = self.singles_and_multiples_after_evolve.new_binary_tree_wrapper()
        roots_and_singles = tree.get_children_subset()
        self.scale_code.scale_particles_to_sphere(roots_and_singles, 1.01 * self.initial_sphere_radius)
        #self.particles_before_scaling = roots_and_singles.copy()
        
    def get_timescale(self):
        self.kepler_orbits 
        self.all_singles_in_evolve
        min_period = None
        for i, iparticle in enumerate(self.all_singles_in_evolve[:-1]):
            for j, jparticle in enumerate(self.all_singles_in_evolve[i+1:]):
                period = self.kepler_orbits.get_period(iparticle, jparticle)
                print("period =", period)
                period = self.kepler_orbits.get_periastron_time(iparticle, jparticle)
                print("time =", period)
                if min_period is None:
                    min_period = period
                else:
                    min_period = min_period.min(period)
        print("tperi =", min_period)
        return min_period
class HandleEncounter(HandleEncounterWithCollisionCode, SelectNeighboursByDistanceMixin):
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code = None,
        G = nbody_system.G
    ):
        HandleEncounterWithCollisionCode.__init__(
            self,
            kepler_code,
            resolve_collision_code,
            interaction_over_code,
            G
        )
        SelectNeighboursByDistanceMixin.__init__(self)
    
        
class StickyHandleEncounter(AbstractHandleEncounter, EmptySelectNeighboursMixin):
    
    def __init__(self, G = nbody_system.G):
        AbstractHandleEncounter.__init__(
            self,
            None,
            G
        )
        
    
    def evolve_singles_in_encounter_until_end_state(self):
        self.scatter_energy_error = quantities.zero
            
        particles = self.all_singles_in_evolve.copy()
        
        working_set = particles.copy()
        parents = Particles()
        
        self.particles_before_scaling = particles.copy()
        counter = len(working_set)
        while counter > 1 :
            
            number_of_particles = len(working_set)
            indices1, indices2 = numpy.triu_indices(number_of_particles, 1)
            dd = lambda x : x[indices1] - x[indices2]
            dx = dd(working_set.x)
            dy = dd(working_set.y)
            dz = dd(working_set.z)
            distances_squared = (dx**2 + dy**2 + dz**2)
            minindex =distances_squared.argmin()
            index1 = indices1[minindex]
            index2 = indices2[minindex]
            partner1 = working_set[index1]
            partner2 = working_set[index2]
            
            mass1 = partner1.mass
            mass2 = partner2.mass
            total_mass = mass1 + mass2
            parent = Particle()
            parent.mass =mass1 + mass1
            parent.position = (partner1.position + partner2.position) / 2.0
            parent.velocity = (
                    (mass1 / total_mass) * partner1.velocity +
                    (mass2 / total_mass) * partner2.velocity
            )
            parent.child1 = partner1
            parent.child2 = partner2
            parents.add_particle(parent)
            
            working_set.remove_particle(partner1)
            working_set.remove_particle(partner2)
            working_set.add_particle(parent)
            
            counter -= 1
        
        result = Particles()
        result.add_particles(particles)
        parents = result.add_particles(parents)
        for x in parents:
            x.child1 = x.child1.as_particle_in_set(result)
            x.child2 = x.child2.as_particle_in_set(result)
        self.singles_and_multiples_after_evolve.add_particles(result)
    
    
    def determine_structure_of_the_evolved_state(self):
        pass
    
    def scale_evolved_state_to_initial_sphere(self):
        pass
    
    def remove_soft_binaries_from_evolved_state(self):
        pass
    
    def scale_up_system_if_two_body_scattering(self):
        pass
    
    def determine_radius_from_components(self, components):
        return components.position.lengths().max() + components.radius.max()

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

    def get_period(self, particle1, particle2):
        
        particles = Particles()
        particles.add_particle(particle1)
        particles.add_particle(particle2)
        
        self.kepler_code.initialize_from_particles(particles)
        
        return self.kepler_code.get_period()
    def get_periastron_time(self, particle1, particle2):
        
        particles = Particles()
        particles.add_particle(particle1)
        particles.add_particle(particle2)
        
        self.kepler_code.initialize_from_particles(particles)
        M,_ = self.kepler_code.get_angles()
        if M < 0:
            self.kepler_code.advance_to_periastron()
        else:
            self.kepler_code.return_to_periastron()
        return self.kepler_code.get_time()


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
            print("true_anomaly:", true_anomaly,  self.kepler_code.get_separation() , scale)
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
        print("MULTI-POSITIONS:")
        print(particles[0].position)
        print(particles[1].position)
        print("scaling:",separation, scale, separation < scale)
        print("scaling2:",  (particles[1].position - particles[0].position).length_squared(), scale ** 2)
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
        
        print("REL POS:", rel_position)
        return self.deltas_to_update_binary(particles, rel_position, rel_velocity)
    
    def expand_binary(self, particles, scale, receeding = False):
        """
        Returns the change in positions and velocities for 
        the two-body system consisting of 'particle1' and 'particle2'.
        After applying the change the particles will lie
        close to distance 'scale' of one another.  
        The particles will be moving towards each other.
        """
        
        separation = (particles[1].position - particles[0].position).length_squared()     
        print("separation:", separation, scale, scale**2,separation > scale**2)
        if separation > scale**2:
            return particles.position * 0, particles.velocity * 0
            
        
        self.kepler_code.initialize_from_particles(particles)

        true_anomaly,   mean_anomaly = self.kepler_code.get_angles()
        semimajor_axis, eccentricity = self.kepler_code.get_elements()
        
        periapsis, apoapsis = self.get_periapsis_and_apoapsis(semimajor_axis, eccentricity)
        # largest distance minus 1% of the distance between peri and apo
        limit = apoapsis - 0.01*(apoapsis-periapsis)
        print("limit:", scale ** 2, scale, limit, apoapsis, periapsis, eccentricity)
        # we cannot scale to larger than the apoapsis distance
        # but if eccentricity > 1 we can!
        # TURNED OFF TO COMPARE WITH MULTIPLES:
        # if eccentricity <= 1 and scale > limit:
        if scale > limit and receeding:
            return particles.position * 0, particles.velocity * 0
        if scale > limit:
            scale = limit
        print("INPUT:", as_vector_quantity(self.kepler_code.get_separation_vector()))
        self.move_binary(scale, true_anomaly, receeding)
    
        rel_position = as_vector_quantity(self.kepler_code.get_separation_vector())
        rel_velocity = as_vector_quantity(self.kepler_code.get_velocity_vector())
        
        print("REL POS:", as_vector_quantity(self.kepler_code.get_separation_vector()))
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
        print("DELTA 1:" , delta_positions[0])
        print("DELTA 2:" , delta_positions[1])
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
        
        tree = trees.BinaryTreeOnParticle(particle)
        
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
        print("scale_particles_to_sphere", radius)
        center_of_mass_position = particles.center_of_mass()
        center_of_mass_velocity = particles.center_of_mass_velocity()

        print(center_of_mass_position)
        print(center_of_mass_velocity)

        particles.position -= center_of_mass_position
        particles.velocity -= center_of_mass_velocity
       
        print(particles.position)
        print(particles.velocity)
        # special case, 1 body
        if len(particles) == 1:
            "The position and velocity of this particle must be zero"
            
            tree = trees.ChildTreeOnParticle(particles[0])
            children = tree.get_children_subset()
            if len(children) == 2:
                "special case, scale the binary"
                scale = 2 * radius
                print("scale:", scale)
                delta_p, delta_v = self.kepler_orbits.compress_binary(children, scale, receeding = True)
                print(delta_p, delta_v) 
                
                for particle, dp, dv in zip(children, delta_p, delta_v):
                    self.move_particle(particle, dp, dv + center_of_mass_velocity)
                particles.velocity += center_of_mass_velocity
            return
            
        kinetic_energy = particles.kinetic_energy()
        potential_energy = particles.potential_energy(G = self.G)
        particle0, particle1 = self.get_particles_with_minimum_separation(particles)
        sphere_radius = particles.position.lengths().max()
        
        distance = (particle0.position - particle1.position).length()
        sum_of_radii = particle0.radius + particle1.radius
        separation = distance - sum_of_radii
        
            
        # special case, 2 bodies, we can use kepler to 
        # do the scaling in a consistent, energy preserving way
        if len(particles) == 2:
            if distance < sum_of_radii:
                scale = max(2*radius, sum_of_radii)
                delta_p, delta_v = self.kepler_orbits.expand_binary(particles, scale, receeding = True)
            elif separation > radius:
                scale = max(2 * radius, sum_of_radii)
                delta_p, delta_v = self.kepler_orbits.compress_binary(particles, scale, receeding = True)
            else:
                print("AA:",separation, 2 * radius,sum_of_radii)
                delta_p, delta_v = self.kepler_orbits.expand_binary(particles, 2 * radius, receeding = True)
            for particle, dp, dv in zip(particles, delta_p, delta_v):
                self.move_particle(particle, dp, dv + center_of_mass_velocity)
            return
        
        
        # for all other situations, we revert to scaling
        # where we preserve energy by scaling
        # the velocities
        print("DD:", distance, sum_of_radii, distance < sum_of_radii, radius, distance < 2 * radius)
        # we need to scale up, as the separation between particles is less than zero
        if distance < sum_of_radii:
            # use the largest scaling factor
            factor_position = max(sum_of_radii / distance, (2 * radius) / distance)
            
        # we need to scale up, as the minimum distance is less than the sphere diameter
        elif distance < 2 * radius:
            factor_position = (2.0 * radius) / distance
        # we need to scale down, the minimum distance is larger than the radius
        else:
            # we have room to scale down
            if distance > sum_of_radii:
                if sum_of_radii > (0.0 * radius):
                    factor_position = sum_of_radii / distance 	 	 
                else: 	 	 
                    factor_position = (2.0 * radius) / distance 
            # we have no room available for any scaling
            else:
                factor_position = 1.0
                
        factor_velocity_squared = 1.0 - (1.0/factor_position-1.0) * potential_energy/kinetic_energy
        if factor_velocity_squared < 0.0:
            from amuse.units import units
            print(particles.position)
            print(particles.velocity)
            print(particles.radius)
            print("radius", radius)
            print("distance", distance.as_quantity_in(units.AU))
            print("sum_of_radii", sum_of_radii)
            print("factor_position", factor_position)
            raise Exception("cannot scale the velocities")
        
        print(particles)
        factor_velocity = numpy.sqrt(factor_velocity_squared)
        delta_position = factor_position*(particles.position-center_of_mass_position) - particles.position
        delta_velocity = center_of_mass_velocity + factor_velocity*(particles.velocity-center_of_mass_velocity) - particles.velocity
        
        for particle, dp, dv in zip(particles, delta_position, delta_velocity):
            self.move_particle(particle, dp, dv + center_of_mass_velocity)
        #print "MINIMUM:",(2 * radius) , sum_of_radii, (particle0.position - particle1.position).length()

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
        
        return super(Binaries, self).remove_particles_from_store(keys)
        
    def get_children_subset(self, binaries, particle):
        return binaries._private.singles._subset(keys = (particle.child1.key, particle.child2.key))
        
        
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
    """
    
    Data model:
    
    1) particles -> multiples (binaries, ternaries etc.) + singles
        are evolved by gravity code
    2) multiples ->  subset of particles with components
        have a list of components
    3) component_singles -> singles part of a multiple 
        are part of a multiple, are stored relative
        to the center of mass position and velocity of the multiple
    4) binaries  -> separate list of particles
        have 2 components (in component_singles list)        
    
    """
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
        self.number_of_collisions = 0
        self.must_handle_one_encounter_per_stopping_condition = True
    
    def reset(self):
    
        self.components_of_multiples = Particles()
        self.multiples = Particles()
        self.singles   = Particles()
        
        self.particles = ParticlesSuperset(
            [self.singles, self.multiples],
            index_to_default_set = 0
        )
        
        self.singles_in_binaries = Particles()
        self.binaries  = Binaries(self.singles_in_binaries)
        self.singles_in_binaries_previous  = None
        
        
        self.gravity_code.reset()
        self.stopping_condition = self.gravity_code.stopping_conditions.collision_detection
        self.stopping_condition.enable()
        
        self.channel_from_code_to_model = self.gravity_code.particles.new_channel_to(self.particles)
        self.channel_from_model_to_code = self.particles.new_channel_to(self.gravity_code.particles)
        
        self.multiples_external_tidal_correction = zero
        self.multiples_internal_tidal_correction = zero
        self.multiples_integration_energy_error = zero
        self.all_multiples_energy = zero
        
        self.stopping_conditions.disable()
    
    def commit_particles(self):
        if len(self.multiples) == 0:
            if not len(self.binaries) == 0:
                for binary in self.binaries:
                    multiple = self.multiples.add_particle(binary)
                    components = self.components_of_multiples.add_particles(binary.components())
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
        
        #if len(self.singles) == 0:
        #    self.singles.add_particles(self.particles)
        # relink te components so these are in the right set
        if len(self.multiples) > 0:
            for x in self.multiples:
                x.components = x.components.get_intersecting_subset_in(self.components_of_multiples)
                
        #if len(self.particles) == 0:
        #    self.particles.add_particles(self.singles)
        #    self.particles.add_particles(self.multiples)
        self.gravity_code.particles.add_particles(self.particles)
        
        
        self.singles_in_binaries_previous = self.singles_in_binaries.copy()
        self.all_multiples_energy = self.get_total_energy_of_all_multiples()
        
        
    def evolve_model(self, time):
        self.stopping_conditions.unset()
        
        attributes_to_update = ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'radius']
        self.channel_from_model_to_code.copy_attributes(attributes_to_update)
        self.particles.synchronize_to(self.gravity_code.particles)
        
        self.model_time = self.gravity_code.model_time
        
        previous_time = None
        while self.model_time < time:
            self.gravity_code.evolve_model(time)
            self.model_time = self.gravity_code.model_time
            self.channel_from_code_to_model.copy()
            
            if self.stopping_condition.is_set():
                
                LOG_ENCOUNTER.info("found collision at time: {0}".format(self.gravity_code.model_time))
                initial_energy = self.gravity_code.get_total_energy()
                
                self.gravity_code.synchronize_model()
                self.channel_from_code_to_model.copy()
                self.handle_stopping_condition()
                self.particles.synchronize_to(self.gravity_code.particles)
                for i,k in enumerate(self.gravity_code.particles.key):
                    print(i, k)
                self.channel_from_model_to_code.copy()
                
                final_energy = self.gravity_code.get_total_energy()
                
                self.update_energy_bookkeeping(initial_energy, final_energy)
                                
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
        
    def get_total_energy_of_all_multiples(self):
        result = zero
        for x in self.multiples:
            result += self.get_energy_of_a_multiple(x)
        return result
        
    def get_energy_of_a_multiple(self, multiple):
        components = multiple.components
        tree = components.new_binary_tree_wrapper()
        singles = Particles()
        for node in tree.iter_descendant_leafs():
            singles.add_particle(node.particle)
        
        # tree is stored in rest state,
        # no energy of central particle
        energy  = singles.kinetic_energy()
        energy += singles.potential_energy(G = self.G)
        
        return  energy
        
    def update_energy_bookkeeping(self, initial_energy, final_energy):
        dE_gravity_code = final_energy - initial_energy
                
        self.local_energy_error = (
            dE_gravity_code 
            - self.handle_encounter_code.delta_energy
            - self.handle_encounter_code.delta_potential_in_field
        )
        self.internal_local_energy_error = (
              dE_gravity_code
            + self.handle_encounter_code.delta_multiple_energy
            - self.handle_encounter_code.delta_potential_in_field
        )
        self.corrected_internal_local_energy_error = (
              dE_gravity_code 
            + self.handle_encounter_code.delta_multiple_energy
            - self.handle_encounter_code.delta_potential_in_field
            + self.handle_encounter_code.delta_internal_potential
            - self.handle_encounter_code.scatter_energy_error
        )
        
        LOG_ENERGY.info('net local error = {0}'.format(self.local_energy_error))
        LOG_ENERGY.info('net local internal error = {0}'.format(self.internal_local_energy_error))
        LOG_ENERGY.info('corrected local internal error = {0}'.format(self.corrected_internal_local_energy_error))
        
        self.multiples_external_tidal_correction += self.handle_encounter_code.delta_potential_in_field
        self.multiples_internal_tidal_correction -= self.handle_encounter_code.delta_internal_potential
        self.multiples_integration_energy_error += self.handle_encounter_code.scatter_energy_error
        
        self.all_multiples_energy = self.get_total_energy_of_all_multiples()
        self.total_energy = final_energy + self.all_multiples_energy
        self.corrected_total_energy = (
            self.total_energy
            - self.multiples_external_tidal_correction
            - self.multiples_internal_tidal_correction
            - self.multiples_integration_energy_error
        )
        
        LOG_ENERGY.info('total energy (top+mul) = {0}'.format(self.total_energy))
        LOG_ENERGY.info('corrected_total energy (top+mul) = {0}'.format(self.corrected_total_energy))
        
        
    @property
    def all_singles(self):
        result = self.singles.copy()
        for multiple in self.multiples:
            components = multiple.components
            tree = components.new_binary_tree_wrapper()
            subset = Particles()
            for node in tree.iter_descendant_leafs():
                subset.add_particle(node.particle)
            subset.position += multiple.position
            subset.velocity += multiple.velocity
            delattr(subset, 'child1')
            delattr(subset, 'child2')
            result.add_particles(subset)
        return result
        
    def handle_stopping_condition(self):
        encounters = self.determine_encounters()
        
        new_binaries = Particles()
        dissolved_binaries = Particles()
        updated_binaries = Particles()
        new_multiples = Particles()
        dissolved_multiples = Particles()
        updated_multiples = Particles()
        
        for particles_in_encounter in encounters:
            self.handle_encounter(
                particles_in_encounter,
                new_binaries,
                dissolved_binaries,
                updated_binaries,
                new_multiples,
                dissolved_multiples,
                updated_multiples
            )
    
        if self.stopping_conditions.multiples_change_detection.is_enabled():
            if len(new_multiples) > 0 or len(dissolved_multiples) > 0 or len(updated_multiples) > 0:
                self.stopping_conditions.multiples_change_detection.set(
                    new_multiples,
                    dissolved_multiples,
                    updated_multiples
                )
                
        if self.stopping_conditions.binaries_change_detection.is_enabled():
            if len(new_binaries) > 0 or len(dissolved_binaries) > 0 or len(updated_binaries) > 0 :
                self.stopping_conditions.binaries_change_detection.set(
                    new_binaries.get_intersecting_subset_in(self.binaries),
                    dissolved_binaries,
                    updated_binaries.get_intersecting_subset_in(self.binaries)                    
                )
                
    def handle_encounter(
            self, 
            particles_in_encounter,
            new_binaries,
            dissolved_binaries,
            updated_binaries,
            new_multiples,
            dissolved_multiples,
            updated_multiples
        ):
        code = self.handle_encounter_code
        code.reset()
        before = particles_in_encounter.copy()
        LOG_ENCOUNTER.info("found encounter with particles {0}".format(particles_in_encounter.key))
        code.particles_in_encounter.add_particles(particles_in_encounter)
        print(self.particles , particles_in_encounter)
        code.particles_in_field.add_particles(self.particles - particles_in_encounter)
        code.existing_binaries.add_particles(self.binaries)
        code.existing_multiples.add_particles(self.multiples)
        LOG_ENCOUNTER.info("handling encounter, {0} particles in encounter".format(len(code.particles_in_encounter)))
        code.execute()
        LOG_ENCOUNTER.info(
            "handling encounter, finished, {0} new multiples, {1} dissolved multiples, {2} updated multiples".format(
                len(code.new_multiples),
                len(code.dissolved_multiples),
                len(code.updated_multiples
            ))
        )
        
        
        new_multiples.add_particles(code.new_multiples)
        
        dissolved_multiples.add_particles(code.dissolved_multiples)
        for x in dissolved_multiples:
            x.components = x.components.copy()
            
        updated_multiples.add_particles(code.updated_multiples)
        for x in updated_multiples:
            x.components = x.components.copy()
        
        
        LOG_ENCOUNTER.info("captured singles: {0}".format(code.captured_singles.key))
        LOG_ENCOUNTER.info("released singles: {0}".format(code.released_singles.key))
        # update the singles (will have singles and multiples)
        self.singles.remove_particles(code.captured_singles)
        self.singles.add_particles(code.released_singles)
        
        LOG_ENCOUNTER.info("dissolved multiples: {0}".format(code.dissolved_multiples.key))
        LOG_ENCOUNTER.info("new multiples: {0}".format(code.new_multiples.key))
        # update multiples
        self.multiples.remove_particles(code.dissolved_multiples)
        for x in code.dissolved_multiples:
            self.components_of_multiples.remove_particles(x.components)
        
        new_multiples = self.multiples.add_particles(code.new_multiples)
        for x in new_multiples:
            x.components = self.components_of_multiples.add_particles(x.components)
        
            
        self.number_of_collisions += 1
        #io.write_set_to_file((code.all_particles_in_encounter, code.particles_after_encounter, code.particles_before_scaling), "encounter-{0}.h5".format(self.number_of_collisions), "amuse", names=('before', 'after', 'after_smalln'), version="2.0", append_to_file=False)
        # code.all_particles_in_encounter
        # update binaries
        for x in code.dissolved_binaries:
            self.singles_in_binaries.remove_particle(x.child1)
            self.singles_in_binaries.remove_particle(x.child2)
        for x in code.new_binaries:
            self.singles_in_binaries.add_particle(x.child1)
            self.singles_in_binaries.add_particle(x.child2)
        for x in code.updated_binaries:
            child1 =  x.child1.as_particle_in_set(self.singles_in_binaries)
            child1.position = x.child1.position
            child1.velocity = x.child1.velocity
            child2 =  x.child2.as_particle_in_set(self.singles_in_binaries)
            child2.position = x.child2.position
            child2.velocity = x.child2.velocity
            
        self.binaries.remove_particles(code.dissolved_binaries)
        self.binaries.add_particles(code.new_binaries)
        
        self.singles_in_binaries_previous = self.singles_in_binaries.copy()
    
        if self.stopping_conditions.encounter_detection.is_enabled():
            model = Particles()
            particles_before_encounter = Particles()
            particles_before_encounter.add_particles(code.all_particles_in_encounter)
            particles_after_encounter = Particles()
            particles_after_encounter.add_particles(code.particles_after_encounter)
            
            particle = Particle()
            particle.particles_before_encounter = particles_before_encounter
            particle.particles_after_encounter = particles_after_encounter
            model.add_particle(particle)
            self.stopping_conditions.encounter_detection.set(model)
        
        channel = code.updated_binaries.new_channel_to(self.binaries)
         
        new_binaries.add_particles(code.new_binaries) 
        dissolved_binaries.add_particles(code.dissolved_binaries) 
        updated_binaries.add_particles(code.updated_binaries) 
        
        if 0:
            print("before:", particles_in_encounter)
            self.particles.remove_particles(particles_in_encounter)
            print("after:", code.particles_after_encounter.key)
            self.particles.add_particles(code.particles_after_encounter)
        elif 0:
            self.gravity_code.particles.remove_particles(particles_in_encounter)
        channel = code.particles_after_encounter.new_channel_to(self.particles)
        channel.copy_attributes(["x","y","z", "vx", "vy","vz"])
        
    
    def determine_encounters(self):
        particles0 = self.stopping_condition.particles(0).copy()
        particles1 = self.stopping_condition.particles(1).copy()
        if self.must_handle_one_encounter_per_stopping_condition:
            particles0 = particles0[:1]
            particles1 = particles1[:1]
            
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
        
    def get_total_energy(self):
        
        self.total_energy = self.gravity_code.get_total_energy() + self.all_multiples_energy
        self.corrected_total_energy = (
            self.total_energy
            - self.multiples_external_tidal_correction
            - self.multiples_internal_tidal_correction
            - self.multiples_integration_energy_error
        )
        return self.corrected_total_energy

    def update_model(self):
        
        # we do all the work on the singles_in_binaries_previous set
        # this makes sure all keys match attributes
        self.singles_in_binaries_previous.delta_position = self.singles_in_binaries.position - self.singles_in_binaries_previous.position
        self.singles_in_binaries_previous.delta_velocity = self.singles_in_binaries.velocity - self.singles_in_binaries_previous.velocity
        
        # take not yet updated positions
        channel = self.components_of_multiples.new_channel_to(self.singles_in_binaries_previous)
        channel.copy_attributes(['x', 'y', 'z', 'vx', 'vy', 'vz'])
        
        # update these
        self.singles_in_binaries_previous.position += self.singles_in_binaries_previous.delta_position
        self.singles_in_binaries_previous.velocity += self.singles_in_binaries_previous.delta_velocity
        
        # copy back
        channel = self.singles_in_binaries_previous.new_channel_to(self.components_of_multiples)
        channel.copy_attributes(['x', 'y', 'z', 'vx', 'vy', 'vz'])
        
        channel = self.singles_in_binaries.new_channel_to(self.components_of_multiples)
        channel.copy_attribute('mass')
        
        for multiple in self.multiples:
            components = self.get_singles_of_a_multiple(multiple)
            multiple.mass = components.mass.sum()
            center_of_mass = components.center_of_mass()
            center_of_mass_velocity = components.center_of_mass_velocity()
            multiple.position += center_of_mass
            multiple.velocity += center_of_mass_velocity
            components.position -= center_of_mass
            components.velocity -= center_of_mass_velocity
    
        attributes_to_update = ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'radius']
        channel = self.particles.new_channel_to(self.gravity_code.particles)
        channel.copy_attributes(attributes_to_update)
        self.singles_in_binaries_previous = self.singles_in_binaries.copy()
        
    def get_singles_of_a_multiple(self, multiple):
        components = multiple.components
        tree = components.new_binary_tree_wrapper()
        singles = Particles()
        for node in tree.iter_descendant_leafs():
            singles.add_particle(node.particle)
        return singles
