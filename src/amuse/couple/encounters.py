"""
This module defines the classe to handle handle close
encounters between particles. 

It is used by the multiples module.
"""

from amuse.datamodel import Particles
from amuse.units import constants
from amuse.units.quantities import as_vector_quantity

#codes to use
from amuse.community.kepler.interface import Kepler

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
    HARD_BINARY_FACTOR=1
    
    def __init__(self, 
        particles_in_encounter, 
        particles_in_field = None, 
        existing_multiples = None, 
        existing_binaries = None,
        kepler_orbits = None,
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
    
        if kepler_orbits is None:
            self.kepler_orbits = KeplerOrbits()
            
    def start(self):
        
        self.determine_scale_of_particles_in_the_encounter()
        
        self.select_neighbours_from_field()
        
        self.determine_singles_from_particles_and_neighbours_in_encounter()
        
        self.determine_initial_sphere_of_singles_in_encounter()
        
        self.move_all_singles_to_initial_sphere_frame_of_reference()
        
        self.evolve_singles_in_encounter_until_end_state()
        
        self.determine_structure_of_the_evolved_state()
        
        self.scale_evolved_state_to_initial_sphere()
        
        self.remove_soft_binaries_from_evolved_state()
        
        self.move_evolved_state_to_original_frame_of_reference()
        
        self.determine_multiples_in_the_evolved_state()
        
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
        
    def determine_singles_from_particles_and_neighbours_in_encounter(self):
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
        
        near_distance = self.large_scale_of_particles_in_the_encounter * self.NEIGHBOURS_FACTOR
        near_particles = self.particles_in_field[distances <= near_distance]
        
        self.particles_close_to_encounter.add_particles(near_particles)
                
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
            key = (binary.child1.key,binary.child2.key)
            binary_lookup_table[key] = binary
        
        # a branch in the tree is a node with two children
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
        
class KeplerOrbits(object):
    
    def __init__(self, nbody_converter = None):
        self.kepler_code = Kepler(nbody_converter)
        self.kepler_code.initialize_code()
        
    def get_semimajor_axis_and_eccentricity_for_binary_components(self, particle1, particle2):
        
        particles = Particles()
        particles.add_particle(particle1)
        particles.add_particle(particle2)
        
        self.kepler_code.initialize_from_particles(particles)
        
        return self.kepler_code.get_elements()

    def compress_binary_components(self, particle1, particle2, scale):
        """
        Returns the change in positions and velocities for 
        the two-body system consisting of 'particle1' and 'particle2'.
        After applying the change the particles will lie
        inside distance 'scale' of one another.  
        The final orbit will be receding (moving away from each other).
        """
        total_mass = particle1.mass + particle2.mass
        rel_position = particle1.position - particle2.position
        rel_velocity = particle1.velocity - particle2.velocity
        separation = rel_position.length()
        
        particles = Particles()
        particles.add_particle(particle1)
        particles.add_particle(particle2)
        
        center_of_mass_position = particles.center_of_mass()
        center_of_mass_velocity = particles.center_of_mass_velocity()
        positions_to_center_of_mass = particles.position - center_of_mass_position
        velocities_to_center_of_mass = particles.velocity - center_of_mass_velocity
        
        if separation <= scale: 
            # particles are already close together, no scaling done
            # AVE is this correct, will the particle(s) be receding?
            #     or should some movement always happen
            return particles.position * 0, particles.velocity * 0
        
        
        self.kepler_code.initialize_from_dyn(
            total_mass,
            rel_position[0], rel_position[1], rel_position[2],
            rel_velocity[0], rel_velocity[1], rel_velocity[2]
        )

        true_anomaly, mean_anomaly = self.kepler_code.get_angles()
        semimajor_axis, eccentricity = self.kepler_code.get_elements()
        
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
            # apoapsis is infinity but this is better
            # as we only use it for the limit


        # closest distance plus 1% of the distance between peri and apo
        limit = periapsis + 0.01*(apoapsis-periapsis)
        
        # we cannot scale to smaller than the periapsis distance
        if scale < limit:
            scale = limit
            
        if true_anomaly < 0:
            self.kepler_code.advance_to_periastron()
            self.kepler_code.advance_to_radius(scale)
        else:
            if self.kepler_code.get_separation() < scale:
                self.kepler_code.advance_to_radius(scale)
            else:
                self.kepler_code.return_to_radius(scale)
        
        # Note: Always end up on an outgoing orbit.  If
        # periastron > scale, we are now just past periapsis.

        rel_position = as_vector_quantity(self.kepler_code.get_separation_vector())
        rel_velocity = as_vector_quantity(self.kepler_code.get_velocity_vector())
        
        f = particle2.mass / total_mass 
        fractions = numpy.asarray([-f, (1-f)]).reshape(2,1)
        
        delta_positions  = (rel_position * fractions) - positions_to_center_of_mass
        delta_velocities = (rel_velocity * fractions) - velocities_to_center_of_mass
        
        return delta_positions, delta_velocities
    
class ScaleSystem(object):
    
    def __init__(self, kepler_orbits, G = constants.G):
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
    
    def minimum_separation(self, particles):
        positions = particles.position
        radii = particles.radius

        result = None
        for i in range(len(particles) - 1):
            i_position = positions[i]
            j_positions = positions[i+1:]
            
            i_radius = radii[i+1:]
            j_radii = radii[i+1:]
            
            delta_positions = i_position - j_positions
            dr = delta_positions.lengths()
            sum_radii = i_radius + j_radii
            
            delta = dr - sum_radii
            if result is None:
                result = delta.min()
            else:
                result = min(result, delta.min())
        return result
            


    
    def scale_particles_to_sphere(particles, radius):
        """
        Rescale the system of particles to lie within a sphere
        of the given radius.
        System may be compressed or expanded.
        """
        

        center_of_mass_position = particles.center_of_mass()
        center_of_mass_velocity = particles.center_of_mass_velocity()

        particles.position -= center_of_mass_position
        particles.velocity -= center_of_mass_velocity
        
        kinetic_energy = particles.kinetic_energy()
        potential_energy = particles.potential_energy(G = self.G)
        minimum_separation = self.minimum_separation(particles)
        sphere_radius = position.lengths().max()
        
        
        for i in range(len(node_list)):
            m = 1
            rad = node_list[i].radius.number
            posi = node_list[i].position
            pos = (posi-cmpos).number
            vel = (node_list[i].velocity-cmvel).number
            kin += m*numpy.inner(vel,vel)
            dpot = 0.0
            for j in range(i+1,len(node_list)):
                mj = node_list[j].mass.number
                radj = node_list[j].radius.number
                dposj = (node_list[j].position-posi).number
                rij = math.sqrt(numpy.inner(dposj,dposj))
                if sepmin > rij-rad-radj:
                    radsum = rad + radj
                    imin = i
                    jmin = j
                    sepmin = rij - radsum
                    rijmin = rij
                dpot -= mj/math.sqrt(numpy.inner(dposj,dposj))
            pot += m*dpot
        size = math.sqrt(size)
        kin /= 2

        #fac = 0.5*scale.number/size	# scale to radius
        #fac = scale.number/rijmin		# scale to distance
        fac = radsum/rijmin			# scale to zero separation

        # Compress (or expand) the system and increase (or decrease) the
        # velocities (relative to the center of mass) to preserve the
        # energy.  If fac > 1, expansion is always OK if E > 0, which it
        # should be at this point (but check anyway...).  May have E < 0
        # if we have a system with small negative energy, stopped because
        # it is too big.

        vfac2 = 1-(1/fac-1)*pot/kin
        if vfac2 < 0:
            print "Can't expand top level system to rjmin > ri+rj"
            print "fac =", fac, " pot =", pot, " kin =", kin
            sys.stdout.flush()
            f = pot/(kin+pot)
            vfac2 = 0.0
        vfac = math.sqrt(vfac2)
        if fac > 0.0:
            for n in node_list:
                n.position = cmpos + fac*(n.position-cmpos)
                n.velocity = cmvel + vfac*(n.velocity-cmvel)
