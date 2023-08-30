import sys
import numpy
import collections
import math
import copy

from amuse.datamodel import particle_attributes
from amuse.datamodel import trees
from amuse.datamodel import Particle, Particles 
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units.quantities import zero
from amuse.support.exceptions import KeysNotInStorageException
from amuse import io

#---------------------------------------------------------------------
#
# Steve's ToDo list of features to be added/improved in the multiples
# module.
#
# 1. Should use only perturbers (within ~100 x interaction scale) in
# computing the tidal energy change, not the entire system.
#
# 2. If the gravity module supports it, only perturbers need be
# synchronized before and reinitialized after the interaction.
#
# 3. Including near neighbors in a 2-body interaction is likely to
# lead to spurious binaries, since the 3- (or more-) body interaction
# will be followed to completion, when in fact it should be stopped
# when the neighbor has left the interaction region.  The result is
# that binaries may form prematurely.  If we want to include
# neighbors, we should also implement item 4 below, to allow long
# interactions to be broken into pieces.  Including neighbors in the
# interaction may also lead to problematic final configurations and
# large internal/external tidal errors.  Alternatively, we can let
# near neighbors simply veto the encounter, moving the work back into
# the gravity module, until a "clean" 2-body scattering can be
# identified.  Vetoing is now the default.
#
# 4. We should seek a better prescription for compressing 3-body and
# higher-order configurations.  Currently we conserve energy, but not
# angular momentum.
#
# 5. There is no provision for physical collisions in the smallN code,
# and no logic in the Multiples module to manage stars having both
# dynamical and physical radii.
#
#---------------------------------------------------------------------

# The following simple CM indexing scheme is OK for N < 1000000.  An
# improved scheme might be desirable, but it must be compatible with
# the integer IDs used in most gravity modules.

root_index = 1000000

def new_root_index():
    global root_index
    root_index += 1
    return root_index

def name_object(tree):
    name = "{ "
    if hasattr(tree, "child1"):
        children = [tree.child1, tree.child2]
    else:
        #children = tree.get_tree_subset()
        children = [tree.particle.child1, tree.particle.child2]
    for child in sorted(children, \
                        key=lambda x: x.id if hasattr(x ,"id") else 1e10):
        if child.child1 is not None:
            name += name_object(child)
        else:
            name += str(child.id)
        name += " "
    name += "}"
    return name

def name_pair(comp1, comp2):
    return '('+str(comp1.id)+','+str(comp2.id)+')'

known_roots = {}
def assign_id_to_root(tree):

    # Determine the object's description, then search to see if we
    # know about it.  If we do, return that ID, otherwise create a new
    # ID.

    global known_roots
    my_name = name_object(tree)
    if my_name in known_roots.keys():
        return known_roots[my_name]
    else:
        new_root_id = new_root_index()
        known_roots[my_name] = new_root_id
        return new_root_id
    
def is_a_parent(child1_key, child2_key):
    return child1_key > 0 or child2_key > 0

def is_not_a_child(is_a_child):
    return is_a_child == 0

def get_component_binary_elements(comp1, comp2, kep, peri = 0):

    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()		# per unit reduced mass, note
    if peri:
        M,th = kep.get_angles()
        if M < 0:
            kep.advance_to_periastron()
        else:
            kep.return_to_periastron()
    t = kep.get_time()

    return mass,a,e,r,E,t

def get_cm_binary_elements(p, kep, peri = 0):
    return get_component_binary_elements(p.child1, p.child2, kep, peri)

class DidNotFinishException(Exception):
    pass

class Multiples(object):
    
    def __init__(self, 
                 gravity_code,
                 resolve_collision_code_creation_function,
                 kepler_code, 
                 gravity_constant = None, **options):

        # Codes to use.
        
        self.gravity_code = gravity_code
        self.resolve_collision_code_creation_function \
            = resolve_collision_code_creation_function
        self.kepler = kepler_code

        # Data structures.

        # Local copy of CM (= root) particles in the gravity code.
        self._inmemory_particles = self.gravity_code.particles.copy()

        # Dictionary connecting center of mass particles with the
        # multiple tree structure lying under them.
        #
        # Syntax:
        # 	root_to_tree[center_of_mass] = binary_tree
        #
        # where center_of_mass is a Particle and binary_tree is created by
        # trees.BinaryTreesOnAParticleSet.
        
        self.root_to_tree = {}

        '''
        # Unecessary code since copy_attribute below does the same thing.
        if len(self.gravity_code.particles) == 0:
            self._inmemory_particles.id = None
        else:
            self._inmemory_particles.id = self.gravity_code.particles.index_in_code
        '''

        self._inmemory_particles.child1 = None
        self._inmemory_particles.child2 = None
        self.channel_from_code_to_memory = \
            self.gravity_code.particles.new_channel_to(self._inmemory_particles)
        self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

        # FLASH interface needs a channel the other way also - Josh.
        self.channel_from_memory_to_code = \
            self._inmemory_particles.new_channel_to(self.gravity_code.particles)

        if gravity_constant is None:		# default is N-body units
            gravity_constant = nbody_system.G
        
        self.gravity_constant = gravity_constant

        # Energy bookkeeping.
        zero_energy = zero * self.gravity_code.kinetic_energy
        self.multiples_external_tidal_correction = zero_energy
        self.multiples_internal_tidal_correction = zero_energy
        self.multiples_integration_energy_error = zero_energy

        # Count the number of collisions we found for comparing with
        # encounters algorithm.
        self.number_of_collisions = 0
    
        # Repeat encounter management data
        self.old_star_1 = 0
        self.old_star_2 = 0
        self.repeat_count = 0

        # The following tunable parameters govern the multiples logic:

        # Nominal size of the top-level encounter, relative to the sum
        # of the radii of the interacting components (no veto) or
        # their separation (veto).
        
        #self.initial_scale_factor = 1.0
        self.initial_scale_factor = 2.0

        # Perturbation above which to include a neighbor, estimated
        # using the neighbor distance and the initial separation of
        # the top-level two-body encounter.  (Previously was a simple
        # distance criterion...)
        
        #self.neighbor_distance_factor = 1.0
        #self.neighbor_distance_factor = 2.0
        self.neighbor_perturbation_limit = 0.02

        # Neighbor veto policy.  True means we allow neighbors to veto
        # a two-body encounter (meaning that don't want to deal with
        # complex initial many-body configurations).  False means we
        # include neighbors in the multiple integration.

        #self.neighbor_veto = False
        self.neighbor_veto = True

        # Size of the rescaled final system, relative to the initial
        # scale.  Should be 1 + epsilon.
        
        self.final_scale_factor = 1.01

        # Initial separation for the scattering experiment, relative
        # to the initial scale.  May be limited by apocenter in case
        # of a bound top-level interaction.
        
        self.initial_scatter_factor = 10.0

        # Final separation for the scattering experiment, relative to
        # the initial scattering scale (meaning 2 x 10 times the
        # initial encounter scale, by default).  May be limited by
        # binary properties in case of a bound top-level interaction.

        #self.final_scatter_factor = 10.0
        self.final_scatter_factor = 2.0

        # Binary retention policy.  Retain a binary if its apocenter
        # (True) or 2*semi-major axis (False) is less than the
        # dynamical radius of its CM.  False is the more conservative
        # choice.
        
        self.retain_binary_apocenter = True

        # Maximum allowed perturbation at apocenter on a wide binary.
        self.wide_perturbation_limit = 0.01
        
        # Turn on/off global debugging output level (0 = no output, 1
        # = minimal, 2 = normal debug, 3 = verbose debug).
        
        self.global_debug = 1

        # Turn on debugging in the encounters code.
        self.debug_encounters = False

        # Turn on/off experimental code to check tidal perturbation.
        self.check_tidal_perturbation = False

    @property
    def particles(self):
        return self.gravity_code.particles
    
    @property 
    def total_mass(self):
        return self.gravity_code.total_mass
    
    @property 
    def kinetic_energy(self):
        return self.gravity_code.kinetic_energy
    
    @property 
    def potential_energy(self):
        return self.gravity_code.potential_energy
    
    @property 
    def parameters(self):
        return self.gravity_code.parameters
        
    @property 
    def model_time(self):
        return self.gravity_code.model_time
        
    @property
    def stars(self):
        result = self._inmemory_particles.copy()
        for root, tree in self.root_to_tree.items():
            root_particle = root.as_particle_in_set(self._inmemory_particles)
            result.remove_particle(root)

            # This returns a pointer to the actual leaves - Josh.

            leaves = tree.get_leafs_subset()
      
            original_star = tree.particle

            dx = root_particle.x - original_star.x
            dy = root_particle.y - original_star.y
            dz = root_particle.z - original_star.z
            dvx = root_particle.vx - original_star.vx
            dvy = root_particle.vy - original_star.vy
            dvz = root_particle.vz - original_star.vz

            # Note that here we add leaves to another particle set,
            # and this becomes its own deep copy of leaves.  So
            # changes to leaves_in_result have no effect on leaves -
            # Josh.

            leaves_in_result = result.add_particles(leaves)
            leaves_in_result.x += dx
            leaves_in_result.y += dy
            leaves_in_result.z += dz
            leaves_in_result.vx += dvx
            leaves_in_result.vy += dvy
            leaves_in_result.vz += dvz

        return result

    def update_leaves_pos_vel(self):

        # The FLASH interface needs a function that updates the
        # properties of the leaves from the properties of the
        # particles in the gravity code - Josh.
    
        # NOTE: Unlike multiples.stars, this actually moves the real
        # leaves, not a copy of the leaf particles.  So tree.particle
        # also then needs to be updated - Josh.

        local_debug = False
        self.channel_from_code_to_memory.copy() # update the copy in memory
                                                # from the gravity code - Josh

        for root, tree in self.root_to_tree.items():
            root_particle = root.as_particle_in_set(self._inmemory_particles)

            leaves = tree.get_leafs_subset()
            original_star = tree.particle

            if (local_debug):
                old_leaves_x = leaves.x

                print("In update_leaves_pos_vel before update.")
                print("Tree pos =", tree.particle.position.in_(units.cm))
                print("Root pos =", root.position.in_(units.cm))
                print("Leaf pos =", leaves.position.in_(units.cm))

            dx = root_particle.x - original_star.x
            dy = root_particle.y - original_star.y
            dz = root_particle.z - original_star.z
            dvx = root_particle.vx - original_star.vx
            dvy = root_particle.vy - original_star.vy
            dvz = root_particle.vz - original_star.vz

            leaves.x  += dx
            leaves.y  += dy
            leaves.z  += dz
            leaves.vx += dvx
            leaves.vy += dvy
            leaves.vz += dvz

            # Update the original particle info stored in
            # tree.particle - Josh.

            original_star.x = root_particle.x
            original_star.y = root_particle.y
            original_star.z = root_particle.z

            original_star.vx = root_particle.vx
            original_star.vy = root_particle.vy
            original_star.vz = root_particle.vz

            if (local_debug):
                new_leaves = tree.get_leafs_subset()
                leaves_dx  = leaves.x - old_leaves_x

                if (leaves_dx[0].number == 0.0):
                    print("These leaves aren't moving!")
                elif (leaves_dx[0].number == dx[0].number):
                    print("These leaves arrived precisely when they meant to!")
                else:
                    print("I have no idea what these damn leaves are doing!")
                    print("leaves_dx =", leaves_dx)
                    print("dx =", dx)

            if (local_debug):
                print("In update_leaves_pos_vel after update.")
                print("Tree pos =", tree.particle.position.in_(units.cm))
                print("Root pos =", root.position.in_(units.cm))
                print("Leaf pos =", leaves.position.in_(units.cm))

        return

    def create_binary(self, star1, star2):

        # Experimental code to include a binary directly into the
        # multiples database.

        M,a,e,r,E,tperi = get_component_binary_elements(star1, star2, 
                                                        self.kepler, 1)

        binary = Particles()
        binary.add_particle(star1)
        binary.add_particle(star2)
        
        cm = Particle(mass=M,
                      position=binary.center_of_mass(),
                      velocity=binary.center_of_mass_velocity())
        binary.add_particle(cm)

        binary.child1 = None
        binary.child2 = None
        cm = binary[2]
        cm.child1 = binary[0]
        cm.child2 = binary[1]
        
        set_radii(binary, self.kepler, self.global_debug)

        self.gravity_code.particles.remove_particle(star1)
        self.gravity_code.particles.remove_particle(star2)

        # Complete the bookkeeping.
        
        tree = trees.BinaryTreesOnAParticleSet(binary,
                                               "child1", "child2")

        for t in tree.iter_binary_trees():	# only one tree...
            t.particle.id = assign_id_to_root(t)
            self.gravity_code.particles.add_particle(t.particle)
            self.root_to_tree[t.particle] = t.copy()
            print('\nCreated binary from', star1.id, 'and', star2.id, \
                  ' CM =', t.particle.id)
            print('M =', M, ' a =', a, ' e =', e, ' E =', E)

        self.gravity_code.particles.synchronize_to(self._inmemory_particles)
        self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    def check_trees(self):

        # Print out some debugging information on multiples in the system.

        print('')
        print('check_trees:', len(self.root_to_tree), 'tree(s)')
        for root, tree in self.root_to_tree.items():
            print(root.position)			# current
            print(tree.particle.position)	# original
            leaves = tree.get_leafs_subset()	# components (original)
            print(leaves.center_of_mass())
        print('')

    def get_gravity_at_point(self, radius, x, y, z):
        return self.gravity_code.get_gravity_at_point(radius, x, y, z)
    
    def get_potential_at_point(self, radius, x, y, z):
        return self.gravity_code.get_potential_at_point(radius, x, y, z)
    
    def commit_particles(self):
        return self.gravity_code.commit_particles()
    
    def get_time(self):
        return self.gravity_code.get_time()
    
    def get_total_energy(self, code):
        try:
            binaries_energy = code.get_binary_energy()	# include binaries if
        except:						# the code understands
            binaries_energy = zero
        total_energy = code.potential_energy + code.kinetic_energy \
	                + binaries_energy

        return total_energy

    #--------------------------------------------------------------
    # Note that the true total energy of a multiple isn't quite the
    # Emul returned below, since the tidal potential of components on
    # one another is not taken into account.

    def get_total_multiple_energy(self):	# uses kepler
        Nbin = 0
        Nmul = 0
        Emul = zero
        for x in self.root_to_tree.values():	# loop over top-level trees
            Nmul += 1
            nb,E = get_multiple_energy(x, self.kepler)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul

    # This version returns the true total energy of all multiples.

    def get_total_multiple_energy2(self):
        Nbin = 0
        Nmul = 0
        Emul = zero
        for x in self.root_to_tree.values():	# loop over top-level trees
            Nmul += 1
            nb,E = get_multiple_energy2(x, self.gravity_constant)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul

    def print_multiples(self):					# uses kepler

        # Print basic information on all multiples in the system,
        # using the root_to_tree database.  This version uses
        # print_multiple_simple() to format the output.

        if self.global_debug > 0:
            for x in self.root_to_tree.values():
                print_multiple_simple(x, self.kepler)

    def print_multiples2(self, pre, kT, dcen):	# uses kepler

        # Print information on all multiples in the system, using the
        # root_to_tree database.  This version uses
        # print_multiple_detailed() to format the output, and returns
        # the numbers and energies of multiples found.

        Nbin = 0
        Nmul = 0
        Emul = zero
        for x in self.root_to_tree.values():
            Nmul += 1
            nb,E = print_multiple_detailed(x, self.kepler, pre, kT, dcen)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul
            
    def print_trees_summary(self):
        if len(self.root_to_tree) > 0:
            print('number of multiples:', len(self.root_to_tree))
            sys.stdout.flush()

    def evolve_model(self, end_time, callback=None):

        stopping_condition = \
            self.gravity_code.stopping_conditions.collision_detection
        #stopping_condition.enable()  # allow user to set this; don't override
        
        time = self.gravity_code.model_time
        print("\nmultiples: evolve model to", end_time, "starting at", time)
        sys.stdout.flush()

        count_resolve_encounter = 0
        count_ignore_encounter = 0
        
        while time <= end_time:		# the <= here allows zero-length steps

            if self.global_debug > 1:
                print('')
                print('calling evolve_model from', \
	    	    self.gravity_code.model_time, 'to', end_time)
                sys.stdout.flush()

            self.gravity_code.evolve_model(end_time)
            newtime = self.gravity_code.model_time
           
            # JB modified this: in Bonsai we can take a zero-length
            # time step to detect multiples. That would cause the
            # newtime == time to evaluate to true when there are
            # multiples detected and break out of the evaluate loop
            # before the time reached end_time.  Same is now possible
            # with ph4 (SLWM, 6/18).

            if newtime == time and (stopping_condition.is_set() == False):
                break
            
            time = newtime

            if stopping_condition.is_set():

                # An encounter has occurred.  Synchronize all stars in
                # the gravity code.  We synchronize everything for
                # now, but it would be better to just synchronize
                # neighbors if gravity_code supports that.  TODO
                
                self.gravity_code.synchronize_model()
                
                star1 = stopping_condition.particles(0)[0]
                star2 = stopping_condition.particles(1)[0]
                ignore = 0
                self.before = Particles()
                self.after = Particles()
                self.after_smalln = Particles()
                #print 'self.before:', self.before

                # Note from Steve, 8/12: We can pick up a lot of
                # encounters that are then ignored here.  I have
                # (temporarily?) duplicated this check in the ph4
                # module (jdata.cc).

                r = (star2.position-star1.position).length()
                v = (star2.velocity-star1.velocity).length()
                vr = ((star2.velocity-star1.velocity) \
                        * (star2.position-star1.position)).sum()

                EPS = 0.001
                if True or vr < EPS*r*v:    # True ==> keep all encounters
                                            # returned by gravity_code

                    if self.global_debug > 1:
                        print('\n'+'~'*60)
                    elif self.global_debug > 0:
                        print('')                        
                    if self.global_debug > 0:
                        print('interaction at time', time)
                
                    # As with synchronize above, we should only copy
                    # over data for the interacting particles and
                    # their neighbors.  TODO

                    self.channel_from_code_to_memory.copy()
                    self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                    
                    initial_energy = self.get_total_energy(self.gravity_code)

                    star1 = star1.as_particle_in_set(self._inmemory_particles)
                    star2 = star2.as_particle_in_set(self._inmemory_particles)

                    cont = True
                    if callback != None:
                        cont = callback(time, star1, star2)

                    if self.global_debug > 0:
                        print('initial top-level:',         \
                            star1.id, '('+str(star1.radius)+')', \
                            star2.id, '('+str(star2.radius)+')')
                    if self.global_debug > 1:
                        print('                   r =', r)
                        print('                   v =', v)
                        print('                   v.r =', vr)
                    sys.stdout.flush()

                    # Do the scattering.

                    veto, dE_top_level_scatter, dphi_top, dE_mul, \
                        dphi_int, dE_int, final_particles \
                        = self.manage_encounter(time, star1, star2, 
                                                self._inmemory_particles,
                                                self.gravity_code.particles,
                                                self.kepler)

                    if cont and not veto:

                        # Recommit is done automatically and reinitializes all
                        # particles.  Later we will just reinitialize a list if
                        # gravity_code supports it. TODO
                        
                        self.gravity_code.particles.synchronize_to(
                            self._inmemory_particles)		# sets _inmemory = gravity
                        self.channel_from_code_to_memory.copy_attribute(
                            "index_in_code", "id")

                        final_energy = self.get_total_energy(self.gravity_code)
                        dE_top_level = final_energy - initial_energy

                        # Local bookkeeping:
                        #
                        #	dE_top_level is the actual energy
                        #	change in the top-level gravity system
                        #	due to this encounter
                        #
                        #	dE_top_level_scatter is the change in
                        #	top-level internal energy of the
                        #	scattering system
                        #
                        #	dphi_top is the top-level tidal error
                        #	(currently unabsorbed) due to the
                        #	change in configuration of the
                        #	scattering system in the top-level
                        #	tidal field
                        #
                        #	dE_mul is the change in stored
                        #	multiple energy associated with the
                        #	encounter
                        #
                        #	dphi_int is the internal tidal energy
                        #	error due to configuration changes in
                        #	the scattering system
                        #
                        #	dE_int is the integration error in the
                        #	scattering calculation
                        #
                        # We *always* expect
                        #
                        #   dE_top_level - dE_top_level_scatter - dphi_top = 0.
                        #
                        # If this is not the case, then there is an
                        # error in the internal bookkeeping of
                        # manage_encounter().

                        if self.global_debug > 2:
                            #print 'top-level initial energy =', initial_energy
                            #print 'top-level final energy =', final_energy
                            print('dE_top_level =', dE_top_level)
                            print('dE_top_level_scatter =', dE_top_level_scatter)
                            print('dphi_top =', dphi_top)
                            print('dphi_int =', dphi_int)
                            print('dE_int =', dE_int)
                            print('dE_top_level-dE_top_level_scatter-dphi_top =',\
                                dE_top_level - dE_top_level_scatter - dphi_top)

                        if self.global_debug > 2:
                            print('net local error =', \
                                  dE_top_level - dE_top_level_scatter - dphi_top)
                            print('scatter integration error =', dE_int)

                        # We also expect
                        #
                        #	dE_top_level_scatter + dE_mul
                        #			= dphi_top + dE_int - dphi_int.
                        #
                        # Monitor this and keep track of the
                        # cumulative value of the right-hand side of
                        # this equation.

                        if self.global_debug > 2:
                            print('dE_mul =', dE_mul)
                            print('internal local error =', \
                                  dE_top_level + dE_mul - dphi_top)
                            print('corrected internal local error =', \
                            	  dE_top_level + dE_mul - dphi_top \
                                  + dphi_int - dE_int)

                        self.multiples_external_tidal_correction += dphi_top
                        self.multiples_internal_tidal_correction -= dphi_int
                        self.multiples_integration_energy_error += dE_int
                        
                        # Doing this energy calculation at every
                        # encounter is expensive when dealing with
                        # hundreds of binaries or more.  It is clearly
                        # problematic when building a whole system of
                        # binaries.
                        
                        #Nmul, Nbin, Emul = self.get_total_multiple_energy2()

                        # Global bookkeeping:
                        #
                        # We don't absorb individual tidal or
                        # integration errors, but instead store their
                        # totals in
                        #
                        #	    self.multiples_external_tidal_correction,
                        #	    self.multiples_internal_tidal_correction,
                        # and
                        #	    self.multiples_inegration_energy_error.
                        #
                        # Then
                        #
                        #   E(top-level) + Emul
                        #	    - self.multiples_external_tidal_correction
                        #	    - self.multiples_internal_tidal_correction
                        #	    - self.multiples_integration_energy_error
                        #
                        # should be constant.  Any non-conservation
                        # represents an error in bookkeeping or
                        # algorithm design.
                        '''
                        print 'total energy (top+mul) =', \
                            final_energy + Emul
                        print 'corrected total energy =', \
                            final_energy + Emul \
                                - self.multiples_external_tidal_correction \
                                - self.multiples_internal_tidal_correction \
                                - self.multiples_integration_energy_error
                        '''
                        
                        # Print info on all multiples associated with
                        # the current interaction.

                        if self.global_debug > 1:
                            for x in final_particles:
                                if hasattr(x, "child1") \
                                   and not (getattr(x, "child1") is None):
                                    print_multiple_simple(
                                        trees.BinaryTreeOnParticle(x),
                                        self.kepler)

                        count_resolve_encounter += 1

                    else:
                        ignore = 1

                    if self.global_debug > 1:
                        print('~'*60)
                    sys.stdout.flush()

                else:
                    ignore = 1
        
                self.number_of_collisions += 1

                '''
                io.write_set_to_file((self.before,
                                     self.after, self.after_smalln),
                                     "multiples-{0}.h5".format(self.number_of_collisions),
                                     "amuse", names=('before', 'after', 'after_smalln'),
                                     version="2.0", append_to_file=False)
                '''
                count_ignore_encounter += ignore

        print('')
        print('Resolved', count_resolve_encounter, 'encounters')
        print('Ignored', count_ignore_encounter, 'encounters')
        sys.stdout.flush()

        self.gravity_code.synchronize_model()
        self.channel_from_code_to_memory.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    def expand_encounter(self, scattering_stars):

        # Create an encounter particle set from the top-level stars.
        # Add stars to the encounter set, add in components when we
        # encounter a binary/multiple.

        particles_in_encounter = Particles(0)
        Emul = zero

        for star in scattering_stars:
            if star in self.root_to_tree:
                tree = self.root_to_tree[star]
                isbin, dEmul = get_multiple_energy2(tree, self.gravity_constant)
                Emul += dEmul
                openup_tree(star, tree, particles_in_encounter)
                del self.root_to_tree[star]
            else:
                particles_in_encounter.add_particle(star)

        return particles_in_encounter, Emul

    def manage_encounter(self, global_time, star1, star2,
                         stars, gravity_stars, kep):

        # Manage an encounter between star1 and star2.  Stars is the
        # python memory dataset (_inmemory_particles).  Gravity_stars
        # is the gravity code data (only used to remove the old
        # components and add new ones).  On entry, stars and
        # gravity_stars should contain the same information.  Return
        # values are the change in top-level energy, the tidal error,
        # and the integration error in the scattering calculation.
        # Steps below follow those defined in the PDF description.

        # print 'in manage_encounter'
        # sys.stdout.flush()

        # Record the state of the system prior to the encounter, in
        # case we need to abort and return without making any changes.
        #
        # 'gravity_stars' is no longer included in the snapshot below,
        # as it would make N individual calls to get_state...

        snapshot = {
            'global_time': global_time,
            'star1': star1.copy(),
            'star2': star2.copy(),
            'stars': stars.copy(),
            #'gravity_stars': gravity_stars.copy(),
            'self.root_to_tree': self.root_to_tree.copy(),
            'particles_in_encounter': Particles(0),
            'scattering_stars': Particles(0)
        }
        snapshot['particles_in_encounter'].add_particle(snapshot['star1'])
        snapshot['particles_in_encounter'].add_particle(snapshot['star2'])

        # find_binaries(stars, self.gravity_constant)

        # self.check_trees()

        #----------------------------------------------------------------
        # 1a. Build a list of stars involved in the scattering.  Start
        # with star1 and star2.

        scattering_stars = Particles(particles = (star1, star2))
        star1 = scattering_stars[0]
        star2 = scattering_stars[1]
        center_of_mass = scattering_stars.center_of_mass()
        other_stars = stars - scattering_stars	# probably only need perturbers?
        
        # Brewer Mod:  Check for a repeat encounter.

        if (star1.id == self.old_star_1 and star2.id == self.old_star_2) \
           or (star1.id == self.old_star_2 and star2.id == self.old_star_1):
            self.repeat_count += 1
        else:
            self.repeat_count = 0
            self.old_star_1 = star1.id
            self.old_star_2 = star2.id

        # 1b. Add neighbors if desired.  Use a perturbation criterion.
        # Also impose a simple neighbor veto, if specified.  Start by
        # sorting all stars by perturbation on the CM.  Later, use
        # neighbors only, if supported.  TODO

        sep12 = ((star1.position-star2.position)**2).sum().sqrt()
        rad12 = star1.radius + star2.radius

        # Sep12 is the separation of the two original components.  It
        # should be slightly less than the sum of their radii, rad12,
        # but it may be much less in unexpected circumstances or if
        # vetoing is in effect.  Initial_scale sets the "size" of the
        # interaction and the distance to which the final products
        # will be rescaled.  Rad12 also ~ the 90 degree scattering
        # distance for two stars, and hence the natural limit on
        # binary scale.

        if not self.neighbor_veto:
            initial_scale = self.initial_scale_factor * rad12
        else:
            initial_scale = self.initial_scale_factor * sep12

        if self.global_debug > 1:
            print('initial_scale =', initial_scale)

        # The basic sort on other_stars is by perturbation, not
        # distance.  Maintain sorted lists of stars, distances (d),
        # and perturbations (actually m/d**3).

        distances = (other_stars.position - center_of_mass).lengths()
        pert = other_stars.mass / distances**3
        indices = numpy.argsort(-pert.number)	# decreasing sort
        sorted_stars = other_stars[indices]
        sorted_distances = distances[indices]
        sorted_perturbations = pert[indices]
        fac12 = 0.5*(star1.mass + star2.mass)/sep12**3

        largest_perturbers = []
        if self.check_tidal_perturbation and len(sorted_stars) > 0:
            
            if self.global_debug > 1:
                print("sorted_stars", sorted_stars[:5])
                print("sorted_distances", sorted_distances[:5])
                print("sorted_perturbations", sorted_perturbations[:5]/fac12)

            max_pert = sorted_perturbations[0]/fac12
            largest_perturbers = [sorted_stars[0]]
            
            # This should be replaced with something faster using
            # numpy, like:
            #
            # largest_perturbers = sorted_stars[np.greater(sorted_perturbations,
            #     0.025*sorted_perturbations[0])] - Josh.

            for i in range(1, len(sorted_stars)):
                if sorted_perturbations[i] > 0.025*sorted_perturbations[0]:
                    largest_perturbers.append(sorted_stars[i])

        # Perturbation limit for identification as a neighbor.
        
        pert_min = self.neighbor_perturbation_limit*fac12
        for i in range(len(sorted_stars)):    # NB no loop if len() = 0
            
            star = sorted_stars[i]

            # Include anything lying "inside" the binary, even if it
            # is a weak perturber.

            if sorted_perturbations[i] > pert_min \
                    or sorted_distances[i] < sep12:
                if not self.neighbor_veto:
                    scattering_stars.add_particle(star)
                    if self.global_debug > 1:
                        print('added', end=' ')
                        if hasattr(star, 'id'):
                            print('star', star.id, end=' ')
                        else:
                            print('unknown star', end=' ')
                        print('to scattering list')
                        sys.stdout.flush()
                    snapshot['scattering_stars'].add_particle(star)
                    #initial_scale = sorted_distances[i]    # don't expand!
                else:
                    if self.global_debug > 0:
                        print('encounter vetoed by', \
                            star.id, 'at distance', \
                            sorted_distances[i], \
                            'pert =', sorted_perturbations[i]/fac12)
                    if self.repeat_count > 0: self.repeat_count -= 1
                    return True, 0., 0., 0., 0., 0., None

        self.before.add_particles(scattering_stars)

        # Note: sorted_stars, etc. are used once more, when checking
        # for wide binaries (at 6b below).

        #----------------------------------------------------------------
        # 2a. Calculate the total internal and external potential
        # energy of stars to remove from the gravity system, using the
        # potential of scattering_stars relative to the other
        # top-level objects in the stars list (later: just use
        # neighbors TODO).

        # Terminology from the PDF description:

        E0 = scattering_stars.kinetic_energy() \
            + scattering_stars.potential_energy(G=self.gravity_constant)
        phi_rem = potential_energy_in_field(scattering_stars,         
                                            stars-scattering_stars,
                                            G=self.gravity_constant)

        if self.global_debug > 2:
            print('E0 =', E0)
            print('phi_rem =', phi_rem)

        # 2b. If there are no neighbors, separate star1 and star2 to
        #     some larger "scattering" radius.  If neighbors exist,
        #     just start the "scattering" interaction in place.

        # First define some basic properties of the top-level
        # interaction.

        M,a,e,r,E,tperi = get_component_binary_elements(star1, star2, 
                                                        self.kepler, 1)

        Etop = E*star1.mass*star2.mass/M
        ttrans = self.gravity_constant*M/(4*abs(E))**1.5

        # Note: transit time = 0.056 * period for a bound orbit.

        if e < 1:
            peri = a*(1-e)
            apo = a*(1+e)
            period = self.kepler.get_period()
        else:
            peri = a*(e-1)
            apo = 1.e9*a	# 1.e9 is large but otherwise arbitrary
            period = 1.e9*ttrans

        initial_scatter_scale = self.initial_scatter_factor * initial_scale

        # Limit initial_scatter_scale (rescale_binary_components will
        # impose a limit, but good to have the limit available at this
        # level).

        if initial_scatter_scale > 0.9*apo:
            initial_scatter_scale = 0.9*apo

        if len(scattering_stars) == 2:
            #print "rescaling in:", (star1.position - star2.position).length()
            rescale_binary_components(star1, star2, kep,
                                      initial_scatter_scale, compress=False)
            #print "rescaling out:", (star1.position - star2.position).length()

        # 2c. Remove the interacting stars from the gravity module.
        
        for s in scattering_stars:
            gravity_stars.remove_particle(s)

        #----------------------------------------------------------------
        # 3. Create a particle set to perform the close encounter
        #    calculation.

        # Note this has to delete the root_to_tree particle in
        # multiples as we have no idea what the end product of the
        # encounter will be.  So deleting in expand_encounter really
        # is a feature, not a bug.  Don't mess with it! - Josh
        
        particles_in_encounter, Emul_init \
                = self.expand_encounter(scattering_stars)

        # Terminology from the PDF description:

        E1 = particles_in_encounter.kinetic_energy() + \
             particles_in_encounter.potential_energy(G=self.gravity_constant)

        dphi_1 = E1 - E0 - Emul_init

        if self.global_debug > 2:
            print('E1 =', E1)
            print('Emul_init =', Emul_init)
            print('dphi_1 =', dphi_1)

        #----------------------------------------------------------------
        # 4. Run the small-N encounter in the center of mass frame.
     
        total_mass = scattering_stars.mass.sum()

        cmpos = scattering_stars.center_of_mass()
        cmvel = scattering_stars.center_of_mass_velocity()
        particles_in_encounter.position -= cmpos
        particles_in_encounter.velocity -= cmvel

        #E1CM = particles_in_encounter.kinetic_energy() + \
        #       particles_in_encounter.potential_energy(G=self.gravity_constant)
        #print 'E1 (CM) =', E1CM

        # Relevant available length and time scales:
        #
        # Encounter:
        #	sep12 =	 actual separation
        #	rad12 =	 sum of radii (should be ~b90)
        #
        # Top-level orbit:
        #	a =	 orbital semimajor axis
        #	peri =	 orbital periastronn
        #	apo =	 orbital apastron
        #	tperi =	 time to pericenter
        #	period = orbital period, if defined
        #	ttrans = transit time
        #
        # Resonance:
        #	rvir =	 viral length scale
        #	tvir =	 virial time scale

        rvir = self.gravity_constant*M/(4*abs(E1/M))
        tvir = self.gravity_constant*M/(4*abs(E1/M))**1.5

        if self.global_debug > 2:
            print('Encounter:')
            print('    sep12 =', sep12)
            print('    rad12 =', rad12)
            print('Top-level:')
            print('    E/mu =', E)
            print('    Etop =', Etop)
            print('    M =', M)
            print('    semi =', a)
            print('    ecc =', e)
            print('    peri =', peri)
            print('    apo =', apo)
            print('    tperi =', tperi)
            print('    ttrans =', ttrans)
            print('    period =', period)
            print('Resonance:')
            print('    rvir =', rvir)
            print('    tvir =', tvir)
        else:
            if self.global_debug > 0:
                print('M =', M, ' Etop =', Etop)
            if self.global_debug > 1:
                print('a =', a, ' e =', e, ' P =', period)

        sys.stdout.flush()

        # The original concept of this module was to follow the
        # encounter as an isolated scattering experiment until it is
        # cleanly resolved.  In this case, the bookkeeping and
        # post-encounter logic are straightforward.  We expect that a
        # clean resolution always eventually occurs, but for a complex
        # interaction this may take a long time.  In addition, the
        # long time scales and large excursions of the intermediate
        # orbit may render the physicality of the scattering approach
        # questionable.

        # The alternative approach pursued below is to try to define
        # limiting length and time scales for the encounter, based on
        # the initial configuration.  Encounters exceeding these limits
        # will be returned to the large-N simulation, possibly to be
        # picked up again later.  This approach leads to significant
        # bookkeeping issues, and the "clean" original concept should
        # always be retained as a fallback option.

        # We need a reliable time scale to set end_time and delta_t
        # for the scattering interaction.  It is possible that we pick
        # up an encounter very close to periastron, so tperi may not
        # be useful.  With reasonable limits on encounter size and a
        # treatment of quasi-stable systems, a time limit on the
        # smallN integration may not be needed, but impose some
        # reasonable value here, just in case.  Note that we have to
        # deal with the possible consequences in resolve_collision().

        # Note that the current check for quasistability requires that
        # the system configuration remain unchanged for 10 outer
        # orbital periods, and is not yet reliable.  TODO

        # If the encounter is a flyby, then the relevant scales are
        # the orbital semimajor axis and transit time (*10, say).  We
        # don't want to follow a wide bound system onto a second orbit
        # unless the size of the orbit is less than a few times the 90
        # degree turnaround distance.  If the encounter is a
        # resonance, then the relative scales are the virial radius
        # (*10) and virial time scale (*100).  If it is bound but wide
        # (and likely a flyby), then the relevant scales are the the
        # orbital semimajor axis and period (*10, say).

        # Also set a limit on the minimum scale, in case of retries.

        end_time = max(2*abs(tperi), 10*ttrans, 100*tvir)
        if E.number < 0: end_time = max(end_time, 10*period)

        delta_t = max(1.5*abs(tperi), tvir)

        if self.global_debug > 1:
            print('end_time =', end_time)
            print('delta_t =', delta_t)

        # Note: radii used here should really be based on
        # perturbation, not simply distance.  TODO

        orbit_scale = 2*a
        if E.number < 0: orbit_scale = 1.1*a*(1+e) # a*(1+0.9*e)

        if self.global_debug > 2:
            print('orbit_scale =', orbit_scale)

        # Final_scatter_scale is the scale at which we will terminate
        # the smallN integration. This is a guess of the scale where,
        # if the system exceeds it, the interacting particles can be
        # decomposed into well separated pieces that can be returned
        # to the N=body code, even if the encounter isn't over.

        final_scatter_scale \
            = max(self.final_scatter_factor * initial_scatter_scale,
                  orbit_scale, 10*rvir)

        # Limit the scatter scale in case of a very wide orbit.

        if orbit_scale > 2*initial_scatter_scale \
           	and final_scatter_scale > orbit_scale:
            final_scatter_scale = orbit_scale

        min_scatter_scale = 2*initial_scale	# never go below this value
        if min_scatter_scale >= 0.5*final_scatter_scale:
            final_scatter_scale = 2*min_scatter_scale

        # The integration ends when any particle is more than
        # final_scatter_scale from the CM of the system (hence the
        # factor of 2).  RECONSIDER - needs a mass scale factor, and
        # still OK for a wide orbit?  TODO

        final_scatter_scale /= 2
        min_scatter_scale /= 2

        if self.global_debug > 1:
            print('final_scatter_scale =', final_scatter_scale)
            print('min_scatter_scale =', min_scatter_scale)

        # NOTE: to revert to the original concept, simply set
        # final_scatter_scale and end_time to very large values.

        if 0:
            print('particles in encounter:')
            print('position:', particles_in_encounter.position)
            print('velocity:', particles_in_encounter.velocity)

        try:
            scatter_energy_error \
                = self.resolve_collision(particles_in_encounter,
                                         final_scatter_scale,
                                         min_scatter_scale,
                                         end_time, delta_t)

        except DidNotFinishException:

            # In this case, simply abort the encounter and continue
            # the main simulation.

            print("*** SmallN encounter did not finish. ", \
                  "Aborting and returning to top level.")

            global_time = snapshot['global_time']
            star1 = snapshot['star1']
            star2 = snapshot['star2']
            stars = snapshot['stars']
            #gravity_stars = snapshot['gravity_stars']
            gravity_stars.add_particle(star1)
            gravity_stars.add_particle(star2)
            gravity_stars.add_particles(snapshot['scattering_stars'])
            self.root_to_tree = snapshot['self.root_to_tree']
            zero_en = 0.0 * E0

            return False, zero_en, zero_en, zero_en, zero_en, zero_en, \
                snapshot['particles_in_encounter']

        # Note that on return, particles_in_encounter contains CM
        # nodes in the list.

        E2CM = get_energy_of_leaves(particles_in_encounter,
                                    G=self.gravity_constant)
        Etop = particles_in_encounter.kinetic_energy() \
             + particles_in_encounter.potential_energy(G=self.gravity_constant)
        if self.global_debug > 1:
            print('E2 (CM) =', E2CM)

        particles_in_encounter.position += cmpos
        particles_in_encounter.velocity += cmvel

        # Terminology from the PDF description:

        E2 = get_energy_of_leaves(particles_in_encounter,
                                  G=self.gravity_constant)
        dE_int = E2 - E1	# should equal scatter_energy_error
        err = (dE_int-scatter_energy_error)/max(E1,E2)
        if abs(err) > 1.e-12:
            if self.global_debug > 0:
                print('*** warning: dE_int mismatch ***')
                if self.global_debug > 1:
                    print('scatter_energy_error =', scatter_energy_error)
                    print('dE_int =', dE_int)
                    #print particles_in_encounter
                    print('E1 =', E1, 'E2 =', E2)

        if self.global_debug > 2:
            print('E2 =', E2)
            print('scatter_energy_error =', scatter_energy_error)
            print('dE_int =', dE_int)

        #----------------------------------------------------------------
        # 5a. Identify multiple structure after the encounter.  First
        #     create an object to handle the new binary information.
        
        #Brewer Mod:  Create the appropriate COM particle for the pseudo-binary
        '''
        if self.repeat_count > 9:
            print "repeat encounter detected; forcing binary creation"
            pseudoCOM = Particles(1)
            pseudoCOM.child1 = star1
            pseudoCOM.child2 = star2
            pseudoCOM.mass = star1.mass + star2.mass
            pseudoCOM.position = cmpos
            pseudoCOM.velocity = cmvel
            pseudoCOM.radius = star1.radius + star2.radius
            print particles_in_encounter
            print pseudoCOM
            particles_in_encounter.add_particles_in_store(pseudoCOM)
        '''
        #End Mod section.

        binaries = trees.BinaryTreesOnAParticleSet(particles_in_encounter,
                                                   "child1", "child2")

        # 5b. Compress the top-level nodes before adding them to the
        #     gravity code.  Also recompute the external potential and
        #     optionally absorb the tidal error into the top-level
        #     nodes of the encounter list.  Finally, add the change in
        #     top-level energy of the interacting subset into dEmult,
        #     so E(ph4) + dEmult should be conserved.

        # Single stars.
        stars_not_in_a_multiple = binaries.particles_not_in_a_multiple()

        #print 'stars_not_in_a_multiple:'
        #print stars_not_in_a_multiple

        # Multiple centers of mass.
        roots_of_trees = binaries.roots()

        #----------------------------------------------------------------
        # 6a. Scale to a radius slightly larger than the initial one.
        # Rescaling does just that -- neither computes nor attempts to
        # absorb the tidal error.  If we want to absorb the tidal
        # error rather than simply recording it, do so after splitting
        # wide binaries below.  TODO

        final_scale = self.final_scale_factor * initial_scale
        
        # Note that stars_not_in_a_multiple and roots_of_trees are
        # simply convenient partitions of particles_in_encounter.
        # They are pointers into the underlying particle set.

        scale_top_level_list(stars_not_in_a_multiple,
                             roots_of_trees,
                             self.kepler,
                             final_scale,
                             self.gravity_constant,
                             self.global_debug)

        # 6b. Break up wide top-level binaries.  Do this after
        #     rescaling because we want to preserve binary binding
        #     energies.  Also place the wide binaries at pericenter to
        #     minimize the tidal error.

        # Number of top-level nodes.
        lt = len(stars_not_in_a_multiple) + len(roots_of_trees)

        # Roots to be deleted after the loop.
        roots_to_remove = []

        for root in roots_of_trees:
            comp1 = root.child1
            comp2 = root.child2

            mass,semi,ecc,r,E,t = \
                get_component_binary_elements(comp1, 
                                              comp2, 
                                              self.kepler)
            apo = semi*(1+ecc)
            if self.retain_binary_apocenter:
                binary_scale = apo
            else:
                binary_scale = 2*semi	    # (the more conservative choice)

            # Estimate the maximum perturbation on this binary due to
            # its current strongest external perturber.

            max_perturbation = 0.0
            if len(sorted_perturbations) > 0:
                max_perturbation = \
                    	2*sorted_perturbations[0]*binary_scale**3/mass
                perturber = sorted_stars[0]
                perturber_distance = sorted_distances[0]
            
            # Check that other stars involved in the encounter but not
            # in this multiple are not the dominant perturbation.

            stars_to_check = Particles()
            for t in binaries.iter_binary_trees():
                if t.particle != root:		# exclude self interaction
                    stars_to_check.add_particles(t.get_leafs_subset())

            #while len(roots_to_check) > 0:
            #    r = roots_to_check.pop()
            #    if r != root:
            #        if hasattr(r, "child1"):
            #            if r not in roots_to_check:
            #                roots_to_check.append(r)
            #        else:
            #            stars_to_check.extend(r)
            try:
                stars_to_check.remove_particle(star1)
            except KeysNotInStorageException:
                #print 'failed to remove star1'
                pass
            try:
                stars_to_check.remove_particle(star2)
            except KeysNotInStorageException:
                #print 'failed to remove star2'
                pass

            # Check perturbation due to stars_to_check on root.

            for s in stars_to_check:
                distance = (s.position - root.position).length()
                pert = s.mass / distance**3
                s_perturbation = 2*pert*binary_scale**3/mass
                if self.global_debug > 1:
                    print("star %s, distance %s, pert %s, s_pert %s, max_pert %s" \
                        % (s.id, distance, pert, s_perturbation,
                           max_perturbation))
                if s_perturbation > max_perturbation:
                    max_perturbation = s_perturbation
                    perturber = s
                    perturber_distance = distance

            #if binary_scale > rad12:
            if max_perturbation < self.wide_perturbation_limit \
               or self.repeat_count > 9:
                if self.global_debug > 0:
                    print('accepting lightly perturbed or repeat binary', \
                    	name_pair(comp1,comp2))
                    if self.global_debug > 1:
                        print('    semi =', semi, 'E/mu =', E)
                        print('    apo =', apo, 'peri =', semi*(1-ecc))
                if max_perturbation > 0:
                    if self.global_debug > 1:
                        print('    strongest perturber is', perturber.id, \
                            'with apo perturbation', max_perturbation)
                        print('    nearest neighbor is', perturber.id, \
                            'at distance', perturber_distance)
                        print('    repeat_count =', self.repeat_count)
                else:
                    if max_perturbation > 0:
                        print('    perturbation = 0')
                self.repeat_count = 0		# probably unnecessary
                sys.stdout.flush()

            else:
                if self.global_debug > 0:
                    if max_perturbation > 0:
                        print('splitting perturbed binary', \
                            name_pair(comp1,comp2))
                if self.global_debug > 1:
                    print('    semi =', semi, 'E/mu =', E)
                    print('    apo =', apo, 'peri =', semi*(1-ecc))
                    print('    strongest perturber is', perturber.id, \
                          'with apocenter perturbation', max_perturbation)
                    print('    nearest neighbor is', perturber.id, \
                          'at distance', perturber_distance)
                sys.stdout.flush()

                # See the "special case" logic in
                # scale_top_level_list().  If this is a sole bound
                # top-level object, it has already been scaled to the
                # desired separation and should *not* be modified
                # here.  Otherwise, move the components past
                # periastron to initial_scatter_scale.

                if lt > 1:

                    # Could use rescale_binary_components() for this,
                    # but code here is more compact, since we have
                    # already initialized the kepler structure.

                    cmpos = root.position
                    cmvel = root.velocity
                    if self.global_debug > 1:
                        print('moving binary to periastron')
                    self.kepler.advance_to_periastron()
                    if self.global_debug > 1:
                        print('advancing binary to', final_scale)
                        sys.stdout.flush()
                    self.kepler.advance_to_radius(final_scale)
                    
                    dx = quantities.AdaptingVectorQuantity()
                    dx.extend(kep.get_separation_vector())
                    dv = quantities.AdaptingVectorQuantity()
                    dv.extend(kep.get_velocity_vector())

                    f1 = comp1.mass/root.mass
                    comp1.position = cmpos - f1*dx
                    comp1.velocity = cmvel - f1*dv
                    comp2.position = cmpos + (1-f1)*dx
                    comp2.velocity = cmvel + (1-f1)*dv
                        
                # Changing the list here would disrupt the loop
                # bookkeeping.  Remove any split-up roots after the
                # loop, then recalculate all data structures and
                # restart the loop.

                #particles_in_encounter.remove_particle(root)
                roots_to_remove.append(root)

                # Note that removing the root will reinstate the
                # children as top-level objects:

        if len(roots_to_remove) > 0:
            
            for r in roots_to_remove:
                particles_in_encounter.remove_particle(r)
    
            # Recompute the tree structure.

            binaries = \
                trees.BinaryTreesOnAParticleSet(particles_in_encounter,
                                                "child1", "child2")
            # Single stars.
            stars_not_in_a_multiple = binaries.particles_not_in_a_multiple()

            # Multiple centers of mass.
            roots_of_trees = binaries.roots()

        #----------------------------------------------------------------
        # 7. Add the new top-level nodes to the gravity module.

        top_level_nodes = stars_not_in_a_multiple + roots_of_trees

        # Terminology from the PDF description:

        KE3 = top_level_nodes.kinetic_energy()
        E3 = KE3 + top_level_nodes.potential_energy(G=self.gravity_constant)

        phi_ins = potential_energy_in_field(top_level_nodes, 
                                            stars - scattering_stars,
                                            G=self.gravity_constant)

        Emul_final = zero
        for tree in binaries.iter_binary_trees():            
            isbin, dEmul = get_multiple_energy2(tree, self.gravity_constant)
            Emul_final += dEmul

        dphi_2 = E2 - Emul_final - E3

        if self.global_debug > 2:
            print('E3 =', E3)
            print('phi_ins =', phi_ins)
            print('Emul_final =', Emul_final)
            print('dphi_2 =', dphi_2)

        # 7a. Set radii to reflect multiple structure.

        set_radii(particles_in_encounter, self.kepler, self.global_debug)

        # Print diagnostics on added particles. Strip dimensions
        # because of numpy problem noted below.

        if self.global_debug > 0:
            print('final top-level:', end=' ')
        r = zero
        v = zero
        vr = zero
        for i in top_level_nodes:
            if self.global_debug > 0:
                print(i.id, '('+str(i.radius)+')', end=' ')
            for j in top_level_nodes:
                if i.id > j.id:
                    rij = ((i.position-j.position)**2).sum().sqrt()
                    if rij > r:
                        r = rij
                        v = ((i.velocity-j.velocity)**2).sum().sqrt()
                        vr = ((j.velocity-i.velocity) \
                               * (j.position-i.position)).sum()

        if self.global_debug > 0:
            print('')
            print('M =', top_level_nodes.mass.sum(), end=' ')
            print('Etop =', Etop)
        if self.global_debug > 1 and len(top_level_nodes) > 1:
            print('                 r =', r)
            print('                 v =', v)
            print('                 v.r =', vr)
            #print 'top_level_nodes:'
            #print top_level_nodes
        sys.stdout.flush()

        # Update the gravity module with the new data.
        self.after.add_particles(stars_not_in_a_multiple)

        # 7b. Add stars not in a binary to the gravity code.
        if len(stars_not_in_a_multiple) > 0:
            #print 'adding stars_not_in_a_multiple:'
            #print stars_not_in_a_multiple
            gravity_stars.add_particles(stars_not_in_a_multiple)

        # 7c. Add the roots to the gravity code
        multiples_particles = Particles()
        multiples_particles.id = None

        for tree in binaries.iter_binary_trees():
            tree.particle.id = assign_id_to_root(tree)	# assign CM ID (was 0)
            #tree.particle.components = subset
            #print 'adding particle:'
            #print tree.particle
            gravity_stars.add_particle(tree.particle)
            self.after.add_particle(tree.particle)
            multiples_particles.add_particle(tree.particle)

        if self.global_debug > 1:            
            print("multiples: interaction products: singles:", \
                stars_not_in_a_multiple.id, "multiples: ", \
                multiples_particles.id) 
            
        # 7d. Store all trees in memory for later reference.

        # Note this is actually where the trees get added to the
        # multiples module, and is the appropriate place to modify any
        # of the leaves / roots in the module.  Also this is what as
        # getting deleted in a call to expand_encounter, but is
        # unmodified in a call to stars - Josh.
        
        for tree in binaries.iter_binary_trees():
            self.root_to_tree[tree.particle] = tree.copy()
            
        # Return enough information to monitor all energy errors.

        dE_top = E3 - E0
        dphi_top = phi_ins - phi_rem
        dEmul = Emul_final - Emul_init
        dphi_int = dphi_2 - dphi_1

        #-------------------------------------------------------
        # Flag (but don't yet correct) large tidal corrections.

        dph = dphi_top/KE3
        if abs(dph) > 1.e-2:		# 1.e-2 is small but otherwise arbitrary
            if self.global_debug > 0:
                print('*** tidal correction =', dph, 'KE ***')
            #print 'initial configuration: phi =', \
            #    potential_energy_in_field(scattering_stars, 
            #                              stars - scattering_stars,
            #                              G=self.gravity_constant)
            pminmin, fminmin, dxminmin \
                = find_nn2(scattering_stars, stars-scattering_stars,
                           self.gravity_constant)
            if pminmin != None:
                if self.global_debug > 1:
                    print('closest field/list pair is', \
                        str(fminmin.id)+'/'+str(pminmin.id), \
                        ' distance/scale =', dxminmin/initial_scale)
            #print 'final configuration: phi =', \
            #    potential_energy_in_field(top_level_nodes, 
            #                              stars - scattering_stars,
            #                              G=self.gravity_constant)
            pminmin, fminmin, dxminmin \
                = find_nn2(top_level_nodes, stars-scattering_stars,
                           self.gravity_constant)
            if pminmin != None:
                if self.global_debug > 1:
                    print('closest field/list pair is', \
                        str(fminmin.id)+'/'+str(pminmin.id), \
                        ' distance/scale =', dxminmin/initial_scale)
        #-------------------------------------------------------

        # Experimental code to try to correct external tidal errors.
        # Compare dphi_top with range of possible quadrupole
        # corrections due to closest perturber.  Start with the
        # simplest case.
        #
        #	tidal potential change is dphi_top
        #	multiple center of mass is cmpos
        #	perturbers are in largest_perturbers
        
        if self.check_tidal_perturbation \
            and len(particles_in_encounter) == 2 and len(top_level_nodes) == 2:

            print('checking quadrupole perturbations')

            # *** Retain unitless code for now (Steve, 4/18). ***
            
            m1 = top_level_nodes[0].mass
            m2 = top_level_nodes[1].mass
            dx = top_level_nodes[1].position - top_level_nodes[0].position
            x = (dx**2).sum().sqrt()
            print('x =', x, 'M =', m1+m2)

            for p in largest_perturbers:
                m3 = p.mass
                id = p.id
                dr = p.position - cmpos
                r = (dr**2).sum().sqrt()
                phi = -self.gravity_constant*M*m3/r
                dphiQ = -(self.gravity_constant*(m1*m2/M)*m3/r)*(x/r)**2
                print(' ', str(id)+':', 'r =', r, 'm =', p.mass, \
                      'dphi_top/dphiQ =', dphi_top/dphiQ)

        return False, dE_top, dphi_top, dEmul, dphi_int, dE_int, \
               particles_in_encounter

    def resolve_collision(self,
                          particles,
                          final_scatter_scale,
                          min_scatter_scale,
                          end_time,
                          delta_t):

        pre = 'encounter:'		# identifier for all output

        # Take the system described by particles and evolve it forward
        # in time until it is over.  Don't update global quantities,
        # don't interpret the outcome.  Return the energy error due to
        # the smallN integration.

        if self.debug_encounters:
            delta_t *= 0.1

        initial_delta_t = delta_t
        if self.global_debug > 1:
            print(pre, 'evolving to time', end_time)
            print(pre, 'initial step =', initial_delta_t)

        # Allow delta_t to increase, with an upper limit.  (The factor
        # of 25 below should permit quasi-stable systems to be
        # detected.)

        delta_t_max = 64*delta_t
        while delta_t_max < end_time/25: delta_t_max *= 2

        # Save some useful initial quantities.

        initial_position = particles.position
        initial_velocity = particles.velocity
        initial_cmvel = particles.center_of_mass_velocity()   # should be 0
        initial_ke = particles.kinetic_energy()
        initial_end_time = end_time

        # Allow the possibility of repeating the encounter if it fails
        # to terminate.

        loop_count = 0
        loop_max = 10
        pert = 0.001					  # retry option 1
        pert_fac = 10.**(1./loop_max)
        scale_fac = (min_scatter_scale
                     / final_scatter_scale)**(2./loop_max)	# option 2
        end_time_fac = 1.5					# option 2
        over = 0

        while loop_count < loop_max:

            loop_count += 1
            #print pre, 'loop_count =', loop_count
            #sys.stdout.flush()

            resolve_collision_code \
              = self.resolve_collision_code_creation_function()

            # Channel to copy values from the code to the set in memory.

            channel = resolve_collision_code.particles.new_channel_to(particles)

            time = 0 * end_time

            resolve_collision_code.set_time(time)
            resolve_collision_code.particles.add_particles(particles)
            resolve_collision_code.commit_particles()

            delta_t = initial_delta_t
            resolve_collision_code.set_break_scale(final_scatter_scale)

            initial_scatter_energy \
              = self.get_total_energy(resolve_collision_code)

            if self.global_debug > 1:
                print(pre, 'number_of_stars =', len(particles), ' ', \
                    particles.id)
                print(pre, 'initial energy =', initial_scatter_energy)
            #print particles

            if self.debug_encounters:
                print(pre, '### START ENCOUNTER ###')
                print(pre, '### snapshot at time %f' % 0.0)
                for p in particles:
                    print(pre, '### id=%d, x=%f, y=%f, z=%f,'\
                        'vx=%f, vy=%f, vz=%f' % \
                        (p.id, p.x.number, p.y.number, p.z.number,
                         p.vx.number, p.vy.number, p.vz.number))
                sys.stdout.flush()

            #------------------------------------------------------------
            #
            # If the encounter fails to terminate within the specified
            # time we have some options:
            #
            # 1. Try perturbing the encounter in various energy
            # conservative ways, starting from the original
            # velocities.
            #
            # 2. Modify the termination conditions.  This is
            # potentially less expensive, but may not lead to a clean
            # outcome.  Increasing end_time simply involves extending
            # the while loop; changing final_scatter_scale requires a
            # new calculation.

            option = 2
            inner_loop = 1

            #############################################################
            # Set this to enable step-by-step debugging output.
            # resolve_collision_code.parameters.outfile='abc.dat'
            #
            # e.g.
            # if self.gravity_code.model_time.number > 31.4159:
            #     resolve_collision_code.parameters.outfile = 'debug.dat'
            #############################################################

            while time < end_time:

                tt = time
                time += delta_t
                # print pre, '...to time', time
                # sys.stdout.flush()

                # Work with internal substeps of initial_delta_t to
                # allow checks for quasi-stable motion.

                while tt < time:

                    tt += initial_delta_t
                    if tt > time: tt = time

                    if 0:
                        print(pre, '    ...', time, tt, \
                            'model_time =', \
                            resolve_collision_code.model_time)
                        sys.stdout.flush()

                    resolve_collision_code.evolve_model(tt)

                    if 0:
                        print(pre, '    ...back:', \
                            ': model_time =', \
                            resolve_collision_code.model_time)
                        sys.stdout.flush()

                    tt = resolve_collision_code.model_time

                    # DEBUGGING:
                    if self.debug_encounters:
                        print(pre, '### snapshot at time %f' \
                            		% time.number)
                        #resolve_collision_code.update_particle_tree()
                        #resolve_collision_code.update_particle_set()
                        resolve_collision_code.particles \
					.synchronize_to(particles)
                        channel.copy()
                        for p in particles:
                            print(pre, '### id=%d, x=%f, y=%f, z=%f,'\
                              'vx=%f, vy=%f, vz=%f' % \
                              (p.id, p.x.number, p.y.number, p.z.number,
                                 p.vx.number, p.vy.number, p.vz.number))
                        sys.stdout.flush()

                    # The argument final_scatter_scale is used to
                    # limit the size of the system.  It has to be
                    # supplied again because the code that determines
                    # if the scattering is over isn't necessarily the
                    # same as resolve_collision_code.  However,
                    # currently only smallN has an "is_over()"
                    # function.
                    #
                    # Return values:	0 - not over
                    #			1 - over
                    #			2 - quasi-stable system
                    #			3 - size exceeded limit
                    #
                    # Note that this is really a stopping condition,
                    # and should eventually be handled that way.  TODO
                    #
                    # If over = 3, if the parameters were properly
                    # chosen, the resulting system should stil be
                    # usable.  The interface function will take steps
                    # to return proper hierarchical structure even if
                    # the inner subsystem is not well resolved.
                    #
                    # Note that we are currently ignoring any
                    # possibility of a physical collision during the
                    # multiples encounter.  TODO

                    over = resolve_collision_code.is_over\
                           		(final_scatter_scale,
                                         0)    # verbose = 0
                                                          
                    if over:
                        final_scatter_energy \
                          = self.get_total_energy(resolve_collision_code)
                        scatter_energy_error \
                          = final_scatter_energy - initial_scatter_energy

                        if self.global_debug > 1:
                            print(pre, 'over =', over, 'at time', tt)
                            #print pre, 'initial energy =', \
                            #      initial_scatter_energy
                            #print pre, 'final energy =', \
                            #      final_scatter_energy
                            #print pre, 'energy error =', \
                            #      scatter_energy_error
                            print(pre, 'fractional energy error =', \
                                scatter_energy_error/initial_scatter_energy)
                        if self.debug_encounters:
                            print(pre, '### END ENCOUNTER ###')
                            sys.stdout.flush()

                        # Create a tree in the module representing the
                        # binary structure.

                        resolve_collision_code.update_particle_tree(over)

                        # Note: A quasi-stable system (over = 2)
                        # should be handled properly, as it will
                        # appear to be a bound top-level binary.  If
                        # over = 3, the top level should be a receding
                        # bound or unbound system, and the tree
                        # structure should still be usable.

                        # Note that center of mass particles are now
                        # part of the particle set.

                        # Return the tree structure to AMUSE.
                        # Children are identified by
                        # get_children_of_particle in interface.??,
                        # and the information is returned in the copy
                        # operation.

                        resolve_collision_code.update_particle_set()
                        resolve_collision_code.particles \
                                            .synchronize_to(particles)
                        #print "resolve_collision_code.particles.radius"\
                        #      , resolve_collision_code.particles.radius
                        channel.copy()
                        #resolve_collision_code.stop()

                        if 1:

                            # Count number of top-level multiples.  Must
                            # be >0 for the post-encounter bookkeeping to
                            # work.

                            binaries = trees.BinaryTreesOnAParticleSet(
                            			particles, "child1", "child2")
                            singles = binaries.particles_not_in_a_multiple()
                            multiples = binaries.roots()
                            if self.global_debug > 0:
                                print('after', pre, len(singles), \
                                      'single(s),', \
                                      len(multiples), 'multiple(s)')

                        return scatter_energy_error

                    if tt >= 0.99999999*time: break	# avoid roundoff

                # -- end of while tt < time: loop --

                time = resolve_collision_code.model_time
                if not self.debug_encounters:
                    if delta_t < delta_t_max \
                      and time > 0.999999*4*delta_t:
                        delta_t *= 2
                        if self.global_debug > 1:
                            print(pre, 'setting delta_t =', delta_t)
                        sys.stdout.flush()

                if time > 0.99999999*end_time:		# avoid roundoff

                    # Encounter has failed to terminate and we are
                    # about to break out of the loop.  If option = 2
                    # and this is an odd-numbered loop, just increase
                    # end_time (once only).  Otherwise, break and
                    # allow other options to take effect.

                    if option == 2 and 2*(loop_count/2) != loop_count \
                      and inner_loop == 1: 

                        # Adjust the bulk scattering parameters.
                        # Simply increase end_time.

                        end_time *= end_time_fac

                        # Same print output as below.

                        if self.global_debug > 1:
                            print(pre, 'loop', loop_count, ' over =', over)
                            print('increasing end_time to', end_time)
                            print('-----')

                        inner_loop = 2
                        loop_count += 1

                    else:
                        break
                        
            # -- end of while time < end_time: loop --
 
            # As currently coded, if we get here we have not
            # overwritten the original particle set, particles.
            # Nevertheless, we restore particle data here prior
            # to a retry.

            particles.position = initial_position
            particles.velocity = initial_velocity

            if self.global_debug > 1:
                print(pre, 'loop', loop_count, ' over =', over)

            if option == 1:

                # Perturbing the encounter can be done in several
                # ways, of increasing intrusiveness and decreasing
                # reasonableness.
                #
                #     1. Randomize the phases of all binary orbits.
                #     2. Randomize the orientations of all binary orbits.
                #     3. Perturb (jiggle) the top-level orbits.
                #     4. Jiggle all velocities.
                #
                # In all cases, the total energy must be preserved and
                # the CM motion must remain at the origin.  However,
                # in case 4, the total multiple energy and hence the
                # bookkeeping will be compromised unless we explicitly
                # correct it -- need an additional return value.
                #
                # TODO: We should implement options 1-3 -- these
                #       require scattering_stars to be passed as an
                #       argument.
                #
                # For now, choose the least desirable but easiest
                # option #4, with increasing pert as loop_count
                # increases.

                ran = 1 + pert*(2*numpy.random.random(len(particles)) - 1)
                for k in range(len(particles)):
                    particles[k].velocity *= ran[k]

                # Preserve momentum and energy.

                final_cmvel = particles.center_of_mass_velocity()
                particles.velocity -= final_cmvel - initial_cmvel
                final_ke = particles.kinetic_energy()
                particles.velocity *= math.sqrt(initial_ke/final_ke)

                pert *= pert_fac
                print('retrying with pert =', pert)

            elif option == 2:

                # Adjust the bulk scattering parameters.  First
                # increase end_time, then reduce and increase
                # final_scatter_scale, etc.  Should only be able to
                # get here if loop_count is even.  End_time has
                # already been increased.  Use the larger version, but
                # decrease final_scatter_scale.

                final_scatter_scale *= scale_fac

                print('retrying with final_scatter_scale =', final_scatter_scale)
                print('              end_time =', end_time)

            print('-----')

        raise DidNotFinishException(
            pre + \
            " Small-N simulation did not finish before end time {0}".
            format(end_time)
        )

def openup_tree(star, tree, particles_in_encounter):

    # List the leaves.

    leaves = tree.get_leafs_subset()
    #print 'leaves:'
    #print leaves
      
    original_star = tree.particle

    # Compare with the position stored when replacing the particles
    # with the root particle, and move the particles accordingly.
    # Note that once the CM is in the gravity module, the components
    # are frozen and the coordinates are absolute, so we need the
    # original coordinates to offset them later.

    # Maybe better just to store relative coordinates?  TODO

    dx = star.x - original_star.x
    dy = star.y - original_star.y
    dz = star.z - original_star.z
    dvx = star.vx - original_star.vx
    dvy = star.vy - original_star.vy
    dvz = star.vz - original_star.vz
    
    leaves.x += dx
    leaves.y += dy
    leaves.z += dz
    leaves.vx += dvx
    leaves.vy += dvy
    leaves.vz += dvz
    
    particles_in_encounter.add_particles(leaves)

def phi_tidal(star1, star2, star3, G): # compute tidal potential of
                                       # (star1,star2) relative to star3
    phi13 = -G*star1.mass*star3.mass/(star1.position-star3.position).length()
    phi23 = -G*star2.mass*star3.mass/(star2.position-star3.position).length()
    m12 = star1.mass + star2.mass
    cm = Particles([star1, star2]).center_of_mass()
    phicm = -G*m12*star3.mass/(star3.position-cm.position).length
    return phi13+phi23-phicm

def find_nn(plist, field, G):

    # Find and print info on the closest field particle (as
    # measured by potential) to any particle in plist.

    pminmin = None
    fminmin = None
    phiminmin = zero
    for f in field:
        dx = (plist.position - f.position).lengths()
        phi = -G*f.mass*plist.mass/dx
        phimin = zero
        dxmin = 1.e30
        pmin = None
        for i in range(len(phi)):
            if phi[i] < phimin:
                phimin = phi[i]
                dxmin = dx[i]
                pmin = plist[i]
        if phimin < phiminmin:
            phiminmin = phimin
            pminmin = pmin
            dxminmin = dxmin
            fminmin = f

    return pminmin, fminmin, dxminmin

def find_nn2(plist, field, G):
    
    # Find and print info on the closest field particle (as
    # measured by potential) to any particle in plist.
    # revised, faster version of find_nn
    
    pminmin = None
    fminmin = None
    phiminmin = zero
    
    for p in plist:
        
        dx = (p.position - field.position).lengths()
        phi = -G*field.mass*p.mass/dx
        #phi = numpy.divide(numpy.prod([-1,G,field.mass,p.mass]),dx)
        phimin = zero
        dxmin = 1.e30
        pmin = None
        j = numpy.argmin(phi.number)
        phimin = phi[j]
        dxmin = dx[j]
        pmin = p
        if phimin < phiminmin:
            phiminmin = phimin
            pminmin = pmin
            dxminmin = dxmin
            fminmin = field[j]

    return pminmin, fminmin, dxminmin

def find_binaries(particles, G):

    # Search for and print out bound pairs using a numpy-accelerated
    # N^2 search.

    for p in particles:
        mu = p.mass*particles.mass/(p.mass+particles.mass)
        dr = (particles.position - p.position).lengths()
        dv = (particles.velocity - p.velocity).lengths()
        E = 0.5*mu*dv*dv - G*p.mass*particles.mass/dr
        indices = numpy.argsort(E.number)
        sorted_E = E[indices]
        Emin = sorted_E[1].number
        if Emin < -1.e-4 and p.id < particles[indices[1]].id:
            print('bound', p.id, particles[indices[1]].id, Emin)

def potential_energy_in_field(particles, field_particles,
                              smoothing_length_squared = zero,
                              G=constants.G):
    """
    Returns the total potential energy of the particles in the particles
    set. argument field_particles: the external field consists of these
    (i.e. potential energy is calculated relative to the field
    particles) argument smoothing_length_squared: the smoothing length
    is added to every distance.  argument G: gravitational constant,
    need to be changed for particles in different units systems
    """

    if len(field_particles) == 0:
        return zero

    sum_of_energies = zero
    for particle in particles:
        dr_squared = (particle.position-field_particles.position).lengths_squared()
        dr = (dr_squared+smoothing_length_squared).sqrt()
        m_m = particle.mass * field_particles.mass
        potentials = -m_m/dr
        energy_of_this_particle = potentials.sum()
        sum_of_energies += energy_of_this_particle
        imin = numpy.argmin(potentials.number)
        imin = numpy.argmin(dr.number)

    return G * sum_of_energies

def offset_particle_tree(particle, dpos, dvel):

    # Recursively offset a particle and all of its descendants by
    # the specified position and velocity.

    if not particle.child1 is None:
        offset_particle_tree(particle.child1, dpos, dvel)
    if not particle.child2 is None:
        offset_particle_tree(particle.child2, dpos, dvel)
    particle.position += dpos
    particle.velocity += dvel
    # print 'offset', int(particle.id), 'by', dpos; sys.stdout.flush()

def rescale_binary_components(comp1, comp2, kep, scale, compress=True):

    # Rescale the two-body system consisting of comp1 and comp2 to lie
    # inside (compress=True) or outside (compress=False) distance
    # scale of one another.  If compress=True, the final orbit will be
    # receding; otherwise it will be approaching.  In a typical case,
    # scale is comparable to the separation at which the interaction
    # started.  It is possible that the input system is very close to
    # periastron.  To avoid problems with very eccentric systems,
    # force the system to be scaled to a separation of at least
    # 0.1*scale (0.1 is ~arbitrary: should be <1, and not too small).

    if compress:
        pre = 'rescale_binary_components(-):'
    else:
        pre = 'rescale_binary_components(+):'

    pos1 = comp1.position
    pos2 = comp2.position
    sep12 = ((pos2-pos1)**2).sum()
    mass1 = comp1.mass
    mass2 = comp2.mass
    total_mass = mass1 + mass2
    vel1 = comp1.velocity
    vel2 = comp2.velocity
    cmpos = (mass1*pos1+mass2*pos2)/total_mass
    cmvel = (mass1*vel1+mass2*vel2)/total_mass

    mass = comp1.mass + comp2.mass
    rel_pos = pos2 - pos1
    rel_vel = vel2 - vel1
    
    if 0:
        print(pre, 'mass =', mass)
        print(pre, 'pos =', rel_pos)
        print(pre, 'vel =', rel_vel)

    kep.initialize_from_dyn(mass,
                            rel_pos[0], rel_pos[1], rel_pos[2],
                            rel_vel[0], rel_vel[1], rel_vel[2])
    M,th = kep.get_angles()
    a,e = kep.get_elements()

    rescale = (compress and sep12 > scale**2) \
                or (not compress and sep12 < scale**2)

    min_scale = 0.1*scale			# see note above

    if 0:
        print(pre, 'M, th, a, e, =', M, th, a, e)
        print(pre, 'compress =', compress)
        print(pre, sep12, scale**2, min_scale**2)

    if compress == True:
        rescale = rescale or sep12 < min_scale**2

    if rescale:

        #print 'rescaling components', int(comp1.id), \
        #      'and', int(comp2.id), 'to separation', scale
        # sys.stdout.flush()
        
        # print pre, 'a, e =', a, e
        if e < 1:
            peri = a*(1-e)
            apo = a*(1+e)
        else:
            peri = a*(e-1)
            apo = peri+a		# OK - used only to reset scale

        if compress:

            # Logic here is to handle special configurations.

            limit = peri + 1.e-4*(apo-peri)		# numbers are
            if limit > 1.1*peri: limit = 1.1*peri	# ~arbitrary
            if limit < min_scale: limit = min_scale
            if scale < limit:
                # print pre, 'changed scale from', scale, 'to', limit
                scale = limit

            if M < 0:
                # print pre, 'advance_to_periastron'
                kep.advance_to_periastron()
                # print pre, 'advance_to_radius', scale
                kep.advance_to_radius(scale)
            else:
                if kep.get_separation() < scale:
                    # print pre, 'advance_to_radius', scale
                    kep.advance_to_radius(scale)
                else:
                    # print pre, 'return_to_radius', scale
                    kep.return_to_radius(scale)

            # Note: Always end up on an outgoing orbit.  If periastron
            # > scale, we are now just past periastron.

        else:
            limit = apo - 0.01*(apo-peri)
            # print pre, "limit:", limit, apo, peri, scale , M, e
            if scale > limit:
                # print pre, 'changed scale from', scale, 'to', limit
                scale = limit
            
            #print "INPUT:", kep.get_separation_vector()
            #print "true_anomaly:", M, kep.get_separation() , scale
            if M > 0:
                kep.return_to_periastron()
                kep.return_to_radius(scale)
            else:
                if kep.get_separation() < scale:
                    kep.return_to_radius(scale)
                else:
                    kep.advance_to_radius(scale)

	    # Note: Always end up on an incoming orbit.  If
	    # apastron < scale, we are now just before apastron.

        #print 'scale =', scale, 'sep =', kep.get_separation(), 'M =', M

        new_rel_pos = kep.get_separation_vector()
        new_rel_vel = kep.get_velocity_vector()

        # Problem: the vectors returned by kepler are lists, not numpy
        # arrays, and it looks as though we can say comp1.position =
        # pos, but not comp1.position[k] = xxx, as we'd like...  Also,
        # Steve doesn't know how to copy a numpy array with units...
        # TODO - help?
        #print "REL POS:", new_rel_pos
        newpos1 = pos1 - pos1	# stupid trick to create zero vectors
        newpos2 = pos2 - pos2	# with the proper form and units...
        newvel1 = vel1 - vel1
        newvel2 = vel2 - vel2

        frac2 = mass2/total_mass
        for k in range(3):
            dxk = new_rel_pos[k]
            dvk = new_rel_vel[k]
            newpos1[k] = cmpos[k] - frac2*dxk
            newpos2[k] = cmpos[k] + (1-frac2)*dxk
            newvel1[k] = cmvel[k] - frac2*dvk
            newvel2[k] = cmvel[k] + (1-frac2)*dvk

        # Perform the changes to comp1 and comp2, and recursively
        # transmit them to the (currently absolute) coordinates of all
        # lower components.
        #print "DP1:", newpos1-pos1, hasattr(comp1, 'child1')
        #print "DP2:", newpos2-pos2, hasattr(comp2, 'child1')
        if hasattr(comp1, 'child1'):
            offset_particle_tree(comp1, newpos1-pos1, newvel1-vel1)
        if hasattr(comp2, 'child1'):
            offset_particle_tree(comp2, newpos2-pos2, newvel2-vel2)

    # print pre, 'done'
    sys.stdout.flush()

    return a

def offset_children(n, dx, dv):
    if n.child1 != None:
        n.child1.position -= dx
        n.child1.velocity -= dv
        offset_children(n.child1, dx, dv)
    if n.child2 != None:
        n.child2.position -= dx
        n.child2.velocity -= dv
        offset_children(n.child2, dx, dv)

def compress_nodes(node_list, scale, G):

    local_debug = False

    # Compress (or expand) the top-level nodes in node_list to lie
    # within diameter scale.  Rescale velocities to conserve total
    # energy (but currently not angular momentum -- TODO).

    pre = 'compress_nodes:'

    # Compute the center of mass position and velocity of the
    # top-level system.

    cmpos = node_list.center_of_mass()
    cmvel = node_list.center_of_mass_velocity()

    # Child positions and velocities will not be explicitly changed by
    # the scaling.  Temporarily store child data as offsets relative
    # to the root.  We will undo this at the end, immediately before
    # returning.

    for n in node_list:
        if n.child1 != None:
            dx = n.position
            dv = n.velocity
            offset_children(n, dx, dv)

    if local_debug:
        print('node_list:')
        print(node_list)
        print('top_level:')
        print_top_level(node_list, G)

    x0 = (node_list[0].position**2).sum().sqrt()
    lunit = x0/x0.number
    v0 = (node_list[0].velocity**2).sum().sqrt()
    vunit = v0/v0.number
    vunit2 = vunit**2
    
    # Compute various measures of the size, potential, and kinetic
    # energy of the system in the center of mass frame.

    size = zero			# max distance(**2) from center of mass

    rijmin = 1.e100*lunit	# minimum separation
    imin = -1
    jmin = -1
    phimin = zero		# minimum potential
    ipmin = -1
    jpmin = -1

    n = len(node_list)
    pot = zero
    kin = zero
    dr = numpy.zeros((n,n))	# unit = lunit
    dv2 = numpy.zeros((n,n))	# unit = vunit2
    for i in range(n):
        m = node_list[i].mass
        posi = node_list[i].position
        pos = posi - cmpos
        veli = node_list[i].velocity
        vel = veli - cmvel
        r2 = (pos**2).sum()
        if r2 > size:
            size = r2
        kin += m*(vel**2).sum()
        dpot = zero
        for j in range(i+1,n):
            mj = node_list[j].mass
            dposj = node_list[j].position - posi
            rij = (dposj**2).sum().sqrt()
            dphij = -G*mj/rij
            dpot += dphij
            phij = m*dphij
            if rij < rijmin:
                rijmin = rij
                imin = i
                jmin = j
            if phij < phimin:
                phimin = phij
                ipmin = i
                jpmin = j
            dvelj = node_list[j].velocity - veli
            dr[i,j] = rij/lunit
            dv2[i,j] = (dvelj**2).sum()/vunit2
        if dpot != zero:
            pot += m*dpot
    size = size.sqrt()
    kin /= 2
    rphmin = -(node_list[ipmin].mass*node_list[jpmin].mass)/phimin

    if local_debug:
        print(pre, 'scale =', scale)
        print(pre, 'size =', size)
        print(pre, 'rijmin =', rijmin, node_list[imin].id, node_list[jmin].id)
        print(pre, 'rphmin =', rphmin, node_list[ipmin].id, node_list[jpmin].id)

    fac = 0.5*scale/size		# scale to radius
    #fac = scale/rijmin			# scale to minimum distance
    #fac = scale/rphmin			# scale to minimum potential distance

    if local_debug:
        print(pre, 'fac =', fac)

    # Compress (or expand) the system and increase (or decrease) the
    # velocities (relative to the center of mass) to preserve the
    # energy.  If fac > 1, expansion is always OK if E > 0, which it
    # should be at this point (but check anyway...).  May have E < 0
    # if we have a system with small negative energy, stopped because
    # it is too big.

    # An additional consideration (Steve, 1/2017) is that all
    # top-level nodes are mutually unbound at the end of the
    # scattering, by construction, but this may not be preserved by
    # a simple uniform rescaling of the system.  In that case, an
    # unphysical extra interaction may follow the scattering we
    # thought was "over."  Currently we check for this possibility,
    # then modify the way in which velocities are scaled to
    # compensate.  NOT guaranteed to work in all cases, and the code
    # is ugly...

    vfac2 = 1-(1/fac-1)*pot/kin
    #print "vfac2 =", vfac2

    if vfac2 < 0:
        print(pre, "Can't expand top level system to rjmin > ri+rj")
        print("fac =", fac, " pot =", pot, " kin =", kin)
        sys.stdout.flush()
        f = pot/(kin+pot)
        vfac2 = 0.0		# ???

    vfac = math.sqrt(vfac2)
    if local_debug:
        print("vfac =", vfac)
        print(pre, 'dr:')
        print(dr)
        print(pre, 'dv2:')
        print(dv2)
    
    bound_pairs = []
    unbound = numpy.ones(n)
    for i in range(n):
        mi = node_list[i].mass
        bound = False
        for j in range(i+1,n):
            mj = node_list[j].mass
            mu = mi*mj/(mi+mj)
            Eijold = 0.5*mu*dv2[i,j]*vunit2 - G*mi*mj/(dr[i,j]*lunit)
            Eijnew = 0.5*mu*vfac2*dv2[i,j]*vunit2 - G*mi*mj/(fac*dr[i,j]*lunit)
            if Eijnew.number <= 0.0:
                #print 'bound', i, j, Eijold, Eijnew
                bound = True
                bound_pairs.append((i,j))
                unbound[i] = 0
                unbound[j] = 0

    print(pre, 'bound pairs:', bound_pairs)
    unbound_nodes = []
    for i in range(n):
        if unbound[i] == 1:
            unbound_nodes.append(i)
    print(pre, 'unbound_nodes:', unbound_nodes)

    if len(unbound_nodes) == 0:

        # Live with unphysical bound pairs for now.  TODO

        print('*** warning: no unbound nodes ***')
        bound_pairs = []

    if len(bound_pairs) > 0:

        # Strategy #1: Scale positions uniformly as planned, but
        # adjust the velocity scaling for pairs whose binding energy
        # would become negative.  Strategy #2 (unimplemented) would be
        # to modify the spatial scaling by scaling bound pairs only to
        # scale, then adjust velocities as in #1.  Strategy #3 (also
        # unimplemented) would be to not scale bound pairs at all.

        energy = pot + kin		# initial energy - conserved

        for n in node_list:
            n.position = cmpos + fac*(n.position-cmpos)
        dr *= fac
        pot /= fac

        if local_debug:
            print('kinetic energies:')
            for n in node_list:
                print('  ', n.id, 0.5*n.mass*((n.velocity-cmvel)**2).sum())

        # First give the bound components enough relative velocity to
        # just unbind them, keeping their center of mass velocity
        # fixed.  Note that, since Eij was > 0 and this prescription
        # leaves Eij close to 0, this transformation should liberate
        # energy for distribution to the rest of the system.

        kin2 = zero
        kinCM = zero
        for p in bound_pairs:
            i = p[0]
            j = p[1]
            ni = node_list[i]
            nj = node_list[j]
            mi = ni.mass
            mj = nj.mass
            newvfac2 = 2.000001*(G*(mi+mj)/(dr[i,j]*lunit))/(dv2[i,j]*vunit2)
            newvfac = math.sqrt(newvfac2)
            massinv = 1./(mi+mj)
            cmv = (mi*ni.velocity + mj*nj.velocity)*massinv
            ni.velocity = cmv + newvfac*(ni.velocity-cmv)
            nj.velocity = cmv + newvfac*(nj.velocity-cmv)
            kin2 += 0.5*mi*mj*massinv*((ni.velocity-nj.velocity)**2).sum()
            kinCM += 0.5*(mi+mj)*((cmv-cmvel)**2).sum()

        if local_debug:
            print('KECM =', kin2+kinCM)
        for i in unbound_nodes:
            ni = node_list[i]
            mi = ni.mass
            kei = 0.5*mi*((ni.velocity-cmvel)**2).sum()
            if local_debug:
                print('KE', ni.id, kei)
            kinCM += kei

        if local_debug:
            print('energy =', energy, 'pot+kin2+kinCM =', pot+kin2+kinCM)
        kin_to_distribute = energy - (pot+kin2+kinCM)

        if kin_to_distribute.number < 0: 
            print('*** warning: not enough kinetic energy ***')	# TODO

        vfac2 = 1+kin_to_distribute/kinCM
        vfac = math.sqrt(vfac2)
        #print 'vfac =', vfac

        # Then apply an overall scaling to unbound nodes and bound CMs
        # to conserve total energy.

        for i in unbound_nodes:
            ni = node_list[i]
            ni.velocity = cmvel + vfac*(ni.velocity-cmvel)
        for p in bound_pairs:
            i = p[0]
            j = p[1]
            ni = node_list[i]
            nj = node_list[j]
            mi = ni.mass
            mj = nj.mass
            massinv = 1./(mi+mj)
            cmv = (mi*ni.velocity + mj*nj.velocity)*massinv
            newcmv = cmvel + vfac*(cmv-cmvel)
            ni.velocity += newcmv - cmv
            nj.velocity += newcmv - cmv

    if len(bound_pairs) == 0:

        # Perform global scaling of position and velocity.

        for n in node_list:
            n.position = cmpos + fac*(n.position-cmpos)
            n.velocity = cmvel + vfac*(n.velocity-cmvel)

    # Child data have not yet been modified.  Do so here.  Note that
    # child positions and velocities were temporarily offset to the
    # top-level center of mass.

    for n in node_list:
        if n.child1 != None:
            dx = -n.position
            dv = -n.velocity
            offset_children(n, dx, dv)

    #print_top_level(node_list, G)

def get_multiple_energy(node, kep):

    # Return the binary status and the total energy Etot of the
    # specified tree node.  The value of Etot is the total pairwise
    # energy of all binary objects in the hierarchy.  It does not
    # include the tidal potential of one component on another (e.g.
    # in a hierarchical triple Etot will be the sum of two binary
    # energies only).

    # Note that kep should have been initialized with the correct
    # converter to return the proper energy units.

    is_bin = 1
    Etot = zero
    for level, x in node.iter_levels():
        particle = x
        if not particle.child1 is None:
            if level > 0: is_bin = 0
            child1 = particle.child1
            child2 = particle.child2
            M,a,e,r,Emu,t = get_component_binary_elements(child1, child2, kep)
            mu = child1.mass*child2.mass/M
            E = Emu*mu
            Etot += E
    return is_bin, Etot

def get_multiple_energy2(node, G):

    # Return the binary status and the total energy of the specified
    # tree node.  Uses a value of G supplied by the caller.  The
    # returned value is the total energy of all leaves in the
    # hierarchy, properly including tidal potentials, but excluding
    # the center of mass energy.

    is_bin = 1
    Ecm = zero

    for level, x in node.iter_levels():
        if level == 0:
            particle = x
            M_comp = 0*particle.mass
            vcm_comp = M_comp*particle.velocity

            if particle.id == 1000074:
                pp = True

            break

    # List the leaves and do some additional work.  Note that
    # node.get_leafs_subset() seems to do the same thing...

    leaves_in_node = Particles(0)

    for level, x in node.iter_levels():
        particle = x
        if level == 0:

            # Want to compute the top-level kinetic energy.  Might
            # expect

            vcm = particle.velocity
            Ecm = 0.5*particle.mass*(vcm**2).sum()

            # but in some circumstances (e.g. a binary created in a
            # many-body process), the position and velocity of the
            # parent may not correctly reflect the center of mass of
            # its children.
            
        if not particle.child1 is None:
            if level > 0:
                is_bin = 0				# not a multiple
        else:
            leaves_in_node.add_particle(particle)	# list leaves
            M_comp += particle.mass
            vcm_comp += particle.mass*particle.velocity

    vcm_comp /= M_comp
    Ecm_comp = 0.5*M_comp*(vcm_comp**2).sum()

    return is_bin, leaves_in_node.kinetic_energy() \
 		           + leaves_in_node.potential_energy(G=G) \
		           - Ecm_comp

def add_leaves(node, leaf_list):
    if node.child1 == None:
        leaf_list.add_particle(node)
    else:
        add_leaves(node.child1, leaf_list)
        add_leaves(node.child2, leaf_list)

def get_energy_of_leaves(particles, G):
    leaves = Particles(0)
    for p in particles:
        if not hasattr(p, 'child1') or p.child1 == None:
            leaves.add_particle(p)
    #print 'get_energy_of_leaves():'
    #print leaves
    ke = leaves.kinetic_energy()
    pe = leaves.potential_energy(G=G)
    e = ke+ pe
    #print ke, pe, e
    return e

def print_energies(stars, G):

    # Brute force N^2 over top level, pure python...

    top_level = stars.select(is_not_a_child, ["is_a_child"])

    mass = zero
    kinetic = zero
    potential = zero
    for t in top_level:
        m = t.mass
        x = t.x
        y = t.y
        z = t.z
        vx = t.vx
        vy = t.vy
        vz = t.vz
        mass += m
        kinetic += 0.5*m*(vx**2+vy**2+vz**2)
        dpot = zero
        for tt in top_level:
            if tt != t:
                mm = tt.mass
                xx = tt.x-x
                yy = tt.y-y
                zz = tt.z-z
                dpot -= G*mm/(xx**2+yy**2+zz**2).sqrt()
        potential += 0.5*m*dpot
            
    print('len(stars) =', len(stars))
    print('len(top_level) =', len(top_level))
    print('mass =', mass)
    print('kinetic =', kinetic)
    print('potential =', potential)
    print('energy =', kinetic+potential)
    sys.stdout.flush()

def scale_top_level_list(singles, multiples, kep, scale,
                         gravity_constant, global_debug):

    pre = 'scale_top_level_list:'

    # The multiples code followed the particles until their
    # interaction could be unambiguously classified as over.  They may
    # now be very far apart.  Input singles and multiples are lists
    # describing the final top-level structure of the interacting
    # particles in the multiples code.  Singles is a list of single
    # stars.  Multiples is a list of multiple centers of mass (with
    # pointers to the internal structure).

    # Scale the positions and velocities of the top-level nodes to
    # bring them within a sphere of diameter scale, conserving energy
    # and angular momentum (if possible).  Also offset all children to
    # reflect changes at the top level -- TODO: will change if/when
    # offsets are implemented...

    # Logic: 1 node   - must be a binary, use kepler to reduce to scale
    #        2 nodes  - use kepler, reduce binary children too?  TODO
    #        3+ nodes - shrink radii and rescale velocities to preserve
    #                   energy, but this doesn't preserve angular
    #                   momentum  TODO - also reduce children?  TODO

    top_level_nodes = singles + multiples

    # Figure out the tree structure.

    ls = len(singles)
    lm = len(multiples)
    lt = ls + lm

    if global_debug > 1:
        print(pre, 'ls =', ls, ' lm =', lm, ' lt =', lt)
        sys.stdout.flush()

    if lt == 1:
        if lm == 1:

            # Special case.  We have a single bound binary node.  Its
            # children are the components we want to transform.  Note
            # that, if the components are binaries (or multiples),
            # they must be stable, so it is always OK to move the
            # components to periastron.

            # Note: Wide binaries will be split and returned to the
            # large-scale dynamics module after return from this
            # function.

            root = multiples[0]

            if global_debug > 1:
                print(pre, 'bound binary node', scale)
            #print '\nunscaled binary node:'
            #print_multiple_recursive(root)
            comp1 = root.child1
            comp2 = root.child2
            if global_debug > 1:
                print(pre, 'scale:', scale)
            semi = rescale_binary_components(comp1, comp2, kep, scale)
            #true, mean = kep.get_angles()
            #print 'true =', true, 'mean =', mean
            #print 'scaled binary node:'
            #print_multiple_recursive(root, kep)

    elif lt == 2:

        # We have two unbound top-level nodes, and we will scale them
        # using kepler to the desired separation, hence conserving
        # both energy and angular momentum of the top-level motion.

        # We might also want to scale the daughter nodes.  Note as
        # above that, if the daughters are binaries (or multiples),
        # they must be stable, so it is always OK to move them to
        # periastron.

        comp1 = top_level_nodes[0]
        comp2 = top_level_nodes[1]

        if global_debug > 1:
            print(pre, 'top-level unbound pair')
            #print pre, '\nunscaled top-level pair:'
            #print_pair_of_stars('pair', comp1, comp2, kep)
            sys.stdout.flush()
        semi = rescale_binary_components(comp1, comp2, kep, scale)
        #print pre, '\nscaled top-level pair:'
        #print_pair_of_stars('pair', comp1, comp2, kep)
        #sys.stdout.flush()

    else:

        # We have three or more unbound top-level nodes.  We don't
        # know how to compress them in a conservative way.  For now,
        # we will conserve energy and think later about how to
        # preserve angular momentum.  TODO

        print(pre, lt, 'top-level nodes, scale =', scale)
        #print lt, 'unscaled top-level nodes'
        #print top_level_nodes
        compress_nodes(top_level_nodes, scale, gravity_constant)
        #print lt, 'scaled top-level nodes'
        #print top_level_nodes

    # print pre, 'done'
    sys.stdout.flush()

    # Don't attempt to correct or even return the tidal energy error.
    # Manage all of this in the calling function, as desired.

    return

def set_radius_recursive(node, kep, global_debug):

    if node.is_leaf(): return		# nothing to be done

    # Propagate child radii upward.  Since dynamical radius scales
    # with mass, the radius of a parent is just the sum of the radii
    # of the children.  If we are simply handling 2-body encounters,
    # that's all we need.  The semi-major axis of a hard binary is
    # less than the dynamical radius, by definition.  However, we must
    # include the size of a soft binary or multiple, which may be
    # significantly larger than the dynamical radius of the center of
    # mass.
    
    # Note that kep should have been initialized with the correct
    # converter to return the proper energy units.

    rsum = zero
    for child in node.iter_children():
        set_radius_recursive(child, kep, global_debug)
        rsum += child.particle.radius

    # Currently rsum is the dynamical radius of the node. Check how it
    # compares to the node's semimajor axis.

    M,semi,e,x,x,x = get_component_binary_elements(node.particle.child1,
                                                   node.particle.child2, kep)

    if rsum < 2*semi:

        # *** Factor of 2 here is ~arbitrary; should probably be set
        # *** in the class definition.

        if global_debug > 0:
            print('increasing radius for', node.particle.id, 'from', \
                rsum, 'to', 2*semi)
        rsum = 2*semi

    node.particle.radius = rsum

# Note: iter_children() lists the leaves lying under a given node of
# type BinaryTreeOnParticle.  The child leaves are objects of type
# ChildTreeOnParticle.  The particle associated with child x is
# x.particle.

def set_radii(top_level_nodes, kep, global_debug):
    for n in top_level_nodes.as_binary_tree().iter_children():
        set_radius_recursive(n, kep, global_debug)

#----------------------------------------------------------------------

def print_elements(s, a, e, r, Emu, E):

    # Print binary elements in standard form.

    print('{0} elements  a = {1}  e = {2}  r = {3}  E/mu = {4}  E = {5}'\
            .format(s, a, e, r, Emu, E))
    sys.stdout.flush()

def print_pair_of_stars(s, star1, star2, kep):
    m1 = star1.mass
    m2 = star2.mass
    M,a,e,r,E,t = get_component_binary_elements(star1, star2, kep)
    print_elements(s, a, e, r, E, E*m1*m2/(m1+m2))
    print_multiple_recursive(star1, kep)
    print_multiple_recursive(star2, kep)
    
def print_multiple_recursive(m, kep, level=0):	  ##### not working? #####

    # Recursively print the structure of (multiple) node m.

    print('    '*level, 'key =', m.key, ' id =', int(m.id))
    print('    '*level, '  mass =', m.mass)
    print('    '*level, '  pos =', m.position)
    print('    '*level, '  vel =', m.velocity)
    sys.stdout.flush()
    if not m.child1 is None and not m.child2 is None:
        M,a,e,r,E,t = get_component_binary_elements(m.child1, m.child2, kep)
        print_elements('   '+'    '*level+'binary', a, e, r, E,
                       (E*m.child1.mass*m.child2.mass/M))
    if not m.child1 is None:
        print_multiple_recursive(m.child1, kep, level+1)
    if not m.child2 is None:
        print_multiple_recursive(m.child2, kep, level+1)

def print_multiple_simple(node, kep):
    for level, x in node.iter_levels():
        output = ''
        if level == 0: output += 'Multiple '
        output += '    ' * level
        particle = x
        output += "{0} id = {1} mass = {2} radius = {3}".format(particle.key,
                                                   particle.id,
                                                   particle.mass.number,
                                                   particle.radius.number)
        if not particle.child1 is None:
            child1 = particle.child1
            child2 = particle.child2
            M,a,e,r,E,t = get_component_binary_elements(child1, child2, kep)
            mu = child1.mass*child2.mass/M
            output += " semi = {0} energy = {1}".format(a.number,
                                                        (mu*E).number)
        print(output)
        sys.stdout.flush()

def print_multiple_detailed(node, kep, pre, kT, dcen):
    is_bin = 1
    Etot = zero
    for level, x in node.iter_levels():
        particle = x
        init = pre
        if level == 0: init += 'Multiple '
        init += '    ' * level
        id = particle.id
        M = particle.mass.number
        if not particle.child1 is None:
            if level > 0: is_bin = 0
            child1 = particle.child1
            child2 = particle.child2
            idlow = min(child1.id, child2.id)
            idhigh = max(child1.id, child2.id)
            print('%s%d (%d,%d) m=%.5f' % (init, id, idlow, idhigh, M), end=' ')
            sys.stdout.flush()
            M,a,e,r,Emu,t = get_component_binary_elements(child1, child2, kep)
            cm = (child1.mass*child1.position + child2.mass*child2.position)/M
            mu = child1.mass*child2.mass/M
            E = Emu*mu
            Etot += E
            D = 0.
            for k in range(3):
                D += (cm[k].number - dcen[k].number)**2
            D = numpy.sqrt(D)
            print('a=%.6f e=%4f r=%6f D=%.4f E/mu=%.5f E=%.5f E/kT=%.5f' % \
                    (a.number, e, r.number, D, Emu.number, E.number, E/kT))
            sys.stdout.flush()
        else:
            print('%s%d m=%.5f' % (init, id, M))
            sys.stdout.flush()

    return is_bin, Etot

def print_top_level(nodes, G):

    # Print various top-level quantities of interest during rescaling.

    print('')
    print('distances:')
    for i in nodes:
        print(i.id, '    ', end=' ')
        for j in nodes:
            if j.id != i.id:
                rij = (j.position-i.position).length()
                print(j.id, rij, '    ', end=' ')
        print('')

    print('radial velocities:')
    for i in nodes:
        print(i.id, '    ', end=' ')
        for j in nodes:
            if j.id != i.id:
                rij = (j.position-i.position).length()
                vdotr = ((j.velocity-i.velocity)*(j.position-i.position)).sum()
                print(j.id, vdotr/rij, '    ', end=' ')
        print('')

    print('potentials:')
    for i in nodes:
        print(i.id, '    ', end=' ')
        mi = i.mass
        for j in nodes:
            if j.id != i.id:
                mj = j.mass
                rij = (j.position-i.position).length()
                print(j.id, -G*mi*mj/rij, '    ', end=' ')
                
        print('')

    print('energies:')
    pot = 0.0
    kin = 0.0
    for i in nodes:
        print(i.id, '    ', end=' ')
        mi = i.mass
        vi = i.velocity
        kin += 0.5*mi*(vi**2).sum()
        for j in nodes:
            if j.id != i.id:
                mj = j.mass
                muij = mi*mj/(mi+mj)
                rij = (j.position-i.position).length()
                vij2 = (j.velocity-i.velocity).length_squared()
                print(j.id, 0.5*muij*vij2 - mi*mj/rij, '    ', end=' ')
                if j.id > i.id:
                    pot -= G*mi*mj/rij
                
        print('')
    print('totals:', pot, kin, -kin/pot, pot+kin)
    print('')
