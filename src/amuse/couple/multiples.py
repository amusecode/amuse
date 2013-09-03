import sys
import numpy
import collections
import math

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.datamodel import trees
from amuse.datamodel import Particles 
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import zero

#---------------------------------------------------------------------
#
# Steve's ToDo list of features to be added/improved in the multiples
# module (last modified 12 Mar 2013):
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
# large internal/external tidal errors.  As an alternative, we can let
# near neighbors simply veto the encounter, moving the work back into
# the gravity module, until a "clean" 2-body scattering can be
# identified.
#
# 4. If a scattering consists of several long excursions, it should be
# broken up into pieces, and returned to the gravity code during each
# excursion.  In that case, the "smallN" code will terminate with a
# configuration that isn't "over," in the sense used in the current
# code.  It is unclear if the structure returned from
# update_particle_tree() is correct in that case.  Function is_over()
# returns a stopping condition -- an integer representing the reason
# for its return -- but there is no code to accommodate a return in an
# intermediate state.
#
# 5. We should seek a better prescription for compressing 3-body and
# higher-order configurations.  Currently we conserve energy, but not
# angular momentum.
#
# 6. There is no provision for physical collisions in the smallN code,
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

    # Determine the object's description, then search to see
    # if we know about it.  If we do, return that ID, otherwise
    # create a new ID.

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

class Multiples(object):

    def __init__(self, 
                 gravity_code,
                 resolve_collision_code_creation_function,
                 kepler_code, 
                 gravity_constant = None, **options):
        self.gravity_code = gravity_code
        self._inmemory_particles = self.gravity_code.particles.copy()
        self._inmemory_particles.id = self.gravity_code.particles.index_in_code
                
        self.channel_from_code_to_memory = \
            self.gravity_code.particles.new_channel_to(self._inmemory_particles)
        
        self.resolve_collision_code_creation_function \
            = resolve_collision_code_creation_function
        self.kepler = kepler_code
        self.multiples_external_tidal_correction \
            = zero * self.gravity_code.kinetic_energy
        self.multiples_internal_tidal_correction \
            = self.multiples_external_tidal_correction
        self.multiples_integration_energy_error \
            = self.multiples_external_tidal_correction
        self.root_to_tree = {}
        if gravity_constant is None:
            gravity_constant = nbody_system.G
        
        self.multiples = datamodel.Particles()
        self.gravity_constant = gravity_constant
        self.debug_encounters = False

        # The following tunable parameters govern the multiples logic:

        # Size of the top-level encounter, relative to the sum of the
        # radii of the interacting components (no veto) or their
        # separation (veto).

        self.initial_scale_factor = 1.0

        # Distance within which to include a neighbor, relative to the
        # initial separation of the top-level two-body encounter.
        
        self.neighbor_distance_factor = 1.0
        #self.neighbor_distance_factor = 2.0

        # Neighbor veto policy.  True means we allow neighbors to veto
        # a two-body encounter (meaning that don't want to deal with
        # complex initial many-body configurations).  False means we
        # include neighbors in the multiple integration.

        self.neighbor_veto = False
        #self.neighbor_veto = True

        # Size of the rescaled final system, relative to the initial
        # scale.  Should be 1 + epsilon.

        self.final_scale_factor = 1.01

        # Initial separation for the scattering experiment, relative
        # to the initial scale.

        self.initial_scatter_factor = 10.0

        # Final separation for the scattering experiment, relative to
        # the initial scattering scale (meaning 10 x 10 times the
        # initial encounter scale, by default).

        self.final_scatter_factor = 10.0

        # Binary retention policy.  Retain a binary if its apocenter
        # (True) or 2*semi-major axis (False) is less than the
        # dynamical radius of its CM.  False is the more conservative
        # choice.

        self.retain_binary_apocenter = True

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
        for root, tree in self.root_to_tree.iteritems():
            root_particle = root.as_particle_in_set(self._inmemory_particles)
            result.remove_particle(root)
            leaves = tree.get_leafs_subset()
      
            original_star = tree.particle

            dx = root_particle.x - original_star.x
            dy = root_particle.y - original_star.y
            dz = root_particle.z - original_star.z
            dvx = root_particle.vx - original_star.vx
            dvy = root_particle.vy - original_star.vy
            dvz = root_particle.vz - original_star.vz
            
            leaves_in_result = result.add_particles(leaves)
            leaves_in_result.x += dx
            leaves_in_result.y += dy
            leaves_in_result.z += dz
            leaves_in_result.vx += dvx
            leaves_in_result.vy += dvy
            leaves_in_result.vz += dvz
        return result
        
    def get_gravity_at_point(self, radius, x, y, z):
        return self.gravity_code.get_gravity_at_point(radius, x, y, z)
    
    def get_potential_at_point(self, radius, x, y, z):
        return self.gravity_code.get_potential_at_point(radius, x, y, z)
    
    def commit_particles(self):
        return self.gravity_code.commit_particles()
    
    def get_time(self):
        return self.gravity_code.get_time()
    
    def get_total_energy(self, code):
        # ??? from Steve: what is get_binary_energy()?
        try:
            binaries_energy = code.get_binary_energy()	# include binaries
        except:						# if code understands
            binaries_energy = zero
        total_energy = code.potential_energy + code.kinetic_energy \
				             + binaries_energy

        return total_energy

    #--------------------------------------------------------------
    # Note that the true total energy of a multiple isn't quite the
    # Emul returned below, since the tidal potential of components on
    # one another is not taken into account.

    def get_total_multiple_energy(self):
        Nbin = 0
        Nmul = 0
        Emul = 0.0 | nbody_system.energy
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
        Emul = 0.0 | nbody_system.energy
        for x in self.root_to_tree.values():	# loop over top-level trees
            Nmul += 1
            nb,E = get_multiple_energy2(x, self.gravity_constant)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul

    def print_multiples(self):
        for x in self.root_to_tree.values():
            print_simple_multiple(x, self.kepler)
    #--------------------------------------------------------------

    def pretty_print_multiples(self, pre, kT, dcen):
        Nbin = 0
        Nmul = 0
        Emul = 0.0 | nbody_system.energy
        for x in self.root_to_tree.values():
            Nmul += 1
            nb,E = another_print_multiple(x, self.kepler, pre, kT, dcen)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul
            
    def print_trees_summary(self):
        if len(self.root_to_tree) > 0:
            print 'number of multiples:', len(self.root_to_tree)
            sys.stdout.flush()

    def evolve_model(self, end_time):

        stopping_condition = \
            self.gravity_code.stopping_conditions.collision_detection
        stopping_condition.enable()
        
        time = self.gravity_code.model_time
        print "Evolve model to", end_time, " starting at", time
        sys.stdout.flush()

        count_resolve_encounter = 0
        count_ignore_encounter = 0
        
        while time < end_time:

            print '\ncalling evolve_model to ', end_time
            sys.stdout.flush()

#             if time > 100 | nbody_system.time:

#                 # No neighbors, allow vetos, but none will occur.
#                 self.neighbor_distance_factor = 0.
#                 self.neighbor_veto = True

#                 # Veto all multiples treatments.
#                 self.neighbor_distance_factor = 100.
#                 self.neighbor_veto = True

#                 self.neighbor_distance_factor = 0.0
#                 self.neighbor_veto = True

#                 print 'neighbor_distance_factor =', \
#                     self.neighbor_distance_factor
#                 print 'neighbor_veto =', \
#                     self.neighbor_veto

            self.gravity_code.evolve_model(end_time)
            newtime = self.gravity_code.model_time
           
            #JB, modified this, in Bonsai we can take a 0 time-step to detect
            #multiples. That would cause the newtime == time to evaluate to true
            #when there are multiples detected and break out of the evaluate
            #loop before the time reached end_time
            if newtime == time and (stopping_condition.is_set() == False):
                break
                
            time = newtime
            
            if stopping_condition.is_set():

                star1 = stopping_condition.particles(0)[0]
                star2 = stopping_condition.particles(1)[0]
                ignore = 0

                # Note from Steve, 8/12: We can pick up a lot of
                # encounters that are then ignored here.  I have
                # temporarily duplicated this check in the ph4
                # module (jdata.cc).

                # r = (star2.position-star1.position).length()
                # v = (star2.velocity-star1.velocity).length()
                # cos_angle = numpy.inner((star2.velocity-star1.velocity)/v,
                #                         (star2.position-star1.position)/r)
                # angle = numpy.arccos(cos_angle)
                # 
                # #if r < 0.5 * (star1.radius + star2.radius) \	# ???
                # #        or angle > (numpy.pi * 0.44):
                # 
                # # Proceed only if the stars are moving parallel or
                # # toward each other.  We assume all angles more than
                # # 80 degrees will need to check binary from.
                # 
                # if angle > (numpy.pi * 0.44):

                r = (star2.position-star1.position).length()
                v = (star2.velocity-star1.velocity).length()
                vr = numpy.inner(star2.velocity-star1.velocity,
                                 star2.position-star1.position)
                EPS = 0.001
                if vr < EPS*r*v:

                    print '\n'+'~'*60
                    print 'interaction at time', time

                    # Synchronize everything for now.  Later we can
                    # just synchronize neighbors if gravity supports
                    # that.  TODO

                    self.gravity_code.synchronize_model()
                
                    # Like synchronize.  We only should copy data from
                    # the particles and their neighbors.  TODO

                    self.channel_from_code_to_memory.copy()
                    
                    initial_energy = self.get_total_energy(self.gravity_code)

                    star1 = star1.as_particle_in_set(self._inmemory_particles)
                    star2 = star2.as_particle_in_set(self._inmemory_particles)

                    print 'initial top-level:', \
                        star1.id, '('+str(star1.radius)+')', \
                        star2.id, '('+str(star2.radius)+')'
                    if 1:
                        print '                   r =', r
                        print '                   v =', v
                        print '                   v.r =', vr
                    sys.stdout.flush()

                    # Do the scattering.

                    veto, dE_top_level_scatter, dphi_top, dE_mul, \
                        dphi_int, dE_int \
                        = self.manage_encounter(star1, star2, 
                                                self._inmemory_particles,
                                                self.gravity_code.particles,
                                                self.kepler)

                    if not veto:

                        # Recommit is done automatically and
                        # reinitializes all particles.  Later we will
                        # just reinitialize a list if gravity supports
                        # it. TODO
                    
                        self.gravity_code.particles.synchronize_to(
                            self._inmemory_particles)    # star = gravity_stars
                    
                        self.channel_from_code_to_memory.copy_attribute \
                            ("index_in_code", "id")
                    
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
                        #	dphi_top is the top-levle tidal error
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
                        #	error due to internal configuration
                        #	changes in the scattering system
                        #
                        #	dE_int is the integration error in the
                        #	scattering calculation
                        #
                        # We *always* expect
                        #
                        #   dE_top_level - dE_top_level_scatter - dphi_top = 0.
                        #
                        # If this is not the case, there is an error
                        # in the internal bookkeeping of
                        # manage_encounter().

                        if 0:
                            #print 'top-level initial energy =', initial_energy
                            #print 'top-level final energy =', final_energy
                            print 'dE_top_level =', dE_top_level
                            print 'dE_top_level_scatter =', dE_top_level_scatter
                            print 'dphi_top =', dphi_top
                            print 'dphi_int =', dphi_int
                            print 'dE_int =', dE_int

                        print 'net local error =', \
                            dE_top_level - dE_top_level_scatter - dphi_top
                        print 'scatter integration error =', dE_int

                        # We also expect
                        #
                        #	dE_top_level_scatter + dE_mul
                        #			= dphi_top + dE_int - dphi_int.
                        #
                        # Monitor this and keep track of the
                        # cumulative value of the right-hand side of
                        # this equation.

                        print 'dE_mul =', dE_mul
                        print 'internal local error =', \
                            dE_top_level + dE_mul - dphi_top
                        print 'corrected internal local error =', \
                            dE_top_level + dE_mul - dphi_top + dphi_int - dE_int

                        self.multiples_external_tidal_correction += dphi_top
                        self.multiples_internal_tidal_correction -= dphi_int
                        self.multiples_integration_energy_error += dE_int
                        Nmul, Nbin, Emul = self.get_total_multiple_energy2()

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

                        print 'total energy (top+mul) =', \
                            final_energy + Emul
                        print 'corrected total energy =', \
                            final_energy + Emul \
                                - self.multiples_external_tidal_correction \
                                - self.multiples_internal_tidal_correction \
                                - self.multiples_integration_energy_error
                        
                        # Print all multiples currently in the system.

                        self.print_multiples()
                        count_resolve_encounter += 1

                    else:
                        ignore = 1

                    print '~'*60
                    sys.stdout.flush()

                else:
                    ignore = 1

                count_ignore_encounter += ignore

        print ''
        print 'Resolved', count_resolve_encounter, 'encounters'
        print 'Ignored', count_ignore_encounter, 'encounters'
        sys.stdout.flush()

        self.gravity_code.synchronize_model()
        self.channel_from_code_to_memory.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    def manage_encounter(self, star1, star2, stars, gravity_stars, kep):

        # Manage an encounter between star1 and star2.  Stars is the
        # python memory data set.  Gravity_stars points to the gravity
        # module data.  Return values are the change in top-level
        # energy, the tidal error, and the integration error in the
        # scattering calculation.  Steps below follow those defined in
        # the PDF description.

        # find_binaries(stars, self.gravity_constant) 

        #----------------------------------------------------------------
        # 1a. Build a list of stars involved in the scattering.  Start
        # with star1 and star2.

        scattering_stars = datamodel.Particles(particles = (star1, star2))

        # 1b. Add neighbors if desired.  Use a simple distance
        # criterion for now -- possibly refine later.  Also consider a
        # simple neighbor veto.  Start by sorting all stars by
        # distance from the CM.  Later, use neighbors, if supported.
        # TODO

        distances = (stars.position
                      - scattering_stars.center_of_mass()).lengths()
        indices = numpy.argsort(distances.number)
        sorted_stars = stars[indices]
        sorted_distances = distances[indices]
        #print "sorted_distances", sorted_distances

        sep12 = ((star1.position-star2.position)**2).sum().sqrt()
        rad12 = star1.radius + star2.radius

        # Sep12 is the separation of the two original components.  It
        # should be slightly less than the sum of their radii, rad12,
        # but it may be less in unexpected circumstances or if vetoing
        # is in effect.  Initial_scale sets the "size" of the
        # interaction and the scale to which the final products will
        # be rescaled.  Rad12 is also the 90 degree scattering
        # sistance for the two stars, and hence the natural limit on
        # binary scale.  Rnnmax sets the distance inside which we
        # check for neighbors.

        if not self.neighbor_veto:
            initial_scale = self.initial_scale_factor * rad12
        else:
            initial_scale = self.initial_scale_factor * sep12

        rnnmax = self.neighbor_distance_factor * sep12

        print 'initial_scale =', initial_scale

        for i in range(2,len(sorted_distances)):
            star = sorted_stars[i]

            # Note that it is possible for star1 and star2 not to be
            # at the start of the sorted list.

            if sorted_distances[i] < rnnmax \
                    and star != star1 and star != star2:
                if not self.neighbor_veto:
                    scattering_stars.add_particle(star)
                    print 'added',
                    if hasattr(star, 'id'):
                        print 'star', star.id,
                    else:
                        print 'unknown star',
                    print 'to scattering list'
                    sys.stdout.flush()
		    #initial_scale = sorted_distances[i]    # don't expand!
                else:
                    print 'Encounter vetoed by neighbor', star.id, \
                          'at distance', sorted_distances[i]
                    return True, 0., 0., 0., 0., 0.

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
                                            stars - scattering_stars,
                                            G=self.gravity_constant)

        print_internal = False

        if print_internal:
            print 'E0 =', E0
            print 'phi_rem =', phi_rem

        # 2b. If there are no neighbors, separate star1 and star2 to
        #     some large "scattering" radius.  If neighbors exist,
        #     just start the "scattering" interaction in place.

        initial_scatter_scale = self.initial_scatter_factor * initial_scale

        if len(scattering_stars) == 2:
            rescale_binary_components(star1, star2, kep,
                                      initial_scatter_scale, compress=False)

        # 2c. Remove the interacting stars from the gravity module.

        for star in scattering_stars:
            gravity_stars.remove_particle(star)

        #----------------------------------------------------------------
        # 3a. Create a particle set to perform the close encounter
        #     calculation.

        particles_in_encounter = datamodel.Particles(0)

        # 3b. Add stars to the encounter set, add in components when we
        #    encounter a binary.

        Emul_init = 0.0 | nbody_system.energy

        for star in scattering_stars:
            if star in self.root_to_tree:
                tree = self.root_to_tree[star]
                isbin, dEmul = get_multiple_energy2(tree, self.gravity_constant)
                Emul_init += dEmul
                openup_tree(star, tree, particles_in_encounter)
                del self.root_to_tree[star]
            else:
                particles_in_encounter.add_particle(star)

        # Terminology from the PDF description:

        E1 = particles_in_encounter.kinetic_energy() + \
	       particles_in_encounter.potential_energy(G=self.gravity_constant)

        dphi_1 = E1 - E0 - Emul_init

        if print_internal:
            print 'E1 =', E1
            print 'Emul_init =', Emul_init
            print 'dphi_1 =', dphi_1

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

        M,a,e,r,E,tperi = get_component_binary_elements(star1, star2, 
                                                        self.kepler, 1)
        print 'semi =', a.number, ' ecc =', e, ' E/mu =', E.number
        print 'tperi =', tperi.number, ' period =', self.kepler.get_period()
        sys.stdout.flush()

        if self.gravity_code.unit_converter is None:
            end_time = 10000 | nbody_system.time
            delta_t = min(10*abs(tperi), 1.0 | nbody_system.time)
        else:
            end_time = 10000.0 * abs(tperi)
            delta_t = abs(tperi)

        final_scatter_scale = self.final_scatter_factor * initial_scatter_scale

        scatter_energy_error = self.resolve_collision(particles_in_encounter,
                                                      final_scatter_scale,
                                                      end_time, delta_t)

        # Note that on return, particles_in_encounter contains CM
        # nodes in the list.

        #E2CM = get_energy_of_leaves(particles_in_encounter,
        #                            G=self.gravity_constant)
        #print 'E2 (CM) =', E2CM

        # Note: final_scatter_scale is used to limit the extent of the
        # smallN integration: the integration ends when any particle
        # is more than final_scatter_scale from the CM of the system.
        # TODO: what if the encounter reaches final_scatter_scale and
        # isn't over?  Currently avoided by ignoring this parameter --
        # see resolve_collision() below.

        # May be premature to leave the CM frame here, as we will
        # reuse CM quantities during rescaling...

        particles_in_encounter.position += cmpos
        particles_in_encounter.velocity += cmvel

        # Terminology from the PDF description:

        E2 = get_energy_of_leaves(particles_in_encounter,
                                  G=self.gravity_constant)
        dE_int = E2 - E1	# should equal scatter_energy_error
        if abs(dE_int - scatter_energy_error).number > 1.e-12:
            print '*** warning: dE_int mismatch ***'
            print 'scatter_energy_error =', scatter_energy_error
            print 'dE_int =', dE_int

        if print_internal:
            print 'E2 =', E2
            print 'scatter_energy_error =', scatter_energy_error
            print 'dE_int =', dE_int

        #----------------------------------------------------------------
        # 5a. Identify multiple structure after the encounter.  First
        #     create an object to handle the new binary information.

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

        # Multiple centers of mass.
        roots_of_trees = binaries.roots()

        #----------------------------------------------------------------
        # 6a. Scale to a radius slightly larger than the initial one.
        # Rescaling does just that -- neither computes nor attempts to
        # absorb the tidal error.  If we want to absorb the tidal
        # error rather than simply recording it, do so after splitting
        # wide binaries below.  TODO

        final_scale = self.final_scale_factor * initial_scale

        scale_top_level_list(stars_not_in_a_multiple,
                             roots_of_trees,
                             self.kepler,
                             final_scale,
                             self.gravity_constant)

        # 6b. Break up wide top-level binaries.  Do this after
        #     rescaling because we want to preserve binary binding
        #     energies.  Also place the wide binaries at pericenter.

        # Number of top-level nodes.
        lt = len(stars_not_in_a_multiple) + len(roots_of_trees)

        modified_list = False

        for root in list(roots_of_trees):
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

            if binary_scale > rad12:

                print 'initial top-level:', \
                    comp1.id, '('+str(comp1.radius)+')', \
                    comp2.id, '('+str(comp2.radius)+')'
                print 'splitting wide binary', name_pair(comp1,comp2)
                print 'semi =', semi
                print 'apo =', apo
                print 'peri =', semi*(1-ecc)
                print 'E/mu =', E
                # print 'initial_scale =', initial_scale
                sys.stdout.flush()

                # See the "special case" logic in
                # scale_top_level_list().  If this is a sole bound
                # top-level object, it has already been scaled to the
                # desired separation and should *not* be modified
                # here.  Otherwise, move the components to periastron
                # to minimize the tidal error when they are
                # reinstated..

                if lt > 1:

                    # Could use rescale_binary_components() for this,
                    # but code here is more compact, since we have
                    # already initialized the kepler structure.

                    print 'moving binary to periastron'
                    sys.stdout.flush()
                    cmpos = root.position
                    cmvel = root.velocity
                    self.kepler.advance_to_periastron()

	            ### Question to Arjen: what is the right syntax to
	            ### do this??  We want to say
                    ###
                    ###     rel_pos = self.kepler.get_separation_vector()
                    ###     comp1.position = cmpos - f1*rel_pos
                    ###     etc.
                    ###
                    ### but this doesn't work...

                    x,y,z = self.kepler.get_separation_vector()
                    vx,vy,vz = self.kepler.get_velocity_vector()
                    f1 = comp1.mass/root.mass
                    rel_pos = numpy.array([x,y,z])
                    rel_vel = numpy.array([vx,vy,vz])
                    for k in range(3):
                        comp1.position[k] = cmpos[k] - f1*rel_pos[k]
                        comp1.velocity[k] = cmvel[k] - f1*rel_vel[k]
                        comp2.position[k] = cmpos[k] + (1-f1)*rel_pos[k]
                        comp2.velocity[k] = cmvel[k] + (1-f1)*rel_vel[k]

                # Removing the CM reinstates the children as top-level
                # objects:

                particles_in_encounter.remove_particle(root)
                modified_list = True

        if modified_list:

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

        Emul_final = 0.0 | nbody_system.energy
        for tree in binaries.iter_binary_trees():            
            isbin, dEmul = get_multiple_energy2(tree, self.gravity_constant)
            Emul_final += dEmul

        dphi_2 = E2 - Emul_final - E3

        if print_internal:
            print 'E3 =', E3
            print 'phi_ins =', phi_ins
            print 'Emul_final =', Emul_final
            print 'dphi_2 =', dphi_2

        # 7a. Set radii to reflect multiple structure.
            
        set_radii(particles_in_encounter, self.kepler)

        # Print diagnostics on added particles.

        print 'final top-level:',
        r = 0.0|nbody_system.length
        v = 0.0|nbody_system.speed
        vr = 0.0|nbody_system.length*nbody_system.speed
        for i in top_level_nodes:
            print i.id, '('+str(i.radius)+')',
            for j in top_level_nodes:
                if i.id > j.id:
                    rij = (i.position-j.position).length()
                    if rij > r:
                        r = rij
                        v = (i.velocity-j.velocity).length()
                        vr = numpy.inner(j.velocity-i.velocity,
                                         j.position-i.position)
        print ''
        if 1:
            print '                 r =', r
            print '                 v =', v
            print '                 v.r =', vr
        sys.stdout.flush()

        # Update the gravity module with the new data.

        # 7b. Add stars not in a binary to the gravity code.
        if len(stars_not_in_a_multiple) > 0:
            gravity_stars.add_particles(stars_not_in_a_multiple)
            
        # 7c. Add the roots to the gravity code
        multiples_particles = Particles()
        multiples_particles.id = None
        for tree in binaries.iter_binary_trees():
            tree.particle.id = assign_id_to_root(tree)
            gravity_stars.add_particle(tree.particle)
            multiples_particles.add_particle(tree.particle)

        # DEBUG
        print "multiples: interaction products: singles:", \
                stars_not_in_a_multiple.id, "multiples: ", \
                multiples_particles.id 
            
        # 7d. Store all trees in memory for later reference.

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
        if abs(dph) > 1.e-2:		# 1.e-2 is arbitrary
            print '*** tidal correction =', dph, 'KE ***'
            #print 'initial configuration: phi =', \
            #    potential_energy_in_field(scattering_stars, 
            #                              stars - scattering_stars,
            #                              G=self.gravity_constant)
            pminmin, fminmin, dxminmin \
		= find_nn(scattering_stars, stars-scattering_stars,
                          self.gravity_constant)
            if pminmin != None:
                print 'closest field/list pair is', \
                    str(fminmin.id)+'/'+str(pminmin.id), \
                    ' distance/scale =', dxminmin/initial_scale
            #print 'final configuration: phi =', \
            #    potential_energy_in_field(top_level_nodes, 
            #                              stars - scattering_stars,
            #                              G=self.gravity_constant)
            pminmin, fminmin, dxminmin \
		= find_nn(top_level_nodes, stars-scattering_stars,
                          self.gravity_constant)
            if pminmin != None:
                print 'closest field/list pair is', \
                    str(fminmin.id)+'/'+str(pminmin.id), \
                    ' distance/scale =', dxminmin/initial_scale
        #-------------------------------------------------------

        return False, dE_top, dphi_top, dEmul, dphi_int, dE_int

    def resolve_collision(self,
                          particles,
                          final_scatter_scale,
                          end_time = 1000 | nbody_system.time,
                          delta_t = 10 | nbody_system.time):

        # Take the system described by particles and evolve it forward
        # in time until it is over.  Don't update global quantities,
        # don't interpret the outcome.  Return the energy error due to
        # the smallN integration.

        # Temporarily avoid "is_over" problems.  If we allow
        # collisions to stop early -- when they become too large or
        # last too long -- then we need will logic to manage the
        # intermediate state that results.  TODO
        final_scatter_scale = 1.e30 | nbody_system.length

        resolve_collision_code = self.resolve_collision_code_creation_function()

        time = 0 * end_time
        sys.stdout.flush()
        resolve_collision_code.set_time(time);
        resolve_collision_code.particles.add_particles(particles)
        resolve_collision_code.commit_particles()

        # Channel to copy values from the code to the set in memory.
        channel = resolve_collision_code.particles.new_channel_to(particles)

        initial_scatter_energy = self.get_total_energy(resolve_collision_code)

        print "multiples: number_of_stars =", len(particles), ' ', particles.id
        print 'multiples: initial energy =', initial_scatter_energy
        #print particles
        print "multiples: evolving to time =", end_time, 
        print "in steps of", delta_t
        if self.debug_encounters:
            print 'multiples: ### START ENCOUNTER ###'
            print 'multiples: ### snapshot at time %f' % 0.0
            for p in particles:
                print 'multiples: ### id=%d, x=%f, y=%f, z=%f,'\
                      'vx=%f, vy=%f, vz=%f' % \
                        (p.id, p.x.number, p.y.number, p.z.number,
                         p.vx.number, p.vy.number, p.vz.number)
        sys.stdout.flush()

        resolve_collision_code.set_break_scale(final_scatter_scale)
        delta_t_max = 64*delta_t

        if self.debug_encounters:
            delta_t *= 0.1

        while time < end_time:

            time += delta_t
            print 'multiples: evolving to time', time
            sys.stdout.flush()

            resolve_collision_code.evolve_model(time)

            # DEBUGGING:
            if self.debug_encounters:
                print 'multiples: ### snapshot at time %f' % time.number
                #resolve_collision_code.update_particle_tree()
                #resolve_collision_code.update_particle_set()
                resolve_collision_code.particles.synchronize_to(particles)
                channel.copy()
                for p in particles:
                    print 'multiples: ### id=%d, x=%f, y=%f, z=%f,'\
                          'vx=%f, vy=%f, vz=%f' % \
                            (p.id, p.x.number, p.y.number, p.z.number,
                             p.vx.number, p.vy.number, p.vz.number)
                sys.stdout.flush()

            # The argument final_scatter_scale is used to limit the
            # size of the system.  It has to be supplied again because
            # the code that determines if the scattering is over isn't
            # necessarily the same as resolve_collision_code.
            # Currently, only smallN has an "is_over()" function.
            # TODO
            #
            # Return values:	0 - not over
            #			1 - over
            #			2 - not over, but size exceeded limit
            #
            # Note that this is really a stopping condition, and
            # should eventually be handled that way.  TODO

            # We are currently ignoring any possibility of a physical
            # collision during the multiples encounter.  TODO

            over = resolve_collision_code.is_over(final_scatter_scale)

            # TODO: what happens if we reach final_scatter_scale but
            # the encounter isn't over?

            if over:
                final_scatter_energy \
                    = self.get_total_energy(resolve_collision_code)
                scatter_energy_error \
                    = final_scatter_energy - initial_scatter_energy
                print 'multiples: over =', over, 'at time', time
                print 'multiples: initial energy =', initial_scatter_energy
                print 'multiples: final energy =', final_scatter_energy
                print 'multiples: energy error =', scatter_energy_error
                if self.debug_encounters:
                    print 'multiples: ### END ENCOUNTER ###'
                sys.stdout.flush()

                # Create a tree in the module representing the binary structure.
                resolve_collision_code.update_particle_tree()

                # Note that center of mass particles are now part of
                # the particle set...

                # Return the tree structure to AMUSE.  Children are
                # identified by get_children_of_particle in interface.??,
                # and the information is returned in the copy operation.

                resolve_collision_code.update_particle_set()
                resolve_collision_code.particles.synchronize_to(particles)
                #print "resolve_collision_code.particles.radius", \
                #       resolve_collision_code.particles.radius
                channel.copy()
                #resolve_collision_code.stop()

                return scatter_energy_error

            if not self.debug_encounters:
                if delta_t < delta_t_max and time > 0.999999*4*delta_t:
                    delta_t *= 2

        raise Exception(
            "Did not finish the small-N simulation before end time {0}".
            format(end_time)
        )
    
def openup_tree(star, tree, particles_in_encounter):

    # List the leaves.

    leaves = tree.get_leafs_subset()
      
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

def sep2(star1, star2):			# squared separation of star1 and star2
    return ((star1.position-star2.position).number**2).sum()

def sep(star1, star2):			# separation of star1 and star2
    return math.sqrt(sep2(star1, star2))

def phi_tidal(star1, star2, star3):	# compute tidal potential of
					# (star1,star2) relative to star3
    phi13 = -star1.mass*star3.mass/sep(star1,star3)
    phi23 = -star2.mass*star3.mass/sep(star2,star3)
    m12 = star1.mass + star2.mass
    f1 = star1.mass/m12
    cmx = (f1*star1.x+(1-f1)*star2.x).number
    cmy = (f1*star1.y+(1-f1)*star2.y).number
    cmz = (f1*star1.z+(1-f1)*star2.z).number
    phicm = -m12*star3.mass/math.sqrt((star3.x.number-cmx)**2
                                       + (star3.y.number-cmy)**2
                                       + (star3.z.number-cmz)**2)
    return (phi13+phi23-phicm).number

def find_nnn(star1, star2, stars):	# print next nearest neighbor
				        # of (star1, star2)
    top_level = stars

    min_dr = 1.e10
    id1 = star1.id
    id2 = star2.id
    for t in top_level:
        tid = t.id
        if tid != id1 and tid != id2:
            dr2 = sep2(t, star1)
            if dr2 > 0 and dr2 < min_dr:
                min_dr = dr2
                nnn = t
    min_dr = math.sqrt(min_dr)
    #print 'star =', int(id1), ' min_dr =', min_dr, \
    #      ' nnn =', int(nnn.id), '(', nnn.mass.number, ')'
    #print '    phi_tidal =', phi_tidal(star1, star2, nnn)
    #print '    nnn pos:', nnn.x.number, nnn.y.number, nnn.z.number
    #sys.stdout.flush()

    return nnn

def find_nn(plist, field, G):

    # Find and print info on the closest field particle (as
    # measured by potential) to any particle in plist.

    pminmin = None
    dminmin = 1.e30
    fminmin = None
    phiminmin = 0.0 | nbody_system.energy
    for f in field:
        dx = (plist.position - f.position).lengths()
        phi = -G*f.mass*plist.mass/dx
        phimin = 0.0 | nbody_system.energy
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
            print 'bound', p.id, particles[indices[1]].id, Emin

def potential_energy_in_field(particles, field_particles,
                              smoothing_length_squared = zero,
                              G=constants.G):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument field_particles: the external field consists of these (i.e. potential energy is calculated relative to the field particles) 
    :argument smoothing_length_squared: the smoothing length is added to every distance.
    :argument G: gravitational constant, need to be changed for particles in different units systems

    >>> from amuse.datamodel import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential_energy()
    quantity<-6.67428e-11 m**2 * kg * s**-2>
    """

    if len(field_particles) == 0:
        return zero * G		# ERROR: this is dimensionally incorrect - Steve

    sum_of_energies = zero
    for particle in particles:
        dx = particle.x - field_particles.x
        dy = particle.y - field_particles.y
        dz = particle.z - field_particles.z
        dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
        dr = (dr_squared+smoothing_length_squared).sqrt()
        m_m = particle.mass * field_particles.mass
        energy_of_this_particle = (m_m / dr).sum()
        sum_of_energies -= energy_of_this_particle

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
    # receding; otherwise it will be approaching.

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
    kep.initialize_from_dyn(mass,
                            rel_pos[0], rel_pos[1], rel_pos[2],
                            rel_vel[0], rel_vel[1], rel_vel[2])
    M,th = kep.get_angles()
    a,e = kep.get_elements()
    
    if (compress and sep12 > scale**2) \
            or (not compress and sep12 < scale**2):

        #print 'rescaling components', int(comp1.id), \
        #      'and', int(comp2.id), 'to separation', scale.number
        #sys.stdout.flush()
        
        #print 'a, e =', a, e
        if e < 1:
            peri = a*(1-e)
            apo = a*(1+e)
        else:
            peri = a*(e-1)
            apo = peri+a		# OK - used only to reset scale

        if compress:
            limit = peri + 1.e-4*(apo-peri)		# numbers are
            if limit > 1.1*peri: limit = 1.1*peri	# ~arbitrary
            if scale < limit:
                #print 'changed scale from', scale, 'to', limit
                scale = limit
            if M < 0:
                kep.advance_to_periastron()
                kep.advance_to_radius(scale)
            else:
                if kep.get_separation() < scale:
                    kep.advance_to_radius(scale)
                else:
                    kep.return_to_radius(scale)

	    # Note: Always end up on an outgoing orbit.  If
	    # periastron > scale, we are now just past periastron.

        else:
            limit = apo - 0.01*(apo-peri)
            if scale > limit:
                #print 'changed scale from', scale, 'to', limit
                scale = limit
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

        if hasattr(comp1, 'child1'):
            offset_particle_tree(comp1, newpos1-pos1, newvel1-vel1)
        if hasattr(comp2, 'child1'):
            offset_particle_tree(comp2, newpos2-pos2, newvel2-vel2)

    return a

def compress_nodes(node_list, scale):

    # Compress the top-level nodes in node_list to lie within diameter
    # scale.  Rescale velocities to conserve total energy (but
    # currently not angular momentum -- TODO).

    # Compute the center of mass position and velocity of the
    # top-level system.

    cmpos = node_list.center_of_mass()
    cmvel = node_list.center_of_mass_velocity()

    # Compute the size, potential, and kinetic energy of the system in
    # the center of mass frame.

    size = 0.0
    sepmin = 1.e10
    rijmin = 1.e10
    radsum = 0.0
    imin = -1
    jmin = -1
    pot = 0.0
    kin = 0.0
    for i in range(len(node_list)):
        m = node_list[i].mass.number
        rad = node_list[i].radius.number
        posi = node_list[i].position
        pos = (posi-cmpos).number
        vel = (node_list[i].velocity-cmvel).number
        r2 = numpy.inner(pos,pos)
        if r2 > size: size = r2
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

def print_elements(s, a, e, r, Emu, E):
    print '{0} elements  a = {1}  e = {2}  r = {3}  E/mu = {4}  E = {5}'\
	.format(s, a, e, r, Emu, E)
    sys.stdout.flush()

def print_multiple(m, kep, level=0):		##### not working? #####

    # Recursively print the structure of (multiple) node m.

    print '    '*level, 'key =', m.key, ' id =', int(m.id)
    print '    '*level, '  mass =', m.mass.number
    print '    '*level, '  pos =', m.position.number
    print '    '*level, '  vel =', m.velocity.number
    sys.stdout.flush()
    if not m.child1 is None and not m.child2 is None:
        M,a,e,r,E,t = get_component_binary_elements(m.child1, m.child2, kep)
        print_elements('   '+'    '*level+'binary', a, e, r, E,
                       (E*m.child1.mass*m.child2.mass/M))
    if not m.child1 is None:
        print_multiple(m.child1, kep, level+1)
    if not m.child2 is None:
        print_multiple(m.child2, kep, level+1)

def print_pair_of_stars(s, star1, star2, kep):
    m1 = star1.mass
    m2 = star2.mass
    M,a,e,r,E,t = get_component_binary_elements(star1, star2, kep)
    print_elements(s, a, e, r, E, E*m1*m2/(m1+m2))
    print_multiple(star1)
    print_multiple(star2)
    
def print_simple_multiple(node, kep):
    for level, x in node.iter_levels():
        output = ''
        if level == 0: output += 'Multiple '
        output += '    ' * level
        particle = x
        output += "{0} id = {1} mass = {2}".format(particle.key,
                                                   particle.id,
                                                   particle.mass.number)
        if not particle.child1 is None:
            child1 = particle.child1
            child2 = particle.child2
            M,a,e,r,E,t = get_component_binary_elements(child1, child2, kep)
            mu = child1.mass*child2.mass/M
            output += " semi = {0} energy = {1}".format(a.number,
                                                        (mu*E).number)
        print output
        sys.stdout.flush()

def another_print_multiple(node, kep, pre, kT, dcen):
    is_bin = 1
    Etot = 0.0 | nbody_system.energy
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
            print '%s%d (%d,%d) m=%.5f' % (init, id, idlow, idhigh, M),
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
            print 'a=%.6f e=%4f r=%6f D=%.4f E/mu=%.5f E=%.5f E/kT=%.5f' % \
		(a.number, e, r.number, D, Emu.number, E.number, E/kT)
            sys.stdout.flush()
        else:
            print '%s%d m=%.5f' % (init, id, M)
            sys.stdout.flush()

    return is_bin, Etot

def get_multiple_energy(node, kep):

    # Return the binary status and the total energy Etot of the
    # specified tree node.  The value of Etot is the total pairwise
    # energy of all binary objects in the hierarchy.  It does not
    # include the tidal potential of one component on another (e.g.
    # in a hierarchical triple Etot will be the sum of two binary
    # energies only).

    is_bin = 1
    Etot = 0.0 | nbody_system.energy
    for level, x in node.iter_levels():
        particle = x
        M = particle.mass.number
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

    # Return the binary status and the total energy Etot2 of the
    # specified tree node.  The returned value is the total energy
    # of all leaves in the hierarchy, properly including tidal
    # potentials, but excluding the center of mass energy.

    is_bin = 1
    Ecm = 0.0 | nbody_system.energy

    # List the leaves, and do some additional work.  Note that
    # node.get_leafs_subset() seems to do the same thing...

    leaves_in_node = datamodel.Particles(0)
    for level, x in node.iter_levels():
        particle = x
        if level == 0:		# top-level kinetic energy
            Ecm = 0.5*particle.mass*(particle.velocity**2).sum()
        if not particle.child1 is None:
            if level > 0: is_bin = 0
        else:			# list leaves
            leaves_in_node.add_particle(particle)
    return is_bin, leaves_in_node.kinetic_energy() \
			+ leaves_in_node.potential_energy(G=G) \
			- Ecm

def add_leaves(node, leaf_list):
    if node.child1 == None:
        leaf_list.add_particle(node)
    else:
        add_leaves(node.child1, leaf_list)
        add_leaves(node.child2, leaf_list)

def get_multiple_energy3(root, G):

    # Return the binary status and the total energy Etot2 of the
    # particles under the specified root, assuming that child pointers
    # exist and are correctly set.  The returned value is the total
    # energy of all leaves in the hierarchy, properly including tidal
    # potentials, but excluding the center of mass energy.

    # Not currently used...

    is_bin = 1
    Ecm = 0.0 | nbody_system.energy
    leaves_in_root = datamodel.Particles(0)
    add_leaves(root, leaves_in_root)
    Ecm = 0.5*root.mass*(root.velocity**2).sum()
    return leaves_in_root.kinetic_energy() \
			+ leaves_in_root.potential_energy(G=G) \
			- Ecm

def get_energy_of_leaves(particles, G):
    leaves = datamodel.Particles(0)
    for p in particles:
        if not hasattr(p, 'child1') or p.child1 == None:
            leaves.add_particle(p)
    return leaves.kinetic_energy() + leaves.potential_energy(G=G)

def print_energies(stars):

    # Brute force N^2 over top level, pure python...

    top_level = stars.select(is_not_a_child, ["is_a_child"])

    mass = 0
    kinetic = 0
    potential = 0
    for t in top_level:
        m = t.mass.number
        x = t.x.number
        y = t.y.number
        z = t.z.number
        vx = t.vx.number
        vy = t.vy.number
        vz = t.vz.number
        mass += m
        kinetic += 0.5*m*(vx**2+vy**2+vz**2)
        dpot = 0
        for tt in top_level:
            if tt != t:
                mm = tt.mass.number
                xx = tt.x.number-x
                yy = tt.y.number-y
                zz = tt.z.number-z
                dpot -= mm/math.sqrt(xx**2+yy**2+zz**2)
        potential += 0.5*m*dpot
            
    print 'len(stars) =', len(stars)
    print 'len(top_level) =', len(top_level)
    print 'mass =', mass
    print 'kinetic =', kinetic
    print 'potential =', potential
    print 'energy =', kinetic+potential
    sys.stdout.flush()

# def set_radius_recursive(node, kep):
#
#     if node.is_leaf(): return		# nothing to be done
#
#     # Propagate child radii upward.
#
#     rmax = zero
#     for child in node.iter_children():
#         set_radius_recursive(child, kep)
#         rmax = max(rmax, child.particle.radius)
#
#     # Include binary information.
#
#     node.particle.radius = rmax
#     try:
#         if not node.particle.child1 == None:
#             mass,a,e,r,E,t = get_cm_binary_elements(node.particle, kep)
#             if e < 1:
#                 node.particle.radius = max(2*a, node.particle.radius)
#     		#			   2 here is ~arbitrary
#     except:
#         pass

def set_radius_recursive(node, kep):

    if node.is_leaf(): return		# nothing to be done

    # Propagate child radii upward.  Since dynamical radius scales
    # with mass, the radius of a parent is the sum of the radii of the
    # children.  If we are handling 2-body encounters, that's all we
    # need.  The semi-major axis of a hard binary is irrelevant (less
    # than the dynamical radius, by definition)...

    rsum = zero
    for child in node.iter_children():
        set_radius_recursive(child, kep)
        rsum += child.particle.radius
    node.particle.radius = rsum

def set_radii(top_level_nodes, kep):
    for n in top_level_nodes.as_binary_tree().iter_children():
        set_radius_recursive(n, kep)

def scale_top_level_list(singles, multiples, kep, scale,
                         gravity_constant):

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

    print 'scale_top_level_list: ls =', ls, ' lm =', lm, ' lt =', lt
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

            print "scale_top_level_list: bound binary node"
            #print '\nunscaled binary node:'
            #print_multiple(root)
            comp1 = root.child1
            comp2 = root.child2
            #print "scale:", scale
            semi = rescale_binary_components(comp1, comp2, kep, scale)
            #true, mean = kep.get_angles()
            #print 'true =', true, 'mean =', mean
            #print 'scaled binary node:'
            #print_multiple(root, kep)

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

        print "scale_top_level_list: top-level unbound pair"
        #print '\nunscaled top-level pair:'
        #print_pair_of_stars('pair', comp1, comp2)
        semi = rescale_binary_components(comp1, comp2, kep, scale)
        #print '\nscaled top-level pair:'
        #print_pair_of_stars('pair', comp1, comp2)

    else:

        # Now we have three or more top-level nodes.  We don't know
        # how to compress them in a conservative way.  For now, we
        # will conserve energy and think later about how to preserve
        # angular momentum.  TODO

        print 'scale_top_level_list:', lt, 'top-level nodes'
        #print lt, 'unscaled top-level nodes'
        #print top_level_nodes
        compress_nodes(top_level_nodes, scale)
        #print lt, 'scaled top-level nodes'
        #print top_level_nodes

    sys.stdout.flush()

    # Don't attempt to correct or even return the tidal energy error.
    # Manage all of this in the calling function, as desired.

    return
