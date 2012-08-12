import sys
import numpy
import collections
import math

from amuse.community.ph4.interface import ph4 as grav
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.datamodel import trees
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import zero

root_index = 1000
def new_root_index():
    global root_index
    root_index += 1
    return root_index
    
def is_a_parent(child1_key, child2_key):
    return child1_key > 0 or child2_key > 0

def is_not_a_child(is_a_child):
    return is_a_child == 0

def print_log(s, gravity, E0 = zero):
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    try:
        Ebin = gravity.get_binary_energy()
    except:
        Ebin = zero 
    Etop = T + U
    E = Etop + Ebin
    if E0 == zero: E0 = E
    Rv = -0.5*M*M/U
    Q = -T/U
    print ""
    print "%s: time = %.10f, mass = %.10f, dE/E0 = %.5e" \
          % (s, gravity.get_time().number, M.number, (E/E0 - 1))
    print "%senergies = %.10f %.10f %.10f" \
          % (' '*(2+len(s)), E.number, U.number, T.number)
    #print '%s %.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f' % \
    #	(s+"%%", time.number, M.number, T.number, U.number, \
    #    E.number, Ebin.number, Rv.number, Q.number)
    sys.stdout.flush()
    return E

def get_component_binary_elements(comp1, comp2, peri = 0):
    try:
        comp1.mass.value_in(units.kg) # see if SI units, throw exception if not
        unit_converter \
            = nbody_system.nbody_to_si(comp1.mass + comp2.mass,
                                       (comp2.position-comp1.position).length())
    except Exception as ex:
        unit_converter = None
        
    kep = Kepler(unit_converter, redirection = "none")
    kep.initialize_code()

    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()	# per unit reduced mass, note
    if peri:
        M,th = kep.get_angles()
        if M < 0:
            kep.advance_to_periastron()
        else:
            kep.return_to_periastron()
    t = kep.get_time()
    kep.stop()

    return mass,a,e,r,E,t

def get_cm_binary_elements(p, peri = 0):
    return get_component_binary_elements(p.child1, p.child2, peri)

class Multiples(object):

    def __init__(self, gravity_code, resolve_collision_code_creation_function,
                 gravity_constant = None, **options):
        self.gravity_code = gravity_code
        self._inmemory_particles = self.gravity_code.particles.copy()
        self._inmemory_particles.id = self.gravity_code.particles.index_in_code
                
        self.channel_from_code_to_memory = \
            self.gravity_code.particles.new_channel_to(self._inmemory_particles)
        
        self.resolve_collision_code_creation_function \
            = resolve_collision_code_creation_function
        self.multiples_energy_correction \
            = zero * self.gravity_code.kinetic_energy
        self.root_to_tree = {}
        if gravity_constant is None:
            gravity_constant = nbody_system.G
        
        self.multiples = datamodel.Particles()
        self.gravity_constant = gravity_constant
        
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
        
    # Add?: def gravity.commit_particles()
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
            binaries_energy = code.get_binary_energy()
        except:
            binaries_energy = zero 
        total_energy = code.potential_energy + code.kinetic_energy \
				             + binaries_energy
        return total_energy
        
    def evolve_model(self, end_time):

        stopping_condition = self.gravity_code.stopping_conditions.collision_detection
        stopping_condition.enable()
        
        initial_energy = self.get_total_energy(self.gravity_code)
        
        time = self.gravity_code.model_time
        print "Evolve model to:", end_time, " starting at time:", time
        sys.stdout.flush()

        count_resolve_encounter = 0
        count_ignore_encounter = 0

        rlimit = 8./len(self.gravity_code.particles) | nbody_system.length  # 4 r_90 - TODO
        print 'rlimit =', rlimit
        
        while time < end_time:

            print 'calling evolve_model to ', end_time
            sys.stdout.flush()
            self.gravity_code.evolve_model(end_time)
            newtime = self.gravity_code.model_time
            
            if newtime == time:
                break
                
            time = newtime
            
            if stopping_condition.is_set():

                star1 = stopping_condition.particles(0)[0]
                star2 = stopping_condition.particles(1)[0]

                # Note from Steve, 8/12: We pick up a lot of
                # encounters that are then ignored here.  I have
                # temporarily duplicated this check in the ph4 module.

                r = numpy.sqrt(numpy.inner(star2.position-star1.position,
                                           star2.position-star1.position))
                v = numpy.sqrt(numpy.inner(star2.velocity-star1.velocity,
                                           star2.velocity-star1.velocity))
                vr = numpy.inner(star2.position-star1.position,
                                 star2.velocity-star1.velocity)

                EPS = 0.01		# proceed only if the stars are
                if vr < -EPS*r*v:	# clearly approaching

                    print '\n'+'~'*60
                    print 'interaction at time', time
                    print 'top-level: r =', r.number, ' v =', v.number, \
			  ' v.r =', vr.number
                    sys.stdout.flush()

                    energy = self.get_total_energy(self.gravity_code)

                    # Synchronize everything for now.  Later we can just
                    # synchronize neighbors if gravity supports that.
                    # TODO
                    self.gravity_code.synchronize_model()
                
                    # Like synchronize.
                    # We only should copy data from the particles and their
                    # neighbors.  TODO
                    self.channel_from_code_to_memory.copy()
                    
                    star1 = star1.as_particle_in_set(self._inmemory_particles)
                    star2 = star2.as_particle_in_set(self._inmemory_particles)
                    print 'star1 =', star1.id, ' star2 =', star2.id
                    sys.stdout.flush()

                    self.manage_encounter(star1, star2, 
                                          self._inmemory_particles,
                                          self.gravity_code.particles,
                                          rlimit)

                    # Recommit reinitializes all particles (redundant
                    # here, since it is done automatically).  Later we
                    # will just recommit and reinitialize a list if
                    # gravity supports it. TODO
                    # self.gravity_code.recommit_particles()
                    
                    self.gravity_code.particles.synchronize_to(
                        self._inmemory_particles) # make star = gravity_stars
                    
                    self.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                    
                    energy = self.get_total_energy(self.gravity_code)
                    print "multiples energy correction =", \
                        self.multiples_energy_correction, \
                        ' dE =', energy - initial_energy \
                                 - self.multiples_energy_correction

                    self.print_multiples()

                    print '~'*60
                    sys.stdout.flush()
                    count_resolve_encounter += 1

                else:
                    count_ignore_encounter += 1

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

    def print_multiples(self):
        for x in self.root_to_tree.values():
            print_simple_multiple(x)
            
    def pretty_print_multiples(self, pre, kT):
        Nbin = 0
        Nmul = 0
        Emul = 0.0 | nbody_system.energy
        for x in self.root_to_tree.values():
            Nmul += 1
            nb,E = another_print_multiple(x, pre, kT)
            Nbin += nb
            Emul += E
        return Nmul, Nbin, Emul
            
    def manage_encounter(self, star1, star2, stars, gravity_stars, rlimit):

        # Manage an encounter between star1 and star2.  stars is the
        # python memory data set.  gravity_stars points to the gravity
        # module data.  Return value is the energy correction due to
        # multiples.

        # Currently we don't include neighbors in the integration
        # and the size limitation on any final multiple is poorly
        # implemented.

        #print '\nin manage_encounter'; sys.stdout.flush()
        #print ''
        #find_nnn(star1_in_memory, star2_in_memory, stars)
        #find_nnn(star2_in_memory, star1_in_memory, stars)

        # 1. Calculate the potential of [star1, star2] relative to the
        #    other top-level objects in the stars list (later: just
        #    use neighbors TODO).  Start the correction of dEmult by
        #    removing the initial top-level energy of the interacting
        #    particles from it.  (Add in the final energy later, on
        #    return from scale_top_level_list.)

        collided_stars = datamodel.Particles(particles = (star1, star2))
               
        total_energy_of_stars_to_remove  = collided_stars.kinetic_energy()
        total_energy_of_stars_to_remove += collided_stars.potential_energy(G=self.gravity_constant)

        phi_in_field_of_stars_to_remove = potential_energy_in_field(
            collided_stars, 
            stars - collided_stars,
            G = self.gravity_constant
        )

        # 2. Create a particle set to perform the close encounter
        #    calculation.
        
        particles_in_encounter = datamodel.Particles(0)

        # 3. Add stars to the encounter set, add in components when we
        #    encounter a binary.

        if star1 in self.root_to_tree:
            openup_tree(star1, self.root_to_tree[star1], particles_in_encounter)
            del self.root_to_tree[star1]
        else:
            particles_in_encounter.add_particle(star1)
            
        if star2 in self.root_to_tree:
            openup_tree(star2, self.root_to_tree[star2], particles_in_encounter)
            del self.root_to_tree[star2]
        else:
            particles_in_encounter.add_particle(star2)

        #print '\nparticles in encounter (flat tree):'
        #for p in particles_in_encounter:
        #    print_multiple(p)

        initial_scale = (star1.position - star2.position).length()

        # 4. Run the small-N encounter in the center of mass frame.

        # May be desirable to make this a true scattering experiment
        # by separating star1 and star2 to a larger distance before
        # starting the multiples code.  TODO

        total_mass = collided_stars.mass.sum()
        cmpos = collided_stars.center_of_mass()
        cmvel = collided_stars.center_of_mass_velocity()
        particles_in_encounter.position -= cmpos
        particles_in_encounter.velocity -= cmvel

        M,a,e,r,E,tperi = get_component_binary_elements(star1, star2, 1)
        print 'semi =', a.number, ' ecc =', e, ' energy =', E.number, \
              ' tperi =', tperi.number
        sys.stdout.flush()
        if self.gravity_code.unit_converter is None:
            end_time = 10000 | nbody_system.time
            delta_t = min(10*abs(tperi), 1.0 | nbody_system.time)
        else:
            end_time = 10000.0 * abs(tperi)
            delta_t = abs(tperi)

        self.resolve_collision(particles_in_encounter, rlimit, end_time, delta_t)
       
        particles_in_encounter.position += cmpos
        particles_in_encounter.velocity += cmvel
            
        #print '\nparticles in encounter (after resolve):'
        #print particles_in_encounter.to_string(['x','y','z',
        #                                        'vx','vy','vz',
        #                                        'mass', 'id',
        #                                        'child1', 'child2'])
        #sys.stdout.flush()

        # 5. Carry out bookkeeping after the encounter and update the
        #    gravity module with the new data.
        
        # 5a. Remove star1 and star2 from the gravity module.

        gravity_stars.remove_particle(star1)
        gravity_stars.remove_particle(star2)
        
        # 5b. Create an object to handle the new binary information.

        binaries = trees.BinaryTreesOnAParticleSet(particles_in_encounter,
                                                   "child1", "child2")

        # 5bb. Compress the top-level nodes before adding them to the
        #      gravity code.  Also recompute the external potential and
        #      absorb the tidal error into the top-level nodes of the
        #      encounter list.  Fially, add the change in top-level energy
        #      of the interacting subset into dEmult, so E(ph4) + dEmult
        #      should be conserved.

        stars_not_in_a_multiple = binaries.particles_not_in_a_multiple()
        roots_of_trees = binaries.roots()
        
        # Set radii to reflect multiple structure.  This is probably not
        # the best place to do it...
            
        set_radii(particles_in_encounter)

        total_energy_of_stars_to_add, phi_correction \
            = scale_top_level_list(stars_not_in_a_multiple,
                                   roots_of_trees,
                                   initial_scale,
                                   stars - collided_stars, 
                                   phi_in_field_of_stars_to_remove,
                                   self.gravity_constant)

        # 5d. Add stars not in a binary to the gravity code.
        if len(stars_not_in_a_multiple) > 0:
            gravity_stars.add_particles(stars_not_in_a_multiple)
            
        # 5e. Add the roots to the gravity code
        for tree in binaries.iter_binary_trees():
            tree.particle.id = new_root_index()
            gravity_stars.add_particle(tree.particle)
            
        # 5f. Store all trees in memory for later reference
        for tree in binaries.iter_binary_trees():            
            self.root_to_tree[tree.particle] = tree.copy()
            
        self.multiples_energy_correction \
            += (total_energy_of_stars_to_add
                 - total_energy_of_stars_to_remove) + phi_correction

        #print '\nleaving manage_encounter'; sys.stdout.flush()
    
    def resolve_collision(
            self,
            particles,
            rlimit,
            end_time = 1000 | nbody_system.time,
            delta_t = 10 | nbody_system.time,
        ):

        #print '\nin resolve_collision'; sys.stdout.flush()

        resolve_collision_code = self.resolve_collision_code_creation_function()
        time = 0  * end_time
        #print "\nadding particles to multiples code"
        sys.stdout.flush()
        resolve_collision_code.set_time(time);
        resolve_collision_code.particles.add_particles(particles)
        #print "committing particles to multiples code"
        resolve_collision_code.commit_particles()

        # Channel to copy values from the code to the set in memory.
        channel = resolve_collision_code.particles.new_channel_to(particles)

        initial_energy = self.get_total_energy(resolve_collision_code)

        print "multiples: number_of_stars =", len(particles), ' ', particles.id
        print 'multiples: initial energy =', initial_energy
        #print particles
        print "multiples: evolving to time =", end_time, 
        print "in steps of", delta_t
        sys.stdout.flush()

        delta_t_max = 64*delta_t
        break_scale = rlimit
        resolve_collision_code.set_break_scale(break_scale)

        while time < end_time:

            time += delta_t
            print 'multiples: evolving to time', time
            sys.stdout.flush()

            resolve_collision_code.evolve_model(time)

            # TODO: currently, only smallN has an "is_over()" function.
            over = resolve_collision_code.is_over(break_scale)

            if over:
                print 'multiples: interaction is over at time', time
                energy = self.get_total_energy(resolve_collision_code)
                print 'multiples: final energy =', energy
                sys.stdout.flush()

                # Create a tree in the module representing the binary structure.
                resolve_collision_code.update_particle_tree()

                # Return the tree structure to AMUSE.  Children are
                # identified by get_children_of_particle in interface.??,
                # and the information is returned in the copy operation.

                resolve_collision_code.update_particle_set()
                resolve_collision_code.particles.synchronize_to(particles)
                #print "resolve_collision_code.particles.radius", \
                #      resolve_collision_code.particles.radius
                channel.copy()
                resolve_collision_code.stop()

                return initial_energy, energy

            if delta_t < delta_t_max and time > 0.999999*4*delta_t:
                delta_t *= 2

        resolve_collision_code.stop()
        raise Exception(
            "Did not finish the small-N simulation before end time {0}".
            format(end_time)
        )
    
    def print_trees_summary(self):
        if len(self.root_to_tree) > 0:
            print "number of multiples: ", len(self.root_to_tree)
            sys.stdout.flush()

def openup_tree(star, tree, particles_in_encounter):

    # List the leaves.

    leaves = tree.get_leafs_subset()
      
    original_star = tree.particle

    # Compare with the position stored when replacing the particles
    # with the root particle, and move the particles accordingly.
    # Note that once the CM is in the gravity module, the components
    # are frozen and the coordinates are absolute, so we need the
    # original coordinates to offset them later.

    # Maybe better just to store relative coordinates.  TODO.

    dx = star.x - original_star.x
    dy = star.y - original_star.y
    dz = star.z - original_star.z
    dvx = star.vx - original_star.vx
    dvy = star.vy - original_star.vy
    dvz = star.vz - original_star.vz
    
    #print "delta positions:", dx, dy, dz
    
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

def potential_energy_in_field(particles, field_particles,
                              smoothing_length_squared = zero, G = constants.G):
    """
    Returns the total potential energy of the particles in the particles set.

    :argument field_particles: the external field consists of these (i.e. potential energy is calculated relative to the field particles) 
    :argument smooting_length_squared: the smoothing length is added to every distance.
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
        return zero * G
        
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

def compress_binary_components(comp1, comp2, scale):

    # Compress the two-body system consisting of comp1 and comp2 to
    # lie within distance scale of one another.

    pos1 = comp1.position
    pos2 = comp2.position
    sep12 = ((pos2-pos1)**2).sum()

    if sep12 > scale*scale:
        #print '\ncompressing components', int(comp1.id), \
        #      'and', int(comp2.id), 'to separation', scale.number
        #sys.stdout.flush()
        mass1 = comp1.mass
        mass2 = comp2.mass
        total_mass = mass1 + mass2
        vel1 = comp1.velocity
        vel2 = comp2.velocity
        cmpos = (mass1*pos1+mass2*pos2)/total_mass
        cmvel = (mass1*vel1+mass2*vel2)/total_mass

        # For now, create and delete a temporary kepler process to
        # handle the transformation.  Would be more efficient to
        # define a single kepler at the start of the calculation and
        # reuse it.
        try:
            comp1.mass.value_in(units.kg) # see if SI units; exception if not
            unit_converter = nbody_system.nbody_to_si(
	                comp1.mass + comp2.mass,
                        (comp2.position - comp1.position).length())
        except Exception as ex:
            unit_converter = None
            
        kep = Kepler(unit_converter = unit_converter, redirection = "none")
        kep.initialize_code()
        mass = comp1.mass + comp2.mass
        rel_pos = pos2 - pos1
        rel_vel = vel2 - vel1
        kep.initialize_from_dyn(mass,
                                rel_pos[0], rel_pos[1], rel_pos[2],
                                rel_vel[0], rel_vel[1], rel_vel[2])
        M,th = kep.get_angles()
        a,e = kep.get_elements()
        #print 'a, e =', a, e
        if e < 1:
            peri = a*(1-e)
            apo = a*(1+e)
        else:
            peri = a*(e-1)
            apo = peri+a		# OK - used only to reset scale
        limit = peri + 0.01*(apo-peri)
        if scale < limit:
            #print 'changed scale from', scale, 'to', limit
            scale = limit

        #print 'scale =', scale, 'sep =', kep.get_separation(), 'M =', M
        #kep.print_all()

        if M < 0:
            kep.advance_to_periastron()
            kep.advance_to_radius(scale)
        else:
            if kep.get_separation() < scale:
                kep.advance_to_radius(scale)
            else:
                kep.return_to_radius(scale)

        # Note: if periastron > scale, we are now just past periastron.

        new_rel_pos = kep.get_separation_vector()
        new_rel_vel = kep.get_velocity_vector()
        kep.stop()

        # Enew = 0
        # r2 = 0
        # for k in range(3):
        #     Enew += 0.5*(new_rel_vel[k].number)**2
        #     r2 += (new_rel_pos[k].number)**2
        # rnew = math.sqrt(r2)
        # Enew -= mass.number/r1
        # print 'E, Enew, rnew =', E.number, E1, r1

        # Problem: the vectors returned by kepler are lists,
        # not numpy arrays, and it looks as though we can say
        # comp1.position = pos, but not comp1.position[k] =
        # xxx, as we'd like...  Also, we don't know how to
        # copy a numpy array with units...  TODO

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
        # transmit them to the (currently absolute) coordinates of
        # all lower components.

        offset_particle_tree(comp1, newpos1-pos1, newvel1-vel1)
        offset_particle_tree(comp2, newpos2-pos2, newvel2-vel2)

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

def print_multiple(m, level=0):		##### not working? #####

    # Recursively print the structure of (multiple) node m.

    print '    '*level, 'm =', m.key
    print '    '*level, ' ', int(m.id), ' mass =', m.mass.number
    print '    '*level, '  pos =', m.position.number
    print '    '*level, '  vel =', m.velocity.number
    sys.stdout.flush()
    if not m.child1 is None and not m.child2 is None:
        M,a,e,r,E,t = get_component_binary_elements(m.child1, m.child2)
        print_elements('   '+'    '*level+'binary', a, e, r, E,
                       (E*m.child1.mass*m.child2.mass/M))
    if not m.child1 is None:
        print_multiple(m.child1, level+1)
    if not m.child2 is None:
        print_multiple(m.child2, level+1)

def print_pair_of_stars(s, star1, star2):
    m1 = star1.mass
    m2 = star2.mass
    M,a,e,r,E,t = get_component_binary_elements(star1, star2)
    print_elements(s, a, e, r, E, E*m1*m2/(m1+m2))
    print_multiple(star1)
    print_multiple(star2)
    
def print_simple_multiple(node):
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
            M,a,e,r,E,t = get_component_binary_elements(child1, child2)
            mu = child1.mass*child2.mass/M
            output += " semi = {0} energy = {1}".format(a.number,
                                                        (mu*E).number)
        print output
        sys.stdout.flush()

def another_print_multiple(node, pre, kT):
    is_bin = 1
    Etot = 0.0 | nbody_system.energy
    for level, x in node.iter_levels():
        particle = x
        init = pre
        if level == 0: init += 'Multiple '
        init += '    ' * level
        id = particle.id
        M = particle.mass.number
        print '%s%d m=%.5f' % (init, id, M),
        if not particle.child1 is None:
            if level > 0: is_bin = 0
            child1 = particle.child1
            child2 = particle.child2
            M,a,e,r,Emu,t = get_component_binary_elements(child1, child2)
            mu = child1.mass*child2.mass/M
            E = Emu*mu
            Etot += E
            print 'a=%.6f e=%4f r=%6f E/mu=%.5f E=%.5f E/kT=%.5f' % \
		(a.number, e, r.number, Emu.number, E.number, E/kT),
        print ''
        sys.stdout.flush()
    return is_bin, Etot

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

def set_radius_recursive(node):
    if node.is_leaf(): return

    # Propagate child radii upward.

    rmax = zero
    for child in node.iter_children():
        set_radius_recursive(child)
        rmax = max(rmax, child.particle.radius)

    # Include binary information.

    node.particle.radius = rmax
    try:
        if not node.particle.child1 == None:
            mass,a,e,r,E,t = get_cm_binary_elements(node.particle)
            if e < 1:
                node.particle.radius = max(3*a, node.particle.radius)
    		# 3 is ~arbitrary
    except:
        pass

def set_radii(top_level_nodes):
    for n in top_level_nodes.as_binary_tree().iter_children():
        set_radius_recursive(n)

def scale_top_level_list(
        singles, multiples, 
        scale,
        field, phi_in_field_of_stars_to_remove, gravity_constant):

    # The multiples code followed the particles until their
    # interaction could be unambiguously classified as over.  They may
    # now be very far apart.  Input binaries is an object describing
    # the final hierarchical structure of the interacting particles in
    # the multiples code.  It consists of a flat tree of binary trees.
    # TODO: this is quite specific to smallN.

    # Scale the positions and velocities of the top-level nodes to
    # bring them within a sphere of diameter scale, conserving energy
    # and angular momentum (if possible).  Also offset all children to
    # reflect changes at the top level -- TODO will change when
    # offsets are implemented...

    # We are currently ignoring any possibility of a physical
    # collision during the multiples encounter.  TODO

    # Logic: 1 node   - must be a binary, use kepler to reduce to scale
    #        2 nodes  - use kepler, reduce binary children too? TODO
    #        3+ nodes - shrink radii and rescale velocities to preserve
    #                   energy, but this doesn't preserve angular
    #                   momentum TODO - also reduce children? TODO

    top_level_nodes = singles + multiples

    # Figure out the tree structure.

    ls = len(singles)
    lm = len(multiples)
    lt = ls + lm

    print 'scale_top_level_list: ls =', ls, ' lm =', lm, ' lt =', lt
    sys.stdout.flush()

#     if lt == 1 and lm == 1:

#         # Check if we want to change the status of a wide binary.

#         a_max = 0.1 | nbody_system.length	# TODO

#         root = multiples[0]
#         comp1 = root.child1
#         comp2 = root.child2
#         M,a,e,r,E,t = get_component_binary_elements(comp1, comp2)
#         if a > a_max:

#             singles.add_particle(comp1)
#             singles.add_particle(comp2)
#             #multiples.delete_particle(root)	# dangerous, but multiples isn't re-used

#             top_level_nodes = singles 		# + multiples
#             ls = len(singles)
#             lm = 0				# len(multiples)
#             lt = ls + lm
#             print 'scale_top_level_list: reset ls =', ls, ' lm =', lm, ' lt =', lt
#             sys.stdout.flush()

    if lt == 1:
        if lm == 1:

            # Special case.  We have a single bound binary node.  Its
            # children are the components we want to transform.  Note
            # that, if the components are binaries (or multiples),
            # they must be stable, so it is always OK to move the
            # components to periastron.

            root = multiples[0]

            print "scale_top_level_list: binary node"
            #print '\nunscaled binary node:'
            #print_multiple(root)
            comp1 = root.child1
            comp2 = root.child2
            compress_binary_components(comp1, comp2, scale)
            #print '\nscaled binary node:'
            #print_multiple(root)

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
        compress_binary_components(comp1, comp2, scale)
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

    # Recompute the external field, compute the tidal error, and
    # absorb it into the top-level energy.  Optional code.
    # Alternatively, we can simply absorb the tidal error into the
    # dEmult correction returned for bookkeeping purposes.

    phi_correction = zero

    phi_in_field_of_stars_to_add = potential_energy_in_field(
        top_level_nodes, 
        field,
        G = gravity_constant
    )
    
    #print 'phi_external_final_before =', phi_in_field_of_stars_to_remove
    #print 'phi_external_final_after =', phi_in_field_of_stars_to_add
    dphi = phi_in_field_of_stars_to_add - phi_in_field_of_stars_to_remove
    #print 'dphi =', dphi

    # Correction code parallels that above, but we must repeat it
    # here, since we have to complete the rescaling before the
    # tidal correction can be computed and applied.

    if lt == 1:

        # Only one top-level node.  Don't apply any correction.
        # Instead, always include the tidal term in dEmult.
        
        phi_correction = dphi

    elif lt == 2:

        # Absorb dphi into the relative motion of the top-level nodes,
        # using kepler.  Alternatively, add dphi to dEmult.

        comp1 = top_level_nodes[0]
        comp2 = top_level_nodes[1]
        phi_correction = dphi

    else:

        # Absorb dphi into the relative motion of the top-level nodes,
        # simply by scaling their velocities.  Need to improve this.
        # TODO  Alternatively, add dphi to dEmult.

        #print lt, 'top-level nodes'; sys.stdout.flush()
        phi_correction = dphi

    # Finally, include the internal top-level energy in dEmult.

    total_energy_of_stars_to_add = top_level_nodes.kinetic_energy()
    total_energy_of_stars_to_add \
        += top_level_nodes.potential_energy(G=gravity_constant)
    
    #print 'final etot =', total_energy_of_stars_to_add

    return total_energy_of_stars_to_add, phi_correction
