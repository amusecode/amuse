import sys, unittest, numpy, random, collections, getopt, os, math

from amuse.units import nbody_system
from amuse.units import units
from amuse.community.hacs64.interface import Hacs64 as grav
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.datamodel import trees
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
def is_a_parent(child1_key, child2_key):
    return child1_key > 0 or child2_key > 0

def is_not_a_child(is_a_child):
    return is_a_child == 0

def print_log(s, gravity, E0 = 0.0 | nbody_system.energy):
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    try:
        Ebin = gravity.get_binary_energy()
    except:
        Ebin = 0 | nbody_system.energy
    Etop = T + U
    E = Etop + Ebin
    if E0 == 0 | nbody_system.energy: E0 = E
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

def get_component_binary_elements(comp1, comp2):
    kep = Kepler(redirection = "none")
    kep.initialize_code()

    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()	# per unit reduced mass, note
    kep.stop()

    return mass,a,e,r,E

def get_cm_binary_elements(p):
    return get_component_binary_elements(p.child1, p.child2)

def run_smallN(
        particles,
        end_time = 1000 | nbody_system.time,
        delta_t = 10 | nbody_system.time,
        accuracy_parameter = 0.1
    ):

    gravity = SmallN(redirection = "none") # , debugger="gdb")
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.cm_index = 2001
    gravity.commit_parameters()

    time = 0 | nbody_system.time

    print "\nadding particles to smallN"
    sys.stdout.flush()
    gravity.set_time(time);
    gravity.particles.add_particles(particles)
    print "committing particles to smallN"
    gravity.commit_particles()

    print "smallN: number_of_stars =", len(particles)
    print "smallN: evolving to time =", end_time.number, 
    print "in steps of", delta_t.number
    sys.stdout.flush()
    
    E0 = print_log('smallN', gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(particles)

    while time < end_time:
        time += delta_t
        print 'evolving smallN to time', time.number
        sys.stdout.flush()
        gravity.evolve_model(time)
        print_log('smallN', gravity, E0)
        over = gravity.is_over()
        if over.number:
            print 'interaction is over\n'; sys.stdout.flush()

            # Create a tree in the module representing the binary structure.

            gravity.update_particle_tree()

            # Return the tree structure to AMUSE.  Children are
            # identified by get_children_of_particle in interface.??,
            # and the information is returned in the copy operation.

            gravity.update_particle_set()
            gravity.particles.synchronize_to(particles)
            channel.copy()
            channel.copy_attribute("index_in_code", "id")

            gravity.stop()

            # Basic diagnostics: BinaryTreesOnAParticleSet creates
            # binary tree structure for all particles in the set; then
            # we loop over roots (top-level nodes) and print data on
            # all binaries below each.

            print "smallN binaries:"; sys.stdout.flush()
            x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
            roots = list(x.iter_roots())
            for r in roots:
                for level, particle in r.iter_levels():
                    print '  '*level, int(particle.id.number),
                    if not particle.child1 is None:
                        M,a,e,r,E = get_cm_binary_elements(particle)
                        print " mass = %.5e" % (M.number)
                        m1 = particle.child1.mass
                        m2 = particle.child2.mass
                        print_elements('      ', a, e, r, E*m1*m2/M)
                    else:
                        print ''
                    sys.stdout.flush()

            return E0
    
        sys.stdout.flush()
    
    gravity.stop()
    raise Exception("Did not finish the small-N simulation "
		    +"before end time {0}".format(end_time))
            
root_index = 1000
def new_root_index():
    global root_index
    root_index += 1
    return root_index
    
def openup_tree(star, stars, particles_in_encounter):

    # Create a binary tree for star.  Note that star is in the memory
    # set, and so star and leaves all have IDs reflecting the naming in
    # the dynamics module.

    root = trees.BinaryTreeOnParticle(star, 'child1', 'child2')

    # List the leaves.

    leaves = root.get_descendants_subset()

    # Compare with the position stored when replacing the particles
    # with the root particle, and move the particles accordingly.
    # Note that once the CM is in the gravity module, the components
    # are frozen and the coordinates are absolute, so we need the
    # original coordinates to offset them later.

    # Better just to store relative coordinates.  TODO.

    dx = star.x - star.original_x
    dy = star.y - star.original_y
    dz = star.z - star.original_z
    dvx = star.vx - star.original_vx
    dvy = star.vy - star.original_vy
    dvz = star.vz - star.original_vz
    
    leaves.x += dx
    leaves.y += dy
    leaves.z += dz
    leaves.vx += dvx
    leaves.vy += dvy
    leaves.vz += dvz
    
    particles_in_encounter.add_particles(leaves)
    
    # Remove the inner (CM) binary tree nodes from the stars list.
    # New ones will be added later as necessary.

    stars.remove_particles(root.get_inner_nodes_subset())

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

    top_level = stars.select(is_not_a_child, ["is_a_child"])

    min_dr = 1.e10
    id1 = star1.id.number
    id2 = star2.id.number
    for t in top_level:
        tid = t.id.number
        if tid != id1 and tid != id2:
            dr2 = sep2(t, star1)
            if dr2 > 0 and dr2 < min_dr:
                min_dr = dr2
                nnn = t
    min_dr = math.sqrt(min_dr)
    print 'star =', int(id1), ' min_dr =', min_dr, \
          ' nnn =', int(nnn.id.number), '(', nnn.mass.number, ')'
    print '    phi_tidal =', phi_tidal(star1, star2, nnn)
    print '    nnn pos:', nnn.x.number, nnn.y.number, nnn.z.number
    sys.stdout.flush()

    return nnn

def total_energy(slist):

    # Return the total energy of the particles in slist, in the center
    # of mass approximation.

    kinetic = 0.0
    potential = 0.0
    for s in slist:
        # print 'slist:', int(s.id.number)
        smass = s.mass.number
        spos = s.position.number
        svel = s.velocity.number
        kinetic  += smass*((svel**2).sum())
        spot = 0.0
        for ss in slist:
            if not ss == s:
                ssmass = ss.mass.number
                sspos = ss.position.number
                r = math.sqrt(((sspos-spos)**2).sum())
                # print int(s.id.number), int(ss.id.number), r
                spot -= ssmass/r
        potential += smass*spot
    kinetic /= 2
    potential /= 2
    # print 'kin, pot:', kinetic, potential
    return kinetic+potential

def relative_potential(slist, klist, stars):

    # Return the total potential energy of slist relative to stars
    # (excluding klist).

    potential = 0.0
    top_level = stars.select(is_not_a_child, ["is_a_child"])
    for t in top_level:
        tpot = 0.0
        skip = 0
        for s in klist:
            if s == t:
                skip = 1
                break
        if skip == 0:
            tpos = t.position.number
            for s in slist:
                dpos = s.position.number - tpos
                tpot -= s.mass.number/ math.sqrt((dpos*dpos).sum())
        potential += t.mass.number*tpot
    return potential

def offset_particle_tree(particle, dpos, dvel):

    # Recursively offset a particle and all of its descendants by
    # the specified position and velocity.

    if not particle.child1 is None:
        offset_particle_tree(particle.child1, dpos, dvel)
    if not particle.child2 is None:
        offset_particle_tree(particle.child2, dpos, dvel)
    particle.position += dpos
    particle.velocity += dvel
    # print 'offset', int(particle.id.number), 'by', dpos; sys.stdout.flush()

def compress_binary_components(comp1, comp2, scale):

    # Compress the two-body system consisting of comp1 and comp2 to
    # lie within distance scale of one another.

    pos1 = comp1.position
    pos2 = comp2.position
    sep12 = ((pos2-pos1)**2).sum()

    if sep12 > scale*scale:
        print '\ncompressing components', int(comp1.id.number), \
              'and', int(comp2.id.number), 'to separation', scale.number
        sys.stdout.flush()
        mass1 = comp1.mass
        mass2 = comp2.mass
        total_mass = mass1 + mass2
        vel1 = comp1.velocity
        vel2 = comp2.velocity
        cmpos = (mass1*pos1+mass2*pos2)/total_mass
        cmvel = (mass1*vel1+mass2*vel2)/total_mass

        # For now, create and delete a temporary kepler
        # process to handle the transformation.  Obviously
        # more efficient to define a single kepler at the
        # start of the calculation and reuse it.

        kep = Kepler(redirection = "none")
        kep.initialize_code()
        mass = comp1.mass + comp2.mass
        rel_pos = pos2 - pos1
        rel_vel = vel2 - vel1
        kep.initialize_from_dyn(mass,
                                rel_pos[0], rel_pos[1], rel_pos[2],
                                rel_vel[0], rel_vel[1], rel_vel[2])
        M,th = kep.get_angles()
        a,e = kep.get_elements()
        if e < 1:
            peri = a*(1-e)
            apo = a*(1+e)
        else:
            peri = a*(e-1)
            apo = 2*a		# OK - used ony to reset scale
        limit = peri + 0.01*(apo-peri)
        if scale < limit: scale = limit

        if M < 0:
            # print 'approaching'
            kep.advance_to_periastron()
            kep.advance_to_radius(limit)
        else:
            # print 'receding'
            if kep.get_separation() < scale:
                kep.advance_to_radius(limit)
            else:
                kep.return_to_radius(scale)

        # a,e = kep.get_elements()
        # r = kep.get_separation()
        # E,J = kep.get_integrals()
        # print 'kepler: a,e,r =', a.number, e.number, r.number
        # print 'E, J =', E, J

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

def print_elements(s, a, e, r, E):
    print '%s elements  a = %.4e  e = %.5f  r = %.4e  E = %.4e' \
	  % (s, a.number, e.number, r.number, E.number)

def print_multiple(m, level=0):

    # Recursively print the structure of (multipe) node m.

    print '    '*level, int(m.id.number), ' mass =', m.mass.number
    print '    '*level, 'pos =', m.position.number
    print '    '*level, 'vel =', m.velocity.number
    if not m.child1 is None and not m.child2 is None:
        M,a,e,r,E = get_component_binary_elements(m.child1, m.child2)
        print_elements(' '+'    '*level+'binary', a, e, r,
                       (E*m.child1.mass*m.child2.mass/M))
    if not m.child1 is None:
        print_multiple(m.child1, level+1)
    if not m.child2 is None:
        print_multiple(m.child2, level+1)

def print_pair_of_stars(s, star1, star2):
    m1 = star1.mass
    m2 = star2.mass
    M,a,e,r,E = get_component_binary_elements(star1, star2)
    print_elements(s, a, e, r, E*m1*m2/(m1+m2))
    print_multiple(star1)
    print_multiple(star2)

def scale_top_level_list(binaries, scale,
                         stars, klist, phi_external_init):

    # The smallN particles were followed until their interaction could
    # be unambiguously classified as over.  They may now be very far
    # apart.  Input binaries is an object describing the final
    # hierarchical structure of the interacting particles in smallN.
    # It consists of a flat tree of binary trees.

    # Scale the positions and velocities of the top-level nodes to
    # bring them within a sphere of diameter scale, conserving energy
    # and angular momentum (if possible).  Also offset all children to
    # reflect changes at the top level -- TODO will change when
    # offsets are implemented...

    # We are currently ignoring any possibility of a physical
    # collision during the smallN encounter.  TODO

    # Logic: 1 node   - must be a binary, use kepler to reduce to scale
    #        2 nodes  - use kepler, reduce binary children too? TODO
    #        3+ nodes - shrink radii and rescale velocities to preserve
    #                   energy, but this doesn't preserve angular
    #                   momentum TODO - also reduce children? TODO

    # TODO

    # Figure out the tree structure.

    singles = binaries.particles_not_in_a_multiple()
    multiples = binaries.roots()
    top_level_nodes = singles + multiples

    ls = len(singles)
    lm = len(multiples)
    lt = ls + lm

    if lt == 1:
        if lm == 1:

            # Special case.  We have a single bound binary node.  Its
            # children are the components we want to transform.  Note
            # that, if the components are binaries (or multiples),
            # they must be stable, so it is always OK to move the
            # components to periastron.

            root = multiples[0]
            print '\nunscaled binary node:'
            print_multiple(root)
            comp1 = root.child1
            comp2 = root.child2
            compress_binary_components(comp1, comp2, scale)
            print '\nscaled binary node:'
            print_multiple(root)

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
        print '\nunscaled top-level pair:'
        print_pair_of_stars('pair', comp1, comp2)
        compress_binary_components(comp1, comp2, scale)
        print '\nscaled top-level pair:'
        print_pair_of_stars('pair', comp1, comp2)

    else:

        # Now we have three or more top-level nodes.  We don't know
        # how to compress them in a conservative way.  For now, we
        # will conserve energy and think later about how to preserve
        # angular momentum.  TODO

        print lt, 'top-level nodes'; sys.stdout.flush()

    # Recompute the external field, compute the tidal error, and
    # absorb it into the top-level energy.  Optional code.
    # Alternatively, we can simply absorb the tidal error into the
    # dEmult correction returned for bookkeeping purposes.

    dEmult = 0.0

    phi_external_final = relative_potential(top_level_nodes, klist, stars)
    print 'phi_external_final =', phi_external_final
    dphi = phi_external_final - phi_external_init
    print 'dphi =', dphi

    # Correction code parallels that above, but we must repeat it
    # here, since we have to complete the rescaling before the
    # tidal correction can be computed and applied.

    if lt == 1:

        # Only one top-level node.  Don't apply any correction.
        # Instead, always include the tidal term in dEmult.

        dEmult += dphi

    elif lt == 2:

        # Absorb dphi into the relative motion of the top-level nodes,
        # using kepler.  Alternatively, add dphi to dEmult.

        comp1 = top_level_nodes[0]
        comp2 = top_level_nodes[1]
        dEmult += dphi

    else:

        # Absorb dphi into the relative motion of the top-level nodes,
        # simply by scaling their velocities.  Need to improve this.
        # TODO  Alternatively, add dphi to dEmult.

        print lt, 'top-level nodes'; sys.stdout.flush()
        dEmult += dphi

    # Finally, include the internal top-level energy in dEmult.

    etot = total_energy(top_level_nodes)
    print 'final etot =', etot
    dEmult += etot

    return dEmult

def manage_encounter(star1, star2, stars, gravity_stars):

    # Manage an encounter between star1 and star2.  stars is the
    # python memory data set.  gravity_stars points to the gravity
    # module data.  Return value is the energy correction due to
    # multiples.

    # Establish child flags for use during the encounter calculation.

    stars.is_a_child = 0|units.none
    parents = stars.select(is_a_parent, ["child1", "child2"])
    for s in stars:
        for p in parents:
            if p.child1 == s or p.child2 == s:
                s.is_a_child = 1|units.none

    # Currently we don't include neighbors in the integration and the
    # size limitation on any final multiple is poorly implemented.

    print '\nin manage_encounter'; sys.stdout.flush()

    # 1. Star1 and star2 reflect the data in the gravity module.  Find
    #    the corresponding particles in the local particle set.

    star1_in_memory = star1.as_particle_in_set(stars)	# pointers
    star2_in_memory = star2.as_particle_in_set(stars)
    
    # 1a. Copy the current position and velocity to mememory (need to
    #     create a better call for this, for example:
    #     star1.copy_to(star1_in_memory)

    star1_in_memory.position = star1.position
    star1_in_memory.velocity = star1.velocity
    star2_in_memory.position = star2.position
    star2_in_memory.velocity = star2.velocity
    
    # Basic diagnostics:

    print_pair_of_stars('**encounter**', star1_in_memory, star2_in_memory)

    print ''
    find_nnn(star1_in_memory, star2_in_memory, stars)
    find_nnn(star2_in_memory, star1_in_memory, stars)

    # 1b. Calculate the potential of [star1, star2] relative to the
    #     other top-level objects in the stars list (later: just use
    #     neighbors TODO).  Start the correction of dEmult by removing
    #     the initial top-level energy of the interacting particles
    #     from it.  (Add in the final energy later, on return from
    #     scale_top_level_list.)

    slist = [star1_in_memory, star2_in_memory]
    etot = total_energy(slist)
    print 'initial etot =', etot
    dEmult = -etot
    klist = [star1_in_memory, star2_in_memory]
    phi_external_init = relative_potential(slist, klist, stars)
    print 'phi_external_init =', phi_external_init

    # 2. Create a particle set to perform the close encounter
    #    calculation.

    particles_in_encounter = datamodel.Particles(0)
    
    # 3. Add stars to the encounter set, add in components when we
    #    encounter a binary.

    if not star1_in_memory.child1 is None:
        openup_tree(star1_in_memory, stars, particles_in_encounter)
    else:
        particles_in_encounter.add_particle(star1_in_memory)
        
    if not star2_in_memory.child1 is None:
        openup_tree(star2_in_memory, stars, particles_in_encounter)
    else:
        particles_in_encounter.add_particle(star2_in_memory)

    # particles_in_encounter.id = -1 | units.none  # need to make this -1 to
                                                   # ensure smallN will set the
                                                   # IDs, or else smallN seems
						   # to fail... TODO

    # *** Desirable to propagate the IDs into smallN, for internal
    # *** diagnostic purposes...

    # print 'particles in smallN encounter:'
    # print particles_in_encounter.to_string(['x','y','z',
    #                                         'vx','vy','vz','mass','id'])

    print '\nparticles in encounter (flat tree):'
    for p in particles_in_encounter:
        print_multiple(p)
    
    initial_scale \
        = math.sqrt(((star1.position
                      -star2.position).number**2).sum())|nbody_system.length
    print 'initial_scale =', initial_scale.number; sys.stdout.flush()

    # 4. Run the small-N encounter in the center of mass frame.

    # Probably desirable to make this a true scattering experiment by
    # separating star1 and star2 to a larger distance before starting
    # smallN.  TODO

    total_mass = star1.mass+star2.mass
    cmpos = (star1.mass*star1.position+star2.mass*star2.position)/total_mass
    cmvel = (star1.mass*star1.velocity+star2.mass*star2.velocity)/total_mass
    print 'CM KE =', 0.5*total_mass.number*((cmvel.number)**2).sum()
    for p in particles_in_encounter:
        p.position -= cmpos
        p.velocity -= cmvel

    run_smallN(particles_in_encounter)

    for p in particles_in_encounter:
        p.position += cmpos
        p.velocity += cmvel
    
    # print 'after smallN:'
    # print particles_in_encounter.to_string(['x','y','z',
    #                                         'vx','vy','vz',
    #                                         'mass', 'id'])
    # sys.stdout.flush()

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
    #      of the interacting subset into dEmult, so E(hacs64) + dEmult
    #      should be conserved.

    dEmult += scale_top_level_list(binaries, initial_scale,
                         stars, klist, phi_external_init)

    # 5c. Update the positions and velocities of the stars (leaves) in
    #     the encounter; copy the position and velocity attributes of
    #     all stars updated during the encounter to the stars particle
    #     set in memory.  Do not copy child or any other attributes.

    tmp_channel = particles_in_encounter.new_channel_to(stars)
    tmp_channel.copy_attributes(['x','y','z', 'vx', 'vy', 'vz'])

    # 5d. Add stars not in a binary to the gravity code.

    stars_not_in_a_multiple = binaries.particles_not_in_a_multiple()
    stars_not_in_a_multiple_in_stars = \
        stars_not_in_a_multiple.get_intersecting_subset_in(stars)
    if len(stars_not_in_a_multiple_in_stars) > 0:
        gravity_stars.add_particles(stars_not_in_a_multiple_in_stars)
        
    # 5e. Add the inner (CM) nodes (root plus subnodes) to the stars
    #     in memory (the descendant nodes are already part of the
    #     set).

    for root in binaries.iter_roots():
        stars_in_a_multiple = root.get_descendants_subset()
        # print 'root.get_inner_nodes_subset():'
        # print root.get_inner_nodes_subset(); sys.stdout.flush()
        stars.add_particles(root.get_inner_nodes_subset())
        
    # 5f. Add roots of the binaries tree to the gravity code.

    # Must set radii to reflect multiple structure.  TODO

    for root in binaries.iter_roots():
        root_in_stars = root.particle.as_particle_in_set(stars)
        root_in_stars.id = new_root_index()
        # print 'root_in_stars:'
        # print root_in_stars; sys.stdout.flush(); sys.stdout.flush()
        gravity_stars.add_particle(root_in_stars)
        
    # 5g. Store the original position and velocity of the root so that
    #     the subparticle position can be updated later.

    for root in binaries.iter_roots():        
        root_in_stars = root.particle.as_particle_in_set(stars)
        
        # Save root position and velocity so we can update the
        # position and velocity of the components when we open up the
        # binary tree.

        root_in_stars.original_x = root_in_stars.x
        root_in_stars.original_y = root_in_stars.y
        root_in_stars.original_z = root_in_stars.z
        root_in_stars.original_vx = root_in_stars.vx
        root_in_stars.original_vy = root_in_stars.vy
        root_in_stars.original_vz = root_in_stars.vz

    return dEmult

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

def xtest_multiples(infile = None,
                   number_of_stars = 64,
                   nmax = 2048,
                   end_time = 0.1   | nbody_system.time,
                   delta_t = 0.125 | nbody_system.time,
                   dt_max  = 0.0625 | nbody_system.time,
                   n_ngb   = 16,
                   eta_irr = 0.6,
                   eta_reg = 0.1,
                   softening_length = 0.0 | nbody_system.length,
                   random_seed = 1234):

    if infile != None: print "input file =", infile
    print "end_time =",  end_time.number
    print "nstars= ",    number_of_stars,
    print "nmax= ",      nmax,
    print "delta_t= ",   delta_t.number
    print "dt_max= ",    dt_max.number
    print "n_ngb= ",     n_ngb,
    print "eta_irr= ",   eta_irr
    print "eta_reg= ",   eta_reg
    print "eps2=    ",   softening_length.number**2
    print "\ninitializing the gravity module"
    sys.stdout.flush()

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed
    
    # Note that there are actually three GPU options to test:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    gravity = grav(number_of_workers = 1, debugger="none", redirection = "none")
    gravity.initialize_code()

    #-----------------------------------------------------------------

    if infile == None:

        print "making a Plummer model"
        stars = new_plummer_model(number_of_stars)

        id = numpy.arange(number_of_stars) 
        stars.id = id+1

        print "setting particle masses and radii"
	#stars.mass = (1.0 / number_of_stars) | nbody_system.mass
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass
        stars.radius = 0.0 | nbody_system.length

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared
                                    = 0 | nbody_system.length ** 2)

        time = 0.0 | nbody_system.time
        sys.stdout.flush()

    else:

        # Read the input data.  Units are dynamical.

        print "reading file", infile; sys.stdout.flush()

        id = []
        mass = []
        pos = []
        vel = []

        f = open(infile, 'r')
        count = 0
        for line in f:
            if len(line) > 0:
                count += 1
                cols = line.split()
                if count == 1: snap = int(cols[0])
                elif count == 2: number_of_stars = int(cols[0])
                elif count == 3: time = float(cols[0]) | nbody_system.time
                else:
                    if len(cols) >= 8:
                #        id.append(int(cols[0]))
                        id.append(count)
                        mass.append(float(cols[1]))
                        pos.append((float(cols[2]),
                                    float(cols[3]), float(cols[4])))
                        vel.append((float(cols[5]),
                                    float(cols[6]), float(cols[7])))
        f.close()

        stars = datamodel.Particles(number_of_stars)
        stars.id = id
        stars.mass = mass | nbody_system.mass
        stars.position = pos | nbody_system.length
        stars.velocity = vel | nbody_system.speed
        stars.radius = 0. | nbody_system.length
        nmax = 2*len(mass)

    # print "IDs:", stars.id.number
    sys.stdout.flush()

    global root_index
    root_index = len(stars) + 10000
    #-----------------------------------------------------------------
    
    gravity.parameters.nmax    = nmax;
    gravity.parameters.dtmax   = dt_max;
#    gravity.parameters.n_ngb   = n_ngb;
    gravity.parameters.eta_irr = eta_irr;
    gravity.parameters.eta_reg = eta_reg;
    gravity.parameters.eps2    = softening_length**2
    gravity.commit_parameters();

    print "adding particles"
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print ''
    print "number_of_stars =", number_of_stars
    print "evolving to time =", end_time.number, \
          "in steps of", delta_t.number
    sys.stdout.flush()

    E0 = print_log('hacs64', gravity)
    dEmult = 0.0
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)
    
    stopping_condition = gravity.stopping_conditions.collision_detection
    #stopping_condition.enable()
    
    # Tree structure on the stars dataset:

    stars.child1 = 0 | units.object_key
    stars.child2 = 0 | units.object_key
    
    while time < end_time:
        time += delta_t

        while gravity.get_time() < time:
            gravity.evolve_model(time)
            if stopping_condition.is_set():
                star1 = stopping_condition.particles(0)[0]
                star2 = stopping_condition.particles(1)[0]
                print '\n--------------------------------------------------'
                print 'stopping condition set at time', \
                    gravity.get_time().number

                E = print_log('hacs64', gravity, E0)
                print 'dEmult =', dEmult, 'dE =', (E-E0).number-dEmult
                # channel.copy()	# need other stars to be current in memory
                # print_energies(stars)

                # Synchronize everything for now.  Later we will just
                # synchronize neighbors if gravity supports that.  TODO
                gravity.synchronize_model()
                gravity.particles.synchronize_to(stars)

                dEmult += manage_encounter(star1, star2, stars,
                                           gravity.particles)

                # Recommit reinitializes all particles (and redundant
                # here, since done automatically).  Later we will just
                # recommit and reinitialize a list if gravity supports
                # it. TODO
                gravity.recommit_particles()

                E = print_log('hacs64', gravity, E0)
                print 'dEmult =', dEmult, 'dE =', (E-E0).number-dEmult
                print '\n--------------------------------------------------'
        
        ls = len(stars)
    
        # Copy values from the module to the set in memory.
        channel.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        if len(stars) != ls:
            if 0:
                print "stars:"
                for s in stars:
                    print " ", s.id.number, s.mass.number, \
			       s.x.number, s.y.number, s.z.number
            else:
                print "number of stars =", len(stars)
            sys.stdout.flush()

        E = print_log('hacs64', gravity, E0)
        print 'dEmult =', dEmult, 'dE =', (E-E0).number-dEmult

    print ''
    gravity.stop()

if __name__ == '__main__':

    infile = None
    
    N = 1024
    dt_max  = 0.0625 | nbody_system.time
    n_ngb = 16
    eta_irr = 0.8
    eta_reg = 0.14

#    eta_irr = 0.6
#    eta_reg = 0.1

    t_end = 1.0 | nbody_system.time
    delta_t = 0.125 | nbody_system.time
    
    random_seed = -1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGn:s:t:w:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-c":
            manage_encounters = int(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-e":
            softening_length = float(a) | nbody_system.length
        elif o == "-f":
            infile = a
        elif o == "-g":
            use_gpu = 0
        elif o == "-G":
            use_gpu = 0
            gpu_worker = 0
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print "unexpected argument", o

    
    Nmax = N*2
    softening_length = 0.0/N  | nbody_system.length

    assert is_mpd_running()
    test_multiples(infile,
                   N,
                   Nmax,
                   t_end, 
                   delta_t, 
                   dt_max,
                   n_ngb,
                   eta_irr,
                   eta_reg,
                   softening_length,
                   random_seed)
