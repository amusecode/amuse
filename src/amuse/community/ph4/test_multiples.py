import sys
import unittest
import numpy
import random
import collections
import getopt
import os

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import core
from amuse.support.data import particle_attributes
from amuse.support.codes.core import is_mpd_running
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import new_salpeter_mass_distribution_nbody

from amuse.support.data import trees
from amuse.community.ph4.interface import ph4 as grav
from amuse.community.newsmallN.interface import SmallN
from amuse.community.kepler.interface import Kepler

def print_log(time, gravity, E0 = 0.0 | nbody_system.energy):
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
    print "time =", time.number, " energy = ", E.number, \
	" dE/E0 = ", (E/E0 - 1).number
    print '%s %.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f' % \
	("%%", time.number, M.number, T.number, U.number, \
         E.number, Ebin.number, Rv.number, Q.number)
    sys.stdout.flush()
    return E

def get_binary_elements(p):
    comp1 = p.child1
    comp2 = p.child2
    kep = Kepler(redirection = "none")
    kep.initialize_code()

    mass = comp1.mass + comp2.mass
    pos = [comp2.x-comp1.x, comp2.y-comp1.y, comp2.z-comp1.z]
    vel = [comp2.vx-comp1.vx, comp2.vy-comp1.vy, comp2.vz-comp1.vz]
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()

    return mass,a,e

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

    time = 0 | nbody_system.time

    print "adding particles to smallN"
    sys.stdout.flush()
    gravity.set_time(time);
    gravity.particles.add_particles(particles)
    print "committing particles"
    gravity.commit_particles()

    print ''
    print "smallN: number_of_stars =", len(particles)
    print "smallN: evolving to time =", end_time.number, 
    print " in steps of", delta_t.number
    sys.stdout.flush()
    
    E0 = print_log(time, gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(particles)

    while time < end_time:
        time += delta_t
        gravity.evolve_model(time)
        print_log(time, gravity, E0)

        print "smallN time =", gravity.get_time().number
        over = gravity.is_over()
        if over.number:
            print 'interaction is over\n'
            gravity.update_particle_tree()
            gravity.update_particle_set()
            gravity.particles.synchronize_to(particles)
            channel.copy()
            gravity.stop()
            channel.copy_attribute("index_in_code", "id")
            print "binaries:"
            x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
            roots = list(x.iter_roots())
            for r in roots:
                for level, particle in r.iter_levels():
                    print '  '*level, int(particle.id.number),
                    if not particle.child1 is None:
                        M,a,e = get_binary_elements(particle)
                        print ' ( M =', M.number, ' a =', a.number, \
                              ' e =', e.number, ')'
                    else:
                        print ''
                    sys.stdout.flush()
            return particles
    
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
    root = trees.BinaryTreeOnParticle(star, 'child1', 'child2')
    leaves = root.get_descendants_subset()
    
    # Compare with the position stored when replacing the particles
    # with the root particle, and move the particles accordingly.
    
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
    
    # Remove the binary tree.  It will be recreated later.

    stars.remove_particles(root.get_inner_nodes_subset())
    
def manage_encounter(star1, star2, stars, gravity_stars):
    
    # 1. Find corresponding particles in memory set.
    star1_in_memory = star1.as_particle_in_set(stars)
    star2_in_memory = star2.as_particle_in_set(stars)
    
    # 2. Create a set to perform the close encounter calculation.
    particles_in_encounter = core.Particles(0)
    
    # 3. Add stars to the close encounter, put substars in when
    #    encountering a binary.
    if not star1_in_memory.child1 is None:
        openup_tree(star1_in_memory, stars, particles_in_encounter)
    else:
        particles_in_encounter.add_particle(star1)
        
    if not star2_in_memory.child1 is None:
        openup_tree(star2_in_memory, stars, particles_in_encounter)
    else:
        particles_in_encounter.add_particle(star2)
    
    particles_in_encounter.id = -1 | units.none # need to make this -1 to
                                                # ensure smallN will set the
                                                # IDs, or else smallN seems
						# to fail
    print particles_in_encounter.to_string(['x','y','z','vx','vy','vz','mass'])
    
    # 4. Run the small-N encounter.
    run_smallN(particles_in_encounter)
    
    # 5. Bookkeeping after encounter.
    
    # 5.a Create object to handle the binary information.
    binaries = trees.BinaryTreesOnAParticleSet(particles_in_encounter,
                                               "child1", "child2")
    
    # 5.b. Remove encountered stars from gravity code.
    gravity_stars.remove_particle(star1)
    gravity_stars.remove_particle(star2)
    
    # 5.c. Update the position and velocity of the stars in the
    #      encounter.
    tmp_channel = particles_in_encounter.new_channel_to(stars)
    tmp_channel.copy_attributes(['x','y','z', 'vx', 'vy', 'vz'])
        
    # 5.d Add stars not in a binary to the gravity code.
    stars_not_in_a_multiple = binaries.particles_not_in_a_multiple()
    stars_not_in_a_multiple_in_stars = \
        stars_not_in_a_multiple.get_intersecting_subset_in(stars)
    if len(stars_not_in_a_multiple_in_stars) > 0:
        gravity_stars.add_particles(stars_not_in_a_multiple_in_stars)
        
    # 5.e. Add the inner nodes (root plus subnodes) to the stars in
    #      memory (the descendant nodes are already part of the set).
    for root in binaries.iter_roots():
        stars_in_a_multiple = root.get_descendants_subset()
        stars.add_particles(root.get_inner_nodes_subset())
        
    # 5.f. Add roots of the binaries tree to the gravity code.
    for root in binaries.iter_roots():
        root_in_stars = root.particle.as_particle_in_set(stars)
        root_in_stars.id = new_root_index() | units.none
        gravity_stars.add_particle(root_in_stars)
        
    # 5.g. Store the original position and velocity of the root so
    #      that the subparticle position can be updated later.
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

def test_multiples(infile = None, number_of_stars = 40,
             end_time = 10 | nbody_system.time,
             delta_t = 1 | nbody_system.time,
             n_workers = 1, use_gpu = 1, gpu_worker = 1,
             accuracy_parameter = 0.1,
             softening_length = -1 | nbody_system.length,
             manage_encounters = 1):

    if infile != None: print "input file =", infile
    print "end_time =", end_time.number
    print "delta_t =", delta_t.number
    print "n_workers =", n_workers
    print "use_gpu =", use_gpu
    print "manage_encounters =", manage_encounters
    print "\ninitializing the gravity module"
    sys.stdout.flush()

    # Note that there are actually three GPU options to test:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    if gpu_worker == 1:
        try:
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none", mode = "gpu")
        except Exception as ex:
            gravity = grav(number_of_workers = n_workers, redirection = "none")
    else:
	gravity = grav(number_of_workers = n_workers, redirection = "none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    #-----------------------------------------------------------------

    if infile == None:

        print "making a Plummer model"
        stars = MakePlummerModel(number_of_stars).result

        id = numpy.arange(number_of_stars) 
        stars.id = id+1 | units.none

        print "setting particle masses and radii"
	#stars.mass = (1.0 / number_of_stars) | nbody_system.mass
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass
        stars.radius = 0.0 | nbody_system.length

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared
                                    = gravity.parameters.epsilon_squared)

        time = 0.0 | nbody_system.time
        sys.stdout.flush()

    else:

        # Read the input data.  Units are dynamical.

        print "reading file", infile

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
                        id.append(int(cols[0]))
                        mass.append(float(cols[1]))
                        pos.append((float(cols[2]),
                                    float(cols[3]), float(cols[4])))
                        vel.append((float(cols[5]),
                                    float(cols[6]), float(cols[7])))
	f.close()

        stars = core.Particles(number_of_stars)
        stars.id = id | units.none
        stars.mass = mass | nbody_system.mass
        stars.position = pos | nbody_system.length
        stars.velocity = vel | nbody_system.speed
        stars.radius = 0. | nbody_system.length

    # print "IDs:", stars.id.number
    sys.stdout.flush()

    global root_index
    rooot_index = len(stars) + 1000
    #-----------------------------------------------------------------

    if softening_length == -1 | nbody_system.length:
        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
			| nbody_system.length**2
    else:
        eps2 = softening_length*softening_length

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu
    gravity.parameters.manage_encounters = manage_encounters

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

    E0 = print_log(time, gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)
    
    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    stars.child1 = 0 | units.object_key
    stars.child2 = 0 | units.object_key
    
    while time < end_time:
        time += delta_t
        gravity.evolve_model(time)

        if stopping_condition.is_set():
            print '\nstopping condition set at time', \
                gravity.get_time().number,'for:\n'
            star1 = stopping_condition.particles(0)[0]
            star2 = stopping_condition.particles(1)[0]
            
            print "index of star 1", star1.index_in_code.number
            print "index of star 2", star2.index_in_code.number
            
            manage_encounter(star1, star2, stars, gravity.particles)
            
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

        print_log(time, gravity, E0)
        sys.stdout.flush()

    print ''
    gravity.stop()

if __name__ == '__main__':

    infile = None
    N = 100
    t_end = 1000.0 | nbody_system.time
    delta_t = 10.0 | nbody_system.time
    n_workers = 1
    use_gpu = 0
    gpu_worker = 0
    accuracy_parameter = 0.1
    softening_length = 0  | nbody_system.length
    random_seed = -1
    manage_encounters = 4

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

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed

    assert is_mpd_running()
    test_multiples(infile, N, t_end, delta_t, n_workers,
             use_gpu, gpu_worker,
             accuracy_parameter, softening_length,
             manage_encounters)
