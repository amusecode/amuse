#!/usr/bin/env python

# Scatter3: binary scattering package replicating most of the
# functionality of the scatter3 tool in Starlab.  See McMillan & Hut
# 1996, ApJ, 467, 348.

import sys, unittest, numpy, random, collections, getopt, os, math

from time import *

from amuse.units import nbody_system
from amuse.units import units
from amuse.community.newsmallN.interface import SmallN
from amuse.community.kepler.interface import Kepler

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.datamodel import trees
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

#-----------------------------------------------------------------------

class Initial_state():		# modeled on the Starlab scatter3 class
    m = 0.5			# no units in init -- add in make_triple
    M = 0.5
    eccentricity = 0.
    impact_parameter = 0.
    v_infinity = 1.0
    planar = 0

class Final_state():		# modeled on the Starlab scatter3 class
    is_over = 0
    time = 0.0
    mbinary = 0.0
    semimajoraxis = 0.0
    eccentricity = 0.0
    mescaper = 0.0
    escaper = 0
    separation = 0.0
    v_rel = 0.0

#-----------------------------------------------------------------------
# Build functions could/should be offloaded to a separate module.

zero = 0.0|units.none
one = 1.0|units.none

def normalized(a):
    amagi = 1./math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    b = [a[0]*amagi, a[1]*amagi, a[2]*amagi]
    return b

def sum3(a, b):
    c = [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    return c

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
    return c

def set_inner_orbit(init, kep):
    semi = 1.0|nbody_system.length
    ecc = init.eccentricity
    kep.set_longitudinal_unit_vector(one, zero, zero)
    kep.set_normal_unit_vector(zero, zero, one)
    mean_an = 2*math.pi*numpy.random.random()
    kep.initialize_from_elements(1.0|nbody_system.mass, semi, ecc, mean_an)
    print 'inner semi, ecc =', semi.number, ecc

def set_outer_orbit(init, kep):

    mtotal = 1 + init.M

    # Multiply the incoming velocity by the critical value (see Hut &
    # Bahcall 1983).

    v_inf = init.v_infinity \
	        * math.sqrt( (1 - init.m) * init.m * mtotal / init.M )
    energy3 = 0.5 * v_inf * v_inf
    print 'm1, m2, m3 =', 1-init.m, init.m, init.M
    print 'v_inf =', v_inf, 'energy3 =', energy3, \
          'rho =', init.impact_parameter
    if energy3 > 0:
        semi = -0.5*mtotal/energy3|nbody_system.length
        ang_mom3 = init.impact_parameter * v_inf
        ecc = math.sqrt( 1 + 2 * energy3 * (ang_mom3/mtotal)**2)
        periastron = -semi*max(ecc-1.0, 0.0)	# not used
    else:
        semi = 0|nbody_system.length		# not used
        ecc = 1
        periastron = init.impact_parameter|nbody_system.length

    # Orientation:

    if init.planar != 0:
        kep.set_longitudinal_unit_vector(one, zero, zero)
        if init.planar == 1:
            kep.set_normal_unit_vector(zero, zero, one)
        else:
            kep.set_normal_unit_vector(zero, zero, -one)
    else:
        costheta = 2*numpy.random.random() - 1
        sintheta = math.sqrt(max(0, 1-costheta**2))
        phi = 2*math.pi*numpy.random.random()
        long = [sintheta*math.cos(phi), sintheta*math.sin(phi), costheta]
        kep.set_longitudinal_unit_vector(long[0],
                                         long[1],
                                         long[2])
        if abs(long[0]) < 0.5:
            temp = [1, 0, 0]
        else:
            temp = [0, 1, 0]
        trans = normalized(cross(long, temp))
        normal = cross(long, trans)
        psi = 2*math.pi*numpy.random.random()
        cospsi = math.cos(psi)
        sinpsi = math.sin(psi)
        normal = [cospsi*trans[0]+sinpsi*normal[0],
                  cospsi*trans[1]+sinpsi*normal[1],
                  cospsi*trans[2]+sinpsi*normal[2]]
        kep.set_normal_unit_vector(normal[0],
                                   normal[1],
                                   normal[2])

    time = 0.0|nbody_system.time
    mean_anomaly = 0				# t = 0 at periastron
    if periastron.number == 0: mean_anomaly = -1.e-3
    print 'mean_anomaly =', mean_anomaly

    print 'outer semi, ecc =', semi.number, ecc
    kep.initialize_from_elements(mtotal|nbody_system.mass, semi, ecc,
                                 mean_anomaly, time, periastron)
    print 'outer normal =', kep.get_normal_unit_vector()
    print 'outer periastron =', kep.get_periastron().number
    #kep.print_all()

def make_triple(init, kep):

    # Create IDs and initial masses, positions, and velocities.
    # Convention: initial binary is (1,2).

    # Inner orbit (1,2).

    set_inner_orbit(init, kep)
    rel_pos = kep.get_separation_vector()
    rel_vel = kep.get_velocity_vector()
    f = init.m
    pos1 = [-f*rel_pos[0], -f*rel_pos[1], -f*rel_pos[2]]
    vel1 = [-f*rel_vel[0], -f*rel_vel[1], -f*rel_vel[2]]
    pos2 = [(1-f)*rel_pos[0], (1-f)*rel_pos[1], (1-f)*rel_pos[2]]
    vel2 = [(1-f)*rel_vel[0], (1-f)*rel_vel[1], (1-f)*rel_vel[2]]

    # Outer orbit ((1,2),3).

    set_outer_orbit(init, kep)
    # print '----------'
    # kep.print_all()
    gamma = gravity.parameters.unperturbed_threshold
    kep.return_to_radius((gamma/init.M)**(-1./3)|nbody_system.length)
    # print '----------'
    # kep.print_all()
    # print '----------'
    time = kep.get_time()

    rel_pos = kep.get_separation_vector()
    rel_vel = kep.get_velocity_vector()
    f = init.M/(1+init.M)
    pos12 = [-f*rel_pos[0], -f*rel_pos[1], -f*rel_pos[2]]
    vel12 = [-f*rel_vel[0], -f*rel_vel[1], -f*rel_vel[2]]

    print 'outer separation =', kep.get_separation().number, \
	  ' time =', kep.get_time().number

    pos1 = sum3(pos1, pos12)
    vel1 = sum3(vel1, vel12)
    pos2 = sum3(pos2, pos12)
    vel2 = sum3(vel2, vel12)
    pos3 = [(1-f)*rel_pos[0], (1-f)*rel_pos[1], (1-f)*rel_pos[2]]
    vel3 = [(1-f)*rel_vel[0], (1-f)*rel_vel[1], (1-f)*rel_vel[2]]

    print 'initial time =', time.number, "(time of outer periastron = 0)"

    # Create the 3-body system.

    id = [1, 2, 3]
    mass = [1-init.m, init.m, init.M] | nbody_system.mass
    pos = [pos1, pos2, pos3]
    vel = [vel1, vel2, vel3]

    return time, id, mass, pos, vel

def make_triple2(init, kep, gamma):

    # Use a kepler built-in to replace the Python functionality of
    # make_triple.

    # Create IDs and initial masses, positions, and velocities.
    # Convention: initial binary is (1,2).

    t,m1,m2,m3,x1,x2,x3,y1,y2,y3,z1,z2,z3, \
        vx1,vx2,vx3,vy1,vy2,vy3,vz1,vz2,vz3 \
	= kep.make_binary_scattering(init.m|nbody_system.mass,
                                     init.eccentricity,
                                     init.M|nbody_system.mass,
                                     init.v_infinity|nbody_system.speed,
                                     init.impact_parameter|nbody_system.length,
                                     gamma, init.planar)

    id = [1, 2, 3]
    mass = [m1,m2,m3]
    pos = [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
    vel = [[vx1,vy1,vz1], [vx2,vy2,vz2], [vx3,vy3,vz3]]

    return t, id, mass, pos, vel

#-----------------------------------------------------------------------

def get_binary_elements(p, kep):
    comp1 = p.child1
    comp2 = p.child2
    m = comp1.mass + comp2.mass
    kep.initialize_from_dyn(m,
		comp2.x-comp1.x, comp2.y-comp1.y, comp2.z-comp1.z,
		comp2.vx-comp1.vx, comp2.vy-comp1.vy, comp2.vz-comp1.vz)
    a,e = kep.get_elements()
    return m,a,e

def get_final_state(stars):

    final = Final_state()
    final.is_over = 1
    final.escaper = -1
    ionized = 1

    tree = stars.as_binary_tree()
    sep = numpy.zeros(3)
    vel = numpy.zeros(3)
    for node in tree.iter_children():
        if node.particle.child1 is None:
            final.escaper = node.particle.id
            final.mescaper = node.particle.mass.number
            sep += numpy.array([node.particle.x.number,
                                node.particle.y.number,
                                node.particle.z.number])
            vel += numpy.array([node.particle.vx.number,
                                node.particle.vy.number,
                                node.particle.vz.number])
        else:
            M,a,e = get_binary_elements(node.particle, kep)
            final.mbinary = M.number
            final.semimajoraxis = a.number
            final.eccentricity = e
            ionized = 0
            sep -= numpy.array([node.particle.x.number,
                                node.particle.y.number,
                                node.particle.z.number])
            vel -= numpy.array([node.particle.vx.number,
                                node.particle.vy.number,
                                node.particle.vz.number])

    if ionized == 0:
        final.separation = math.sqrt((sep*sep).sum())
        final.v_rel = math.sqrt((vel*vel).sum())
    else:
        final.mbinary = -1.0
        final.semimajoraxis = -1.0
        final.eccentricity = -1.0
        final.escaper = -1
        final.mescaper = -1.0
        final.separation = -1.0
        final.v_rel = -1.0

    return final

def scatter3(init, kep, gravity,
             delta_t = 0 | nbody_system.time,
             t_end = 1.e4 | nbody_system.time):

    t1 = clock()		# <----------------- t1 -----------------

    # Create the 3-body system.

    #time, id, mass, pos, vel = make_triple(init, kep)  # time = 0 at outer peri
    time, id, mass, pos, vel = \
        make_triple2(init, kep, gravity.parameters.unperturbed_threshold)
    
    stars = datamodel.Particles(3)
    stars.id = id
    stars.mass = mass
    stars.position = pos
    stars.velocity = vel
    stars.radius = 0. | nbody_system.length

    #print stars

    print "adding particles"
    sys.stdout.flush()
    gravity.set_time(time)
    gravity.particles.add_particles(stars)
    #print "committing particles"
    #gravity.commit_particles()
    sys.stdout.flush()

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    # Don't have a proper unperturbed termination criterion in smallN.
    # For now, insist that the final time exceed minus the initial
    # time (recall that outer peri is at t = 0).

    t_crit = -time
    if delta_t.number <= 0.0: delta_t = 2*t_crit	# for efficiency

    print "evolving triple to completion in steps of", delta_t.number
    sys.stdout.flush()

    t2 = clock()		# <----------------- t2 -----------------

    dt_init = t2 - t1
    dt_evolve = 0.0
    dt_over = 0.0
    dt_tree = 0.0

    final = Final_state()
    over = 0

    while time < t_end and over == 0:

        tt3 = clock()		# <----------------- tt3 ----------------

        time += delta_t
        gravity.evolve_model(time)

        energy = gravity.get_kinetic_energy()+gravity.get_potential_energy()
        print "time =", time.number, "energy =", energy.number

        tt4 = clock()		# <----------------- tt4 ----------------

        if time > t_crit:

            ttt5 = clock()	# <---------------- ttt5 ----------------

            over = gravity.is_over()

            ttt6 = clock()	# <---------------- ttt6 ----------------

            if over:
                #print '\nscatter3: interaction is over'

                gravity.update_particle_tree()
                gravity.update_particle_set()
                gravity.particles.synchronize_to(stars)
                channel.copy()
                channel.copy_attribute("index_in_code", "id")

		# Determine the final state.

                final = get_final_state(stars)
                final.time = time.number
 
            ttt7 = clock()	# <---------------- ttt7 ----------------

            dt_over += ttt6 - ttt5
            dt_tree += ttt7 - ttt6

        dt_evolve += tt4 - tt3
 
    if not over:
        #print '\nscatter3: interaction is not over'
        final.is_over = 0
        final.mbinary = -1.0
        final.semimajoraxis = -1.0
        final.eccentricity = -1.0
        final.escaper = -1
        final.mescaper = -1.0
        final.separation = -1.0
        final.v_rel = -1.0

    # Clean up internal data for recycling.

    #gravity.particles.remove_particles(gravity.particles)
    gravity.cleanup_code()
    gravity.initialize_code()

    return final,numpy.array([dt_init, dt_evolve, dt_over, dt_tree])

if __name__ == '__main__':

    nscatter = 1
    accuracy_parameter = 0.1
    gamma = 1.e-6
    delta_t = 0.0 | nbody_system.time
    t_end = 1.e5 | nbody_system.time
    random_seed = -1

    init = Initial_state()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "A:d:e:g:m:M:n:pPr:s:t:v:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    # Command-line arguments are modeled on those in starlab.
    # Units: G = 1, binary mass = 1, binary semi-major axis = 1.

    for o, a in opts:
        if o == "-A":
            accuracy_parameter = float(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-e":
            init.eccentricity = float(a)
            if init.eccentricity >= 1:
                raise Exception("e < 1 required")
        elif o == "-g":
            gamma = float(a)
        elif o == "-m":
            init.m = float(a)
            if init.m <= 0 or init.m >= 1:
                raise Exception("0 < m < 1 required")
        elif o == "-M":
            init.M = float(a)
            if init.M <= 0:
                raise Exception("M > 0 required")
        elif o == "-n":
            nscatter = int(a)
        elif o == "-p":
            init.planar = 1
        elif o == "-P":
            init.planar = -1
        elif o == "-r":
            init.impact_parameter = float(a)
            if init.impact_parameter < 0:
                raise Exception("r >= 0 required")
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a)|nbody_system.time
        elif o == "-v":
            init.v_infinity = float(a)
            if init.v_infinity < 0:
                raise Exception("v >= 0 required")
        else:
            print "unexpected argument", o

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print "random seed =", random_seed, numpy.random.random()

    #-----------------------------------------------------------------

    assert is_mpd_running()

    # Instantiate workers once only and pass to scatter3 as arguments.

    gravity = SmallN(redirection = "none")
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.unperturbed_threshold = gamma

    kep = Kepler(redirection = "none") #, debugger="gdb")
    kep.initialize_code()
    kep.set_random(random_seed)	  # ** Note potential conflict between C++
				  # ** and Python random generators.

    # Timing:

    cpu = numpy.zeros(4)

    for i in range(nscatter):

        final,dcpu = scatter3(init, kep, gravity, delta_t, t_end)
        cpu += dcpu

        print ''
        if final.is_over == 0:
            print 'interaction is not over'
        else:
            if final.escaper > 0:
                print 'escaper =', final.escaper, \
                      'sep =', final.separation, \
                      'vel =', final.v_rel
                print 'binary (mass =', final.mbinary, \
                      'semi =', final.semimajoraxis, \
                      'ecc =', final.eccentricity, ')'
            else:
                print 'ionization'

        if nscatter > 1: print '\n--------------------\n'

    print 'timing:  init', cpu[0], 'evol', cpu[1], 'over', cpu[2], \
          'tree', cpu[3]
    print ''

    gravity.stop()
    kep.stop()
