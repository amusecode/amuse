#!/usr/bin/env python

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
zero = 0.0|units.none
one = 1.0|units.none

class Initial_state():		# modeled on the scatter3 class
    eccentricity = 0		# no units in init -- add in make_triple
    gamma = 1.e-6
    m = 0.5
    M = 0.5
    planar = 0
    random_seed = -1
    impact_parameter = 0
    v_infinity = 0

#-----------------------------------------------------------------------
# Build functions could/should be offloaded to a separate module.

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

def set_inner_orbit(kep, init):
    semi = 1.0|nbody_system.length
    ecc = init.eccentricity
    kep.set_longitudinal_unit_vector(one, zero, zero)
    kep.set_normal_unit_vector(zero, zero, one)
    mean_an = 2*math.pi*numpy.random.random()
    kep.initialize_from_elements(1.0|nbody_system.mass, semi, ecc, mean_an)
    print 'inner semi, ecc =', semi.number, ecc

def set_outer_orbit(kep, init):

    mtotal = 1 + init.M

    # Multiply the incoming velocity by the critical value (see Hut &
    # Bahcall 1983).

    v_inf = init.v_infinity \
	        * math.sqrt( (1 - init.m) * init.m * mtotal / init.M )
    energy3 = 0.5 * v_inf * v_inf
    print 'm1, m2, m3 =', 1-init.m, init.m, init.M
    print 'energy3 =', energy3
    if energy3 > 0:
        semi = -0.5*mtotal/energy3|nbody_system.length
        ang_mom3 = init.impact_parameter * v_inf
        ecc = math.sqrt( 1 + 2 * energy3 * (ang_mom3/mtotal)**2)
        periastron = 0|nbody_system.length	# irrelevant
    else:
        semi = 0|nbody_system.length		# irrelevant
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

    print 'outer semi, ecc =', semi.number, ecc
    kep.initialize_from_elements(mtotal|nbody_system.mass, semi, ecc,
                                 mean_anomaly, time, periastron)
    print 'outer normal =', kep.get_normal_unit_vector()
    print 'outer periastron =', kep.get_periastron().number

def make_triple(init):

    # Create IDs and initial masses, positions, and velocities.
    # Convention: initial binary is (1,2).

    global kep

    # Inner orbit (1,2).

    set_inner_orbit(kep, init)
    rel_pos = kep.get_separation_vector()
    rel_vel = kep.get_velocity_vector()
    f = init.m
    pos1 = [-f*rel_pos[0], -f*rel_pos[1], -f*rel_pos[2]]
    vel1 = [-f*rel_vel[0], -f*rel_vel[1], -f*rel_vel[2]]
    pos2 = [(1-f)*rel_pos[0], (1-f)*rel_pos[1], (1-f)*rel_pos[2]]
    vel2 = [(1-f)*rel_vel[0], (1-f)*rel_vel[1], (1-f)*rel_vel[2]]

    # Outer orbit ((1,2),3).

    set_outer_orbit(kep, init)
    kep.return_to_radius((init.gamma/init.M)**(-1./3)|nbody_system.length)
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

    stars = datamodel.Particles(3)
    stars.id = id
    stars.mass = mass
    stars.position = pos
    stars.velocity = vel
    stars.radius = 0. | nbody_system.length
    return time,stars

#-----------------------------------------------------------------------

def get_binary_elements(p):
    global kep
    comp1 = p.child1
    comp2 = p.child2
    m = comp1.mass + comp2.mass
    kep.initialize_from_dyn(m,
		comp2.x-comp1.x, comp2.y-comp1.y, comp2.z-comp1.z,
		comp2.vx-comp1.vx, comp2.vy-comp1.vy, comp2.vz-comp1.vz)
    a,e = kep.get_elements()
    return m,a,e

def scatter3(init,
             accuracy_parameter = 0.1,
             delta_t = 0 | nbody_system.time,
             t_end = 1.e4 | nbody_system.time):

    t0 = clock()

    gravity = SmallN(redirection = "none")
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    gravity.parameters.timestep_parameter = accuracy_parameter

    t1 = clock()
    dtime1 = t1 - t0

    time,stars = make_triple(init)	# time = 0 at outer periastron
    number_of_stars = len(stars)

    #print stars

    print "adding particles"
    sys.stdout.flush()
    gravity.set_time(time);
    gravity.particles.add_particles(stars)
    print "committing particles"
    sys.stdout.flush()
    gravity.commit_particles()

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    # Don't have a proper unperturbed termination criterion in smallN.
    # For now, insist that the final time exceed minus the initial
    # time (recall that outer peri is at t = 0).

    t_crit = -time
    if delta_t.number <= 0.0: delta_t = 2*t_crit	# for efficiency

    print "evolving triple to completion in steps of", delta_t.number
    sys.stdout.flush()

    t2 = clock()
    dtime2 = t2 - t1

    dtime3 = 0.0
    dtime4 = 0.0

    while time < t_end:

        t3 = clock()

        time += delta_t
        gravity.evolve_model(time)

        t4 = clock()
        dtime3 += t4 - t3

        energy = gravity.get_kinetic_energy()+gravity.get_potential_energy()
        print "time =", time.number, "energy =", energy.number
        over = gravity.is_over()
        if time > t_crit and over:
            print '\ninteraction is over'
            gravity.update_particle_tree()
            gravity.update_particle_set()
            gravity.particles.synchronize_to(stars)
            channel.copy()
            channel.copy_attribute("index_in_code", "id")
            print "binaries:"
            x = trees.BinaryTreesOnAParticleSet(stars, "child1", "child2")
            roots = list(x.iter_roots())
            for r in roots:
                for level, particle in r.iter_levels():
                    print '  '*level, int(particle.id),
                    if not particle.child1 is None:
                        M,a,e = get_binary_elements(particle)
                        print ' ( M =', M.number, ' a =', a.number, \
                              ' e =', e, ')'
                    else:
                        print ''
            break
    
        sys.stdout.flush()
        dtime4 += clock() - t4

    if not over:
        print '\ninteraction is not over'
    gravity.stop()

    global time1, time2, time3, time4
    time1 += dtime1
    time2 += dtime2
    time3 += dtime3
    time4 += dtime4

if __name__ == '__main__':

    nscatter = 1
    accuracy_parameter = 0.1
    delta_t = 0.0 | nbody_system.time
    t_end = 1.e5 | nbody_system.time
    init = Initial_state();

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
            init.gamma = float(a)
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
            init.random_seed = int(a)
        elif o == "-t":
            t_end = float(a)
        elif o == "-v":
            init.v_infinity = float(a)
            if init.v_infinity < 0:
                raise Exception("v >= 0 required")
        else:
            print "unexpected argument", o

    # ** Note potential conflict between C++ and Python random generators. **

    if init.random_seed <= 0:
        numpy.random.seed()
        init.random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(init.random_seed)
    print "random seed =", init.random_seed, numpy.random.random()

    assert is_mpd_running()

    kep = Kepler(redirection = "none") #, debugger="gdb")
    kep.initialize_code()

    time1 = 0.0		# instantiation
    time2 = 0.0		# initialization
    time3 = 0.0		# integration
    time4 = 0.0		# investigation

    for i in range(nscatter):
        scatter3(init, accuracy_parameter, delta_t, t_end)
        if nscatter > 1: print '\n--------------------\n'

    print 'timing:  inst', time1, ' init', time2, \
	  ' evol', time3, ' over', time4

    kep.stop()
