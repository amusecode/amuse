#!/bin/env python

import math
import collections
import getopt
import numpy
import os
import random
import sys
import unittest 
from time import clock as cputime
from time import time as wallclocktime

from amuse.community.ph4.interface import ph4 as grav
from amuse.units import nbody_system
from amuse.units import units

from amuse import datamodel
from amuse.datamodel import particle_attributes as pa
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

def print_log(pre, time, gravity, E0 = 0.0 | nbody_system.energy,
              cpu0 = 0.0, wall0 = 0.0):

    # Standard log output.

    cpu = cputime()
    wall = wallclocktime()
    N = len(gravity.particles)
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    Etop = T + U
    E = Etop
    if E0 == 0 | nbody_system.energy: E0 = E
    Rvir = -0.5*M*M/U
    Q = -T/U
    com = pa.center_of_mass(gravity.particles)
    comv = pa.center_of_mass_velocity(gravity.particles)
    if N >= 100:
        dcen,rcore,rhocore \
            = pa.densitycentre_coreradius_coredens(gravity.particles)
        cmx,cmy,cmz = dcen
        lagr,mf = pa.LagrangianRadii(gravity.particles, cm=dcen)  # no units!

    print('')
    print(pre+"time=", time.number)
    print(pre+"cpu=", cpu-cpu0)
    print(pre+"wall=", wall-wall0)
    print(pre+"Ntot=", N)
    print(pre+"mass=", M.number)
    print(pre+"Etot=", E.number)
    print(pre+"dE/E=", E/E0 - 1)
    print(pre+"Rvir=", Rvir.number)
    print(pre+"Qvir=", Q)
    cmx,cmy,cmz = com
    print(pre+"cmpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    cmx,cmy,cmz = comv
    print(pre+"cmvel[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number))
    if N >= 100:
        cmx,cmy,cmz = dcen
        print(pre+"dcpos[3]= %.8f %.8f %.8f" \
            		% (cmx.number, cmy.number, cmz.number))
        print(pre+"Rcore=", rcore.number)
        print(pre+"Mlagr[9]=", end=' ')
        for m in mf: print("%.4f" % (m), end=' ')
        print('')
        print(pre+"Rlagr[9]=", end=' ')
        for r in lagr.number: print("%.8f" % (r), end=' ')
        print('')

    sys.stdout.flush()
    return E,cpu,wall

def run_ph4(infile = None, number_of_stars = 40,
             end_time = 10 | nbody_system.time,
             delta_t = 1 | nbody_system.time,
             n_workers = 1, use_gpu = 1, gpu_worker = 1, gpu_id = -1,
             accuracy_parameter = 0.1,
             softening_length = -1 | nbody_system.length,
             manage_encounters = 1):

    if infile != None: print("input file =", infile)
    print("end_time =", end_time.number)
    print("delta_t =", delta_t.number)
    print("n_workers =", n_workers)
    print("use_gpu =", use_gpu)
    print("manage_encounters =", manage_encounters)
    print("initializing the gravity module")
    sys.stdout.flush()

    # Note that there are actually really three GPU options to test:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    #print "1"; sys.stdout.flush()
    
    gpu = 0
    if gpu_worker == 1:
        try:
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none", mode = "gpu")
            #              debugger='valgrind')
            gpu = 1
        except Exception as ex:
            print('*** GPU worker code not found. Reverting to non-GPU code. ***')
            gpu = 0

    if gpu == 0:
        gravity = grav(number_of_workers = n_workers,
                       redirection = "none")
        #              debugger='valgrind')

    #print "2"; sys.stdout.flush()
    gravity.initialize_code()

    #print "3"; sys.stdout.flush()
    gravity.parameters.set_defaults()

    gravity.parameters.gpu_id = gpu_id

    #-----------------------------------------------------------------

    #print "4"; sys.stdout.flush()

    print("making a Plummer model")
    stars = new_plummer_model(number_of_stars)

    id = numpy.arange(number_of_stars)
    stars.id = id+1

    print("setting particle masses and radii")
    stars.mass = (1.0 / number_of_stars) | nbody_system.mass
    if 0:
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
        stars.mass = scaled_mass
    stars.radius = 0.0 | nbody_system.length

    print("centering stars")
    stars.move_to_center()
    if 0:
        print("scaling stars to virial equilibrium")
        stars.scale_to_standard(smoothing_length_squared
                                = gravity.parameters.epsilon_squared)

    time = 0.0 | nbody_system.time
    sys.stdout.flush()

    #-----------------------------------------------------------------

    #print "5"; sys.stdout.flush()
    if softening_length == -1 | nbody_system.length:
        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
			| nbody_system.length**2
    else:
        eps2 = softening_length*softening_length

    #print "6"; sys.stdout.flush()
    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu
    gravity.parameters.manage_encounters = manage_encounters

    print("adding particles")
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print('')
    print("number_of_stars =", number_of_stars)
    sys.stdout.flush()

    E0,cpu0,wall0 = print_log('', time, gravity)
    
    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    #-----------------------------------------------------------------

    cpu0 = cputime()
    t0 = 0.
    pi = math.pi
    times = [1., 2., pi, 4*pi/3, 5., 2*pi, 2*pi + pi/100,
             2*pi + pi/5, 7., 8., 3*pi, 10.]

    gravity.parameters.force_sync = 1	# stays set until explicitly unset

    for t in times:
        time = t|nbody_system.time

        print("\nEvolving to time", time)
        sys.stdout.flush()

        gravity.parameters.block_steps = 0
        gravity.parameters.total_steps = 0
        gravity.evolve_model(time)

        dt = t - t0
        t0 = t
        cpu = cputime()
        dcpu = cpu - cpu0
        cpu0 = cpu

        # Ensure that the stars list is consistent with the internal
        # data in the module.

        ls = len(stars)

	# Update the bookkeeping: synchronize stars with the module data.

        try:
            gravity.update_particle_set()
            gravity.particles.synchronize_to(stars)
        except:
            pass
    
        # Copy values from the module to the set in memory.

        channel.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        if stopping_condition.is_set():
            star1 = stopping_condition.particles(0)[0]
            star2 = stopping_condition.particles(1)[0]
            print('\nstopping condition set at time', \
                gravity.get_time().number,'for:\n')
            print(star1)
            print('')
            print(star2)
            print('')
            raise Exception("no encounter handling")

        if len(stars) != ls:
            if 0:
                print("stars:")
                for s in stars:
                    print(" ", s.id.number, s.mass.number, \
			       s.x.number, s.y.number, s.z.number)
            else:
                print("number of stars =", len(stars))
            sys.stdout.flush()

        print_log('', time, gravity, E0, cpu0, wall0)

        print('@@@')
        print('@@@ t =', time.number, ' dt =', dt)
        print('@@@ sync_time =', gravity.parameters.sync_time.number)
        print('@@@ dcpu/dt =', dcpu/dt)
        nb = gravity.parameters.block_steps
        ns = gravity.parameters.total_steps
        print('@@@ d(block_steps) =', nb, ' #/dt =', nb/dt)
        print('@@@ d(total steps) =', ns, ' #/dt =', ns/dt)

        #print stars
        sys.stdout.flush()

    #-----------------------------------------------------------------

    print('')
    gravity.stop()

if __name__ == '__main__':

    infile = None
    N = 100
    t_end = 5.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 1
    use_gpu = 1
    gpu_worker = 1
    gpu_id = -1
    accuracy_parameter = 0.1
    softening_length = 0.01  | nbody_system.length
    random_seed = -1
    manage_encounters = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGi:n:s:t:w:")
    except getopt.GetoptError as err:
        print(str(err))
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
        elif o == "-i":
            gpu_id = int(a)
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print("unexpected argument", o)

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print("random seed =", random_seed)


    #os.system('env')


    assert is_mpd_running()
    run_ph4(infile, N, t_end, delta_t, n_workers,
             use_gpu, gpu_worker, gpu_id,
             accuracy_parameter, softening_length,
             manage_encounters)
