import collections
import getopt
import numpy
import os
import random
import sys
import pickle
import math

import unittest
from time import clock

from amuse.community.ph4.interface import ph4 as grav
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.couple import multiples

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import quantities

from amuse import datamodel
from amuse.datamodel import particle_attributes as pa
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

from amuse import io
from utils import *

def handle_callback (time, star1, star2):
    print('')
    print('    callback')
    print('   ', time)
    print('   ', star1)
    print('   ', star2)
    print('')
    return True		# NOTE: returning False will skip this encounter

def run_ph4(initial_file = None,
            end_time = 0 | nbody_system.time,
            input_delta_t = 0.0 | nbody_system.time,
            input_Delta_t = 1.0 | nbody_system.time,
            input_timestep_parameter = 0.0,
            input_softening_length = -1.0 | nbody_system.length,
            n_workers = 1, use_gpu = 1, gpu_worker = 1,
            use_multiples = True,
            save_restart = False,
            strict_restart = False
        ):

    # Read an N-body system from a file and run it to the specified
    # time using the specified steps.  Print log information and
    # optionally save a restart file after every step.  If the
    # specified time is less than the time in the initial file, don't
    # take a step, but still print out the log info.  (Hence run_ph4
    # also functions like Starlab sys_stats.)

    print("initial_file =", initial_file)
    print("end_time =", end_time.number)
    print("n_workers =", n_workers)
    print("use_gpu =", use_gpu)
    print("use_multiples =", use_multiples)
    print("save_restart =", save_restart)
    print("strict_restart =", strict_restart)
    print("\ninitializing the gravity module")
    sys.stdout.flush()

    init_smalln()
    
    # Note that there are actually three GPU options:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    if gpu_worker == 1:
        try:
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none", mode = "gpu")
        except Exception as ex:
            gravity = grav(number_of_workers = n_workers,
                           redirection = "none")
    else:
        gravity = grav(number_of_workers = n_workers,
                       redirection = "none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    kep = Kepler(None, redirection = "none")
    kep.initialize_code()

    stars, time, delta_t, E0, cpu0, multiples_code \
        	= read_state_from_file(initial_file, gravity, kep)

    # Allow overrides of the restored data (OK for delta_t, NOT
    # recommended for timestep_parameter or softening_length).  Note
    # that reading the state also commits the particles, and hence
    # calculates the initial time steps.  Probably should reinitialize
    # if timestep_parameter or softening_length are changed.  TODO

    if input_delta_t.number > 0:
        if input_delta_t != delta_t:
            print('modifying delta_t from stored', delta_t, \
		  'to input', input_delta_t)
            delta_t = input_delta_t
    else:
        print("using stored delta_t =", delta_t)

    print(input_timestep_parameter)
    print(gravity.parameters.timestep_parameter)

    if input_timestep_parameter > 0:
        if input_timestep_parameter != gravity.parameters.timestep_parameter:
            print('modifying timestep_parameter from stored', \
            	  gravity.parameters.timestep_parameter, \
		  'to input', input_timestep_parameter)
            gravity.parameters.timestep_parameter \
		= input_timestep_parameter
    else:
        print('timestep_parameter =', gravity.parameters.timestep_parameter)

    if input_softening_length.number >= 0:
        if input_softening_length*input_softening_length \
		!= gravity.parameters.epsilon_squared:
            print('modifying softening_length from stored', \
            	  gravity.parameters.epsilon_squared.sqrt(), \
		  'to input', input_softening_length)
            gravity.parameters.epsilon_squared \
                = softening_length*softening_length
    else:
        print('softening length =', gravity.parameters.epsilon_squared.sqrt())

    gravity.parameters.use_gpu = use_gpu
    gravity.parameters.begin_time = time

    if 0:
        print('')
        print(gravity.parameters.begin_time)
        print(stars.mass)
        #print stars.position
        for s in stars:
            print('%.18e %.18e %.18e' % (s.x.number, s.y.number, s.z.number))
        print(stars.velocity)

    channel = gravity.particles.new_channel_to(stars)

    if use_multiples:
        stopping_condition = gravity.stopping_conditions.collision_detection
        stopping_condition.enable()

    gravity.parameters.force_sync = 1	# end exactly at the specified time

    pre = "%%% "
    print_log(pre, time, multiples_code, E0, cpu0)

    tsave = time + Delta_t
    save_file = ''

    while time < end_time:

        time += delta_t
        multiples_code.evolve_model(time) #, callback=handle_callback)

        # Copy values from the module to the set in memory.

        channel.copy()
    
        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        # Write log information.

        print_log(pre, time, multiples_code, E0, cpu0)
        sys.stdout.flush()

        # Optionally create a restart file.

        if save_restart and time >= tsave:
            #save_file = 't='+'{:07.2f}'.format(time.number) # not in Python 2.6
            save_file = 't=%07.2f'%time.number
            write_state_to_file(time, stars, gravity, multiples_code,
                                save_file, delta_t, E0, cpu0)
            sys.stdout.flush()
            tsave += Delta_t
            if strict_restart: break

    gravity.stop()
    kep.stop()
    stop_smalln()

    return time, save_file


def print_help():
    print("Options:")
    print("    -d set log output interval [1.0]")
    print("    -h print this help message")
    print("    -i set initial file name [t=0000.0]")
    print("    -m suppress multiples  [False]")
    print("    -s save restart files [False]")
    print("    -t set final time [0.0]")

if __name__ == '__main__':

    initial_file = 't=0000.0'
    t_end = 0.0 | nbody_system.time	# default is to print log and exit
    delta_t = 0.0 | nbody_system.time	# log output time scale
    Delta_t = 1.0 | nbody_system.time	# restart output time scale
    Delta_t_set = False
    timestep_parameter = -1
    softening_length = -1  | nbody_system.length
    n_workers = 1
    use_gpu = 1
    gpu_worker = 1
    use_multiples = True
    save_restart = False
    strict_restart = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:D:hi:msSt:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    for o, a in opts:
        if o == "-d":
            delta_t = float(a) | nbody_system.time
        elif o == "-D":
            Delta_t = float(a) | nbody_system.time
            Delta_t_set = True
        elif o == "-h":
            print_help()
            sys.exit(1)
        elif o == "-i":
            initial_file = a
        elif o == "-m":
            use_multiples = False
        elif o == "-s":
            save_restart = True
        elif o == "-S":
            strict_restart = True
            save_restart = True
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        else:
            print("unexpected argument", o)
            print_help()
            sys.exit(1)

    if not save_restart: strict_restart = False
    if strict_restart: save_restart = True
    if not Delta_t_set: Delta_t = delta_t

    # Code implicitly assumes that delta_t and Delta_t are
    # commensurate. Check that here.

    if math.fmod(Delta_t.number, delta_t.number) != 0:
        x = math.floor(Delta_t/delta_t)
        Delta_t = x*delta_t
        print('reset Delta_t to', Delta_t)

    assert is_mpd_running()

    # In non-strict_mode, OK to let run_ph4 loop over steps. If
    # strict_restart is True, handle the loop here and force a new
    # restart at every step -- seems substantially slower, but
    # should be reproducible.

    if (not strict_restart):
        run_ph4(initial_file, t_end, delta_t, Delta_t,
                timestep_parameter, softening_length,
                n_workers, use_gpu, gpu_worker,
                use_multiples,
                save_restart, strict_restart)
    else:
        t = -1.0 | nbody_system.time
        while t < t_end:
            t, initial_file = run_ph4(initial_file, t_end, delta_t, Delta_t,
                                      timestep_parameter, softening_length,
                                      n_workers, use_gpu, gpu_worker,
                                      use_multiples,
                                      save_restart, strict_restart)
            delta_t = 0.0 | nbody_system.time
            timestep_parameter = -1
            softening_length = -1  | nbody_system.length
            print('')
