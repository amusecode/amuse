import collections
import getopt
import numpy
import os
import random
import sys
import unittest
import pickle

from time import clock
from time import gmtime
from time import mktime

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

def make_nbody(number_of_stars = 100, time = 0.0,
               n_workers = 1, use_gpu = 1, gpu_worker = 1,
               salpeter = 0,
               delta_t = 1.0 | nbody_system.time,
               timestep_parameter = 0.1,
               softening_length = 0.0 | nbody_system.length,
               random_seed = 1234):

    # Make an N-body system, print out some statistics on it, and save
    # it in a restart file.  The restart file name is of the form
    # 't=nnnn.n.xxx', where the default time is 0.0.

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2,31)-1)
    numpy.random.seed(random_seed)
    print("random seed =", random_seed)

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

    #-----------------------------------------------------------------

    # Make a standard N-body system.

    print("making a Plummer model")
    stars = new_plummer_model(number_of_stars)

    id = numpy.arange(number_of_stars)
    stars.id = id+1

    print("setting particle masses and radii")
    if salpeter == 0:
        print('equal masses')
        total_mass = 1.0 | nbody_system.mass
        scaled_mass = total_mass / number_of_stars
    else:
        print('salpeter mass function')
        mmin = 0.5 | nbody_system.mass
        mmax = 10.0 | nbody_system.mass
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars,
                                                           mass_min = mmin,
                                                           mass_max = mmax)
        
    stars.mass = scaled_mass

    print("centering stars")
    stars.move_to_center()
    print("scaling stars to virial equilibrium")
    stars.scale_to_standard(smoothing_length_squared
                            = gravity.parameters.epsilon_squared)

    # Set dynamical radii (assuming virial equilibrium and standard
    # units).  Note that this choice should be refined, and updated
    # as the system evolves.  Probably the choice of radius should be
    # made entirely in the multiples module.  TODO.  In these units,
    # M = 1 and <v^2> = 0.5, so the mean 90-degree turnaround impact
    # parameter is
    #
    #		b_90 = G (m_1+m_2) / vrel^2
    #		     = 2 <m> / 2<v^2>
    #		     = 2 / N			for equal masses
    #
    # Taking r_i = m_i / 2<v^2> = m_i in virial equilibrium means
    # that, approximately, "contact" means a 90-degree deflection (r_1
    # + r_2 = b_90).  A more conservative choice with r_i less than
    # this value will isolate encounters better, but also place more
    # load on the large-N dynamical module.

    stars.radius = 0.5*stars.mass.number | nbody_system.length

    time = 0.0 | nbody_system.time
    # print "IDs:", stars.id.number

    print("recentering stars")
    stars.move_to_center()
    sys.stdout.flush()

    #-----------------------------------------------------------------

    if softening_length < 0.0 | nbody_system.length:

        # Use ~interparticle spacing.  Assuming standard units here.  TODO

        softening_length = 0.5*float(number_of_stars)**(-0.3333333) \
				| nbody_system.length
    print('softening length =', softening_length)

    gravity.parameters.timestep_parameter = timestep_parameter
    gravity.parameters.epsilon_squared = softening_length*softening_length
    gravity.parameters.use_gpu = use_gpu

    print('')
    print("adding particles")
    # print stars
    sys.stdout.flush()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print('')
    print("number_of_stars =", number_of_stars)
    sys.stdout.flush()

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    # -----------------------------------------------------------------
    # Create the coupled code and integrate the system to the desired
    # time, managing interactions internally.

    kep = init_kepler(stars[0], stars[1])
    multiples_code = multiples.Multiples(gravity, new_smalln, kep)

    multiples_code.neighbor_perturbation_limit = 0.1
    multiples_code.neighbor_veto = True

    print('')
    print('multiples_code.initial_scale_factor =', \
        multiples_code.initial_scale_factor)
    print('multiples_code.neighbor_perturbation_limit =', \
        multiples_code.neighbor_perturbation_limit)
    print('multiples_code.neighbor_veto =', \
        multiples_code.neighbor_veto)
    print('multiples_code.final_scale_factor =', \
        multiples_code.final_scale_factor)
    print('multiples_code.initial_scatter_factor =', \
        multiples_code.initial_scatter_factor)
    print('multiples_code.final_scatter_factor =', \
        multiples_code.final_scatter_factor)
    print('multiples_code.retain_binary_apocenter =', \
        multiples_code.retain_binary_apocenter)
    print('multiples_code.wide_perturbation_limit =', \
        multiples_code.wide_perturbation_limit)

    # Take a dummy step, just in case...

    multiples_code.evolve_model(time)

    # Copy values from the module to the set in memory.

    channel.copy()
    
    # Copy the index (ID) as used in the module to the id field in
    # memory.  The index is not copied by default, as different
    # codes may have different indices for the same particle and
    # we don't want to overwrite silently.

    channel.copy_attribute("index_in_code", "id")

    pre = "%%% "
    E0,cpu0 = print_log(pre, time, multiples_code)
    sys.stdout.flush()

    # file = 't='+'{:07.2f}'.format(time.number)        # fails in Python 2.6
    file = 't=%07.2f'%time.number
    write_state_to_file(time, stars, gravity, multiples_code,
                        file, delta_t, E0, cpu0)
    tree_copy = multiples_code.root_to_tree.copy()

    del multiples_code
    sys.stdout.flush()

    gravity.stop()
    kep.stop()
    stop_smalln()
    print('')


if __name__ == '__main__':

    # Defaults:

    N = 1000
    time = 0.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 1
    use_gpu = 1
    gpu_worker = 1
    salpeter = 0
    timestep_parameter = 0.1
    softening_length = 0.0 | nbody_system.length
    random_seed = -1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "n:st:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    for o, a in opts:
        if o == "-n":
            N = int(a)
        elif o == "-s":
            salpeter = 1
        elif o == "-t":
            time = float(a) | nbody_system.time
        else:
            print("unexpected argument", o)

    assert is_mpd_running()

    make_nbody(N, time, n_workers, use_gpu, gpu_worker,
               salpeter, delta_t, timestep_parameter, softening_length,
               random_seed)
