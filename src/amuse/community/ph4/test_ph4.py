import sys
import unittest
import numpy
import random
import collections
import getopt
import os

# Test code borrowed from test_hermite0.py.

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import particle_attributes

from amuse.community.ph4.interface import ph4 as grav
#from amuse.community.phiGRAPE.interface import PhiGRAPE as grav
#from amuse.community.hermite0.interface import Hermite as grav

from amuse.support.codes.core import is_mpd_running
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import new_salpeter_mass_distribution_nbody

# Control unit output:

def mm(mass): return mass.value_in(nbody_system.mass)
def ll(length): return length.value_in(nbody_system.length)
def tt(time): return time.value_in(nbody_system.time)
def ee(energy): return energy.value_in(nbody_system.energy)
def nn(none): return none.value_in(units.none)

def print_log(time, gravity, particles, initial_total_energy, total_energy):
    print "time =", tt(time), " energy = ", ee(total_energy), " dE/E0 = ", \
        nn((total_energy - initial_total_energy) / initial_total_energy)

def test_ph4(number_of_stars = 50,
             end_time = 10 | nbody_system.time,
             delta_t = 1 | nbody_system.time,
             n_workers = 1, use_gpu = 1,
             accuracy_parameter = 0.1,
             softening_length = -1 | nbody_system.length):

    print "number_of_stars =", number_of_stars
    print "end_time =", end_time
    print "delta_t =", delta_t
    print "n_workers =", n_workers
    print "use_gpu =", use_gpu

    print "initializing the gravity module"
    gravity = grav(number_of_workers = n_workers,
                   redirection = "none",
                   mode = "gpu")
    gravity.initialize_code()
    gravity.parameters.set_defaults()

    if softening_length == -1 | nbody_system.length:
        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
			| nbody_system.length**2
    else:
        eps2 = softening_length*softening_length

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu

    print "making a Plummer model"
    particles = MakePlummerModel(number_of_stars).result
    # print particles

    print "setting particle masses and radii"
    # particles.mass = (1.0 / number_of_stars) | nbody_system.mass
    scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars) 
    particles.mass = scaled_mass
    particles.radius = 0.0 | nbody_system.length

    print "centering particles"
    particles.move_to_center()
    print "scaling particles to virial equilibrium"
    particles.scale_to_standard(smoothing_length_squared
                                = gravity.parameters.epsilon_squared)

    print "adding particles"
    gravity.particles.add_particles(particles)
    gravity.commit_particles()

    time = 0.0 | nbody_system.time
    kinetic_energy = gravity.kinetic_energy
    potential_energy = gravity.potential_energy
    total_energy = initial_total_energy = kinetic_energy + potential_energy

    print ''
    print_log(time, gravity, particles, initial_total_energy, total_energy)
    print "initial virial ratio =", nn(-2*kinetic_energy/potential_energy)
    print "evolving to time =", tt(end_time), "in steps of", tt(delta_t)

    while time < end_time:
        time += delta_t
        gravity.evolve_model(time)
        total_energy = gravity.kinetic_energy + gravity.potential_energy
        print_log(time, gravity, particles, initial_total_energy, total_energy)

    print ''
    gravity.stop()

if __name__ == '__main__':

    N = 100
    t_end = 5.0 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 2
    use_gpu = 1
    accuracy_parameter = 0.1
    softening_length = -1  | nbody_system.length

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:d:e:gn:t:w:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time 
        elif o == "-e":
            softening_length = float(a) | nbody_system.length
        elif o == "-g":
            use_gpu = 0
        elif o == "-n":
            N = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print "unexpected argument", o

    assert is_mpd_running()
    test_ph4(N, t_end, delta_t, n_workers, use_gpu,
             accuracy_parameter, softening_length)
