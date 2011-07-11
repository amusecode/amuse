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

from amuse.community.kepler.interface import Kepler

def test1():
    mass = 1 | nbody_system.mass
    semi = 1 | nbody_system.length
    ecc = 0.5 | units.none
    time = 5.0 | nbody_system.time
    run_kepler(mass, semi, ecc, time)
    
def run_kepler(mass, semi, ecc, time):

    kep = Kepler(redirection = "none")
    kep.initialize_code()

    kep.set_longitudinal_unit_vector(1.0|units.none, 1.0|units.none,
                                     0.0|units.none)
    kep.initialize_from_elements(mass, semi, ecc)
    a,e = kep.get_elements()
    p = kep.get_periastron()
    print "elements:", a, e, p
    kep.transform_to_time(time)
    x,y,z = kep.get_separation_vector()
    print "separation:", x,y,z
    x,y,z = kep.get_longitudinal_unit_vector()
    print "longitudinal:", x,y,z

    pos = [1, 0, 0] | nbody_system.length
    vel = [0, 0.5, 0] | nbody_system.speed
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    p = kep.get_periastron()
    print "elements:", a, e, p
    kep.transform_to_time(time)
    x,y,z = kep.get_separation_vector()
    print "separation:", x,y,z
    x,y,z = kep.get_longitudinal_unit_vector()
    print "longitudinal:", x,y,z

    kep.stop()

if __name__ == '__main__':

    mass = 1 | nbody_system.mass
    semi = 1 | nbody_system.length
    ecc = 0.5 | units.none
    time = 5.0 | nbody_system.time

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:e:m:t:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            semi = float(a) | nbody_system.length
        elif o == "-e":
            ecc = float(a)
        elif o == "-m":
            time = float(a) | nbody_system.mass
        elif o == "-t":
            time = float(a) | nbody_system.time
        else:
            print "unexpected argument", o

    assert is_mpd_running()
    run_kepler(mass, semi, ecc, time)
