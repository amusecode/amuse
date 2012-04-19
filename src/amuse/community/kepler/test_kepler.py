import sys
import unittest
import numpy
import random
import collections
import getopt
import os



from amuse.units import nbody_system
from amuse.units import units
from amuse.community.kepler.interface import Kepler

from amuse import datamodel
from amuse.datamodel import particle_attributes
from amuse.rfi.core import is_mpd_running

def test1():
    mass = 1 | nbody_system.mass
    semi = 1 | nbody_system.length
    ecc = 0.5 
    time = 5.0 | nbody_system.time
    run_kepler(mass, semi, ecc, time)
    
def run_kepler(mass, semi, ecc, time):

    kep = Kepler()
    kep.initialize_code()

    kep.set_longitudinal_unit_vector(1.0, 1.0,
                                     0.0)
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
    x,y,z = kep.get_velocity_vector()
    print "velocity:", x,y,z
    x,y,z = kep.get_longitudinal_unit_vector()
    print "longitudinal:", x,y,z

    kep.stop()

if __name__ == '__main__':

    mass = 1 | nbody_system.mass
    semi = 1 | nbody_system.length
    ecc = 0.5 
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
