import collections
import numpy
import operator
import os
import random
import sys
import unittest

from math import sqrt

from amuse.units import nbody_system
from amuse.units import units

from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
MassFraction = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] \
    | units.none

# This could be more elegant in amuse as an attribute
# for example as: stars.position.length_sq()
# at the moment (April 2010) this is not yet implemented, but the alternative:
# star.position.length()
# is slow, due to the implicite square-root
def distance_sq(stars) :
    return stars.x**2  + stars.y**2  +stars.z**2

def LagrangianRadii(stars, verbose = 0, massf = MassFraction) :

    com = stars.center_of_mass()
    stars.position = stars.position - com

    vcom = stars.center_of_mass_velocity()
    stars.velocity = stars.velocity - vcom

    n_stars = len(stars)
    # Next line is a potential performance bottleneck, becasue
    # the for loop is not in numpy but explicit
    # old but slow: d2 = numpy.array([distance_sq(star) for star in stars])
    d2 = distance_sq(stars)
    m = stars.mass/stars.mass.sum()
    d2m = zip(d2, m)
    d2m.sort(key = operator.itemgetter(0))

    iL = 0
    mt = 0 | units.none
    Lagrad = []
    for d2i,mi in d2m :
        mt += mi
        while mt >= massf[iL] :
            Lagrad.append(d2i.sqrt())
            # Previous line is preferable above here below, beacuse
            #  the former preserves the units
            #  Lagrad.append(sqrt(d2[ni].value_in(nbody_system.length**2)))
            if verbose==1 :
                print "Lagrangian Radius M= ", mt, \
                      "(iL=", iL, ") at d= ", Lagrad[-1]
            iL += 1
            if iL >= len(massf) :
                break
    return Lagrad

if __name__ == '__main__' :

    assert is_mpd_running()
    seed = None
#   seed = numpy.random.RandomState([1,1,1])

    nstars = 128
    if len(sys.argv)>1 :
        stars = int(sys.argv[1])
    with_units = len(sys.argv) > 2

    if not with_units :
        mass_unit = nbody_system.mass
        length_unit = nbody_system.length
    else :
        mass_unit = units.MSun
        length_unit = units.parsec

    m_min = 0.1 | mass_unit
    m_max = 100 | mass_unit
    alpha = -2.35

    r_vir = 1 | length_unit
    masses = new_salpeter_mass_distribution(nstars, m_min, m_max, alpha)
    m_tot = masses.sum()

    if not with_units :
        convert_nbody = None
        masses /= m_tot.value_in(nbody_system.mass)     # scale to unit mass 
        m_tot = 1 | nbody_system.mass
    else :
        convert_nbody = nbody_system.nbody_to_si(m_tot, r_vir)
        convert_nbody.set_as_default()
        print m_tot

    stars = new_plummer_sphere(nstars, convert_nbody, random_state = seed);
    stars.mass = masses 
    
    LagrangianRadii(stars, verbose=1)
