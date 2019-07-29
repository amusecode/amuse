"""
The Astrophysical Multipurpose Software Environment

The aim of AMUSE is to provide a software framework, in which existing codes for
dynamics, stellar evolution, hydrodynamics and radiative transfer can easily be
coupled, in order to perform state-of-the-art simulations of a wide range of
different astrophysical phenomena.

It contains several packages, most notably:
units     - AMUSE uses quantities, i.e. a number (or array) with a unit attached,
            instead of just numbers (vital when coupling different codes!)
datamodel - defines particles and grids, on which data (quantities) can be stored
ic        - a collection of routines to generate initial conditions
io        - how to read and write data in several file formats
community - a variety of existing codes of different physical domains, each with
            a uniform interface to enable coupling

Help is available for each of these packages, e.g.:
> python
>>> import amuse.ic
>>> help(amuse.ic) # doctest: +ELLIPSIS
Help on package amuse.ic in amuse:
...

or (directly from the terminal):
> pydoc amuse.ic
"""
import numpy

from amuse.io import read_set_from_file, write_set_to_file
from amuse.units import (
    units, nbody_system, constants, generic_unit_system,
)
from amuse.units.quantities import (
    ScalarQuantity, VectorQuantity, AdaptingVectorQuantity,
    zero
)
from amuse.datamodel import (
    Particle, Particles, ParticlesSuperset, Grid,
    particle_attributes,
)
from amuse.support.console import (
    set_printing_strategy,
    get_current_printing_strategy,
    set_preferred_units,
)

def numpy_fix():
    """
    Require 1.13 style printing mode
    see https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
    """
    try:
        numpy.set_printoptions(legacy='1.13')
    except TypeError:
        pass

numpy_fix()

try:
    from . import config
except:
    raise ImportError("amuse is not configured correctly")
