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

def numpy_fix():
    try:
      numpy.set_printoptions(legacy='1.13')
    except TypeError:
      pass
      
numpy_fix()

class NoConfig(object):
    def __init__(self, message):
        self._message=message
    def __getattr__(self, attr):
        raise AttributeError(self._message)

try:
    from . import config
except Exception as ex:
    message="Configuration not read in - or configuration invalid, exception:\n"+str(ex)
    config=NoConfig(message)
    
