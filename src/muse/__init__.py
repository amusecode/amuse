"""
MUSE framework: The MUlti-physics Software Environment

The aim of MUSE is to provide a software framework in which existing codes 
can be accessed in a homogeneous way in order to perform state-of-the-art 
simulations of a wide range of different physical systems. It serves as the 
base framework for AMUSE, OMUSE etc

It contains several packages, most notably:
units     - MUSE uses quantities, i.e. a number (or array) with a unit attached, 
            instead of just numbers (vital when coupling different codes!)
datamodel - defines particles and grids, on which data (quantities) can be stored
io        - how to read and write data in several file formats

"""
import numpy

def numpy_fix():
    try:
      numpy.set_printoptions(legacy='1.13')
    except TypeError:
      pass
      
numpy_fix()
