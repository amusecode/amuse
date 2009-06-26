import numpy as numpy
import random
cimport numpy as numpy
import math as m


DTYPE = numpy.float64

ctypedef numpy.float64_t DTYPE_t

def calculate_radius(int n, DTYPE_t mass_cutoff):
    cdef numpy.ndarray result = numpy.zeros(n, dtype=DTYPE)
    cdef DTYPE_t mass_min
    cdef DTYPE_t mass_max
    for x in range(n):
        mass_min = (x * mass_cutoff) / n
        mass_max = ((x+1) * mass_cutoff) / n
        random_mass_fraction = random.uniform(mass_min, mass_max)
        result[x] = 1.0 / m.sqrt( m.pow (random_mass_fraction, -2.0/3.0) - 1.0)
    return result	
    
def calculate_potential_energy(numpy.ndarray[DTYPE_t, ndim=2] mass, numpy.ndarray[DTYPE_t, ndim=2] pos):
    cdef numpy.ndarray[DTYPE_t] result = numpy.zeros(pos.shape[0])
    cdef numpy.ndarray[DTYPE_t, ndim=1] x
    cdef numpy.ndarray[DTYPE_t, ndim=2] delta
    cdef numpy.ndarray[DTYPE_t, ndim=1] r
    cdef numpy.ndarray[DTYPE_t, ndim=1] p
    cdef DTYPE_t m_
    cdef unsigned long n = pos.shape[0]
    cdef unsigned long i = 0
    for i from 0 <= i < n:
        x = pos[i]
        delta = x - pos
        r = numpy.sqrt(numpy.add.reduce(delta * delta,axis=1))
        m_ = mass[i]
        r[i] = 1.0
        p = m_ / r
        p[i] = 0.0
        result[i] = numpy.sum(p)
    return result
