import numpy as np
import random
cimport numpy as np
import math as m


DTYPE = np.float64

ctypedef np.float64_t DTYPE_t

def calculate_radius(int n, DTYPE_t mass_cutoff):
    cdef np.ndarray result = np.zeros(n, dtype=DTYPE)
    cdef DTYPE_t mass_min
    cdef DTYPE_t mass_max
    for x in range(n):
        mass_min = (x * mass_cutoff) / n
        mass_max = ((x+1) * mass_cutoff) / n
        random_mass_fraction = random.uniform(mass_min, mass_max)
        result[x] = 1.0 / m.sqrt( m.pow (random_mass_fraction, -2.0/3.0) - 1.0)
    return result	
