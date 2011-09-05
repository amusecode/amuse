from amuse.support.data.base import *
from amuse.support.data.memory_storage import *
from amuse.support.data.particles import *
from amuse.support.data import particle_attributes
from amuse.support.data.grids import *
from amuse.support.data import grid_attributes

"""
This module provides access to all set handling 
in AMUSE. The actual implementation is in the
base, storage and particle modules.
"""

import warnings

warnings.warn("amuse.support.data.core is deprecated, use amuse.support.data instead", DeprecationWarning)


