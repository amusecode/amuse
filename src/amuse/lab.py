"""
This is the public interface to the AMUSE 
*Astrophysical Multipurpose Software Environment* framework.


"""

from amuse.support.core import late

from amuse.support.units import units
from amuse.support.units import core
from amuse.support.units import si
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.units import nbody_system as nbody

from amuse.support.data import particle_attributes
from amuse.support.data.core import Particle, Particles
from amuse.support.data.values import zero, ScalarQuantity, VectorQuantity,  AdaptingVectorQuantity, new_quantity

from amuse.support.io import write_set_to_file, read_set_from_file, get_options_for_format

from amuse.ext.plummer import new_plummer_sphere
from amuse.ext.salpeter import new_salpeter_mass_distribution, new_salpeter_mass_distribution_nbody

from amuse.legacy.bhtree.interface import BHTree, BHTreeInterface
from amuse.legacy.hermite0.interface import Hermite, HermiteInterface
from amuse.legacy.phiGRAPE.interface import PhiGRAPE, PhiGRAPEInterface
from amuse.legacy.octgrav.interface import Octgrav, OctgravInterface
from amuse.legacy.twobody.twobody import TwoBody, TwoBodyInterface
from amuse.legacy.fi.interface import Fi, FiInterface

from amuse.legacy.sse.interface import SSE, SSEInterface
from amuse.legacy.bse.interface import BSE, BSEInterface
from amuse.legacy.evtwin.interface import EVtwin, EVtwinInterface
from amuse.legacy.mesa.interface import MESA, MESAInterface



def vector(value = [], unit = None):
    if unit is None:
        if isinstance(value, core.unit):
            return VectorQuantity([], unit = value)
        elif isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            return AdaptingVectorQuantity(value)
    else:
        if isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            return VectorQuantity(value, unit)
            
