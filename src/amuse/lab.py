"""
This is the public interface to the AMUSE 
*Astrophysical Multipurpose Software Environment* framework.


"""

from amuse.support.core import *

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system as nbody

from amuse.support.data import particle_attributes

from amuse.ext.plummer import new_plummer_sphere
from amuse.ext.salpeter import new_salpeter_mass_distribution, new_salpeter_mass_distribution_nbody

from amuse.legacy.bhtree.interface import BHTree, BHTreeInterface
from amuse.legacy.hermite0.interface import Hermite, HermiteInterface
from amuse.legacy.phiGRAPE.interface import PhiGRAPE, PhiGRAPEInterface
from amuse.legacy.octgrav.interface import Octgrav, OctgravInterface
from amuse.legacy.twobody.twobody import TwoBody, TwoBodyInterface
from amuse.legacy.fi.interface import Fi, fi

from amuse.legacy.sse.interface import SSE, SSEInterface
from amuse.legacy.bse.interface import BSE, BSEInterface
from amuse.legacy.evtwin.interface import EVtwin, EVtwinInterface
