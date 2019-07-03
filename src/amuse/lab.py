"""
This is the public interface to the AMUSE 
*Astrophysical Multipurpose Software Environment* framework.


"""

from amuse.support.core import late


from amuse.units import units
from amuse.units import core
from amuse.units import si
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units import generic_unit_converter
from amuse.units import generic_unit_system as generic
from amuse.units import nbody_system as nbody
from amuse.units.quantities import zero
from amuse.units.quantities import ScalarQuantity
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units.quantities import new_quantity

from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.io import get_options_for_format
from amuse.io import ReportTable

from amuse.ext.solarsystem import new_solar_system_for_mercury, new_solar_system
from amuse.ext.halogen_model import new_halogen_model
from amuse.ext.galactics_model import new_galactics_model
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution, new_spherical_particle_distribution
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH, pickle_stellar_model

_community_codes=[
    "BHTree", 
    "Hermite", 
    "PhiGRAPE",
    "Octgrav",
    "TwoBody",
    "Huayno",
    "ph4",
    "Bonsai",
    "Pikachu",
    "AarsethZare",
    "Adaptb",
    "Hacs64",
    "HiGPUs",
    "Kepler",
    "Mercury",
    "MI6",
    "Mikkola",
    "SmallN",
    "Rebound",
    "Brutus",
    "Fi",
    "Gadget2",
    "Athena",
    "Capreole",
    "MpiAmrVac",
    "SimpleX",
    "Mocassin",
    "SPHRay",
    "SSE",
    "BSE",
    "MOSSE",
    "MOBSE",
    "SeBa",
    "EVtwin",
    "MESA",
    "MMAMS",
    "Hop",
    ]


def _placeholder(name):
    class _placeholder(object):
        def __init__(self, *arg,**kwargs):
            raise Exception("failed import, code {0} not available, maybe it needs to be (pip) installed?".format(name))
    return _placeholder

for _name in _community_codes:
    _interfacename=_name+"Interface"
    _packagename=_name.lower()
    try:
        _interface=__import__("amuse.community."+_packagename+".interface", fromlist=[_name, _interfacename])
        locals()[_name]=getattr(_interface,_name)
        locals()[_interfacename]=getattr(_interface,_interfacename)
    except ImportError:
        locals()[_name]=_placeholder(_packagename)
        locals()[_interfacename]=_placeholder(_packagename)

from amuse.support.console import set_printing_strategy
from amuse.support.console import get_current_printing_strategy
from amuse.datamodel import particle_attributes
from amuse.datamodel import Particle
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.datamodel import Grid

from amuse.ic.plummer import new_plummer_model, new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.ic.salpeter import new_powerlaw_mass_distribution
from amuse.ic.salpeter import new_powerlaw_mass_distribution_nbody
from amuse.ic.brokenimf import new_broken_power_law_mass_distribution, new_scalo_mass_distribution
from amuse.ic.brokenimf import new_miller_scalo_mass_distribution, new_kroupa_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution, new_flat_mass_distribution_nbody
from amuse.ic.kingmodel import new_king_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.molecular_cloud import new_ism_cube

def vector(value = [], unit = None):
    if unit is None:
        if isinstance(value, core.unit):
            return VectorQuantity([], unit = value)
        elif isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            result = AdaptingVectorQuantity()
            result.extend(value)
            return result
    else:
        if isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            return VectorQuantity(value, unit)
            

