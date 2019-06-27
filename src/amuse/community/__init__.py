import os
import sys
import imp

from amuse.community.interface.stopping_conditions import StoppingConditionInterface
from amuse.community.interface.stopping_conditions import StoppingConditions

from amuse.support.options import option
from amuse.support.options import OptionalAttributes

from amuse.units import units
from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support import exceptions

from amuse.support.interface import *

from amuse.datamodel import parameters
from amuse.datamodel import attributes
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface
from amuse.rfi.core import CodeWithDataDirectories
from amuse.rfi.core import legacy_function,remote_function
from amuse.rfi.core import LegacyFunctionSpecification
from amuse.rfi.core import is_mpd_running

from amuse.support import get_amuse_root_dir

ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
INDEX = MethodWithUnitsDefinition.INDEX
LINK = MethodWithUnitsDefinition.LINK

"""
Existing, production codes

Contains the source code of production codes and software to embed these codes into AMUSE
"""

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
    if not os.path.isfile(os.path.join(os.path.dirname(__file__),_packagename,"interface.py")):
        _package=imp.new_module(_packagename)
        _interface=imp.new_module("interface")

        _interface.__dict__[_name]=_placeholder(_packagename)
        _interface.__dict__[_interfacename]=_placeholder(_packagename)
        _package.__dict__["interface"]=_interface

        sys.modules[__name__+'.'+_packagename]=_package
        sys.modules[__name__+'.'+_packagename+'.interface']=_interface
    
