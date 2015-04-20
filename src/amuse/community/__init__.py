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

import os

from amuse.support import get_amuse_root_dir

ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
INDEX = MethodWithUnitsDefinition.INDEX
LINK = MethodWithUnitsDefinition.LINK

"""
Existing, production codes

Contains the source code of production codes and software to embed these codes into AMUSE
"""

    
