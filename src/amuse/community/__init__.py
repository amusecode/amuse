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
from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification
from amuse.rfi.core import is_mpd_running

import os


ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
INDEX = MethodWithUnitsDefinition.INDEX
LINK = MethodWithUnitsDefinition.LINK

"""
Existing, production codes

Contains the source code of production codes and software to embed these codes into AMUSE
"""

class _Defaults(OptionalAttributes):
    
    @option(sections=['data'])
    def amuse_root_dir(self):
        if 'AMUSE_DIR' in os.environ:
            return os.environ['AMUSE_DIR']    
        previous = None
        result = os.path.abspath(__file__)
        while not os.path.exists(os.path.join(result,'build.py')):
            result = os.path.dirname(result)
            if result == previous:
                raise exceptions.AmuseException("Could not locate AMUSE root directory!")
            previous = result
        return result

def get_amuse_root_dir():
    return _Defaults().amuse_root_dir
    
