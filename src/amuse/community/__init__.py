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
    
def ensure_data_directory_exists(directory):
    directory = os.path.expanduser(directory)
    directory = os.path.expandvars(directory)
    
    if os.path.exists(directory):
        if os.path.isdir(directory):
            return
        else:
            raise exceptions.AmuseException("Path exists but is not a directory {0}".format(directory))
    
    stack_to_make = [directory]
    previous = None
    current = os.path.dirname(directory)
    while previous != current:
        if not os.path.exists(current):
            stack_to_make.append(current)
        else:
            if not os.path.isdir(current):
                raise exceptions.AmuseException("Path exists but is not a directory {0}".format(current))
            break
        previous = current
        current = os.path.dirname(current)
    
    for x in reversed(stack_to_make):
        if(x):
            os.mkdir(x)
    
