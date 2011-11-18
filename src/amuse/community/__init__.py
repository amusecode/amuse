from amuse.support.codes.stopping_conditions import StoppingConditionInterface, StoppingConditions

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
from amuse.rfi.core import legacy_function
from amuse.rfi.core import legacy_global
from amuse.rfi.core import LegacyFunctionSpecification
from amuse.rfi.core import is_mpd_running

import os

"""
Existing, production codes

Contains the source code of production codes and software to embed these codes into AMUSE
"""

def get_amuse_root_dir():
    if not 'amuse_root_dir' in locals():
        import os
        try:
            amuse_root_dir = os.environ['AMUSE_DIR']
        except KeyError:
            amuse_root_dir = os.path.abspath(__file__)
            while not os.path.exists(os.path.join(amuse_root_dir, 'build.py')):
                (amuse_root_dir, subdir) = os.path.split(amuse_root_dir)
                if not subdir:
                    raise exceptions.AmuseException("Could not locate AMUSE root directory!")
    return amuse_root_dir
    
def ensure_data_directory_exists(directory):
    directory = os.path.expanduser(directory)
    directory = os.path.expandvars(directory)
    
    print __file__
    print "will make dir:", directory
    print "will make dir:", os.path.abspath(directory)
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
    
    print stack_to_make
    for x in reversed(stack_to_make):
        if(x):
            os.mkdir(x)
    
