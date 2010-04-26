"""
Exisiting, production coodes

Contains the source code of production codes and software to embed these codes into AMUSE
"""

from amuse.legacy.support.core import LegacyInterface, legacy_function, legacy_global
from amuse.legacy.support.core import LegacyFunctionSpecification, is_mpd_running
from amuse.support.data import parameters
from amuse.support.data import attributes
from amuse.support.units import units
from amuse.support.units import nbody_system

from amuse.support.interface import *

def get_amuse_root_dir():
    import os
    try:
        amuse_root_dir = os.environ['AMUSE_ROOT_DIR']
    except KeyError:
        amuse_root_dir = os.path.abspath(__file__)
        while not os.path.exists(os.path.join(amuse_root_dir, 'build.py')):
            amuse_root_dir = os.path.dirname(amuse_root_dir)
    return amuse_root_dir
