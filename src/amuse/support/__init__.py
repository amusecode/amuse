import os

from amuse.support.options import option, GlobalOptions
from amuse.support.options import OptionalAttributes

import warnings
"""
Support Code


"""

GlobalOptions=GlobalOptions.instance()

class _Defaults(OptionalAttributes):
    
    @option(sections=['data'])
    def amuse_root_dir(self):
        if 'AMUSE_DIR' in os.environ:
            return os.environ['AMUSE_DIR']    

        return GlobalOptions.amuse_data_location

def get_amuse_root_dir():
    return _Defaults().amuse_root_dir

def get_amuse_data_dir():
    return _Defaults().amuse_root_dir
