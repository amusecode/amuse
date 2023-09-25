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

def get_amuse_package_dir():
    filename_of_this_script = __file__
    directory_of_this_script = os.path.dirname(os.path.dirname(filename_of_this_script))
    if os.path.isabs(directory_of_this_script):
        return directory_of_this_script
    else:
        return os.path.abspath(directory_of_this_script)
