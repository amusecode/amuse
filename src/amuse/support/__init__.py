import os

from amuse.support.options import option
from amuse.support.options import OptionalAttributes


"""
Support Code


"""

class _Defaults(OptionalAttributes):
    
    @option(sections=['data'])
    def amuse_root_dir(self):
        if 'AMUSE_DIR' in os.environ:
            return os.environ['AMUSE_DIR']    

        this = os.path.dirname(os.path.abspath(__file__))
        
        # installed
        result=os.path.abspath(os.path.join(this, "..","..","..","..","..","share", "amuse"))
        if os.path.exists(os.path.join(result,'build.py')):
            return result
        
        # in-place
        result=os.path.abspath(os.path.join(this, "..","..",".."))        
        if os.path.exists(os.path.join(result,'build.py')):
            return result

        raise exceptions.AmuseException("Could not locate AMUSE root directory! set the AMUSE_DIR variable")

def get_amuse_root_dir():
    return _Defaults().amuse_root_dir
