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
