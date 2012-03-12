import sys
import os.path
from StringIO import StringIO

from amuse.support import options

amuserc_template = """
[data]
input_data_root_directory={0}
"""

def setup_directories():
    if not options.GlobalOptions.INSTANCE is None:
        options.GlobalOptions.INSTANCE = None
    
    input_data_root_directory = os.path.join(sys.prefix, "share", "amuse", "data")
    amuserc = amuserc_template.format(input_data_root_directory)
    fp = StringIO(amuserc)
    global_options = options.GlobalOptions.instance(fp)

setup_directories()