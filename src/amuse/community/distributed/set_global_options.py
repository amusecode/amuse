from amuse.lab import *

from amuse.support import options
from amuse.community.seba.interface import SeBa

options.GlobalOptions.instance().override_value_for_option("channel_type", "blaa")

#stellar_evolution = SeBa()

del options.GlobalOptions.instance().overriden_options["channel_type"]

stellar_evolution = SeBa()


