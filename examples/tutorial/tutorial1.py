# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# The AMUSE framework is divided over several python packages and modules. All these have one common root module: _amuse_. This root defines the namespace of all sub-packages and modules. You can load the amuse module with:

# <codecell>

import amuse

# <markdowncell>

# However the _amuse_ module is empty, as you can see when you try 'Tab' completion on the amuse module or run dir(amuse)

# <codecell>

dir(amuse)

# <markdowncell>

# The main body of the AMUSE framework is divided over 4 packages, each having subdivided into more packages and modules.
# 
# ### community ###
# 
# This packages contains all the communtity codes. All software integrated into AMUSE, (N-body, stellar evolution, hydrodynamics and radiative transfer codes) is called _community_ code in AMUSE. Each community code is defined in a separate sub-package and every sub-package contains at least one module called `interface`. You can load a community code with: `from amuse.community.<codename>.interface import <codeclass>`. In later tutorials we will learn more about the codes and how to interact with each.

# <codecell>

from amuse.community.bhtree.interface import BHTree
BHTree

# <markdowncell>

# ### units ###
# 
# A package to work with quantities and units. All calculations in AMUSE are done with quantities having units. These quantities and their units are implemented as python classes and can be used almost everywere you would normaly use a number (or a `numpy` array). In the next tutorial we will come back to the units, for now we will show a simple example 

# <codecell>

from amuse.units import units
from amuse.units import constants

# <codecell>

constants.G * ( 5.972e24 | units.kg) /  (6378.1 | units.km )**2

# <markdowncell>

# ### datamodel ###
# 
# All astrophysical bodies (stars, clouds, planets etc) are modelled with sets of particles or on grids. These sets and grids are defined in the _datamodel_ package. You will see these used everywhere in amuse and several tutorial cover them in more detail. 

# <codecell>

from amuse.datamodel import Particles
solar_system_planets = Particles(7)
solar_system_planets.mass = [641.85, 4868.5, 5973.6, 102430, 86832, 568460, 1898600] | (1e21 * units.kg)
print solar_system_planets

# <markdowncell>

# ### rfi ###
# 
# The AMUSE framework is written in Python, most codes are written in C or Fortran. In AMUSE the Remote Function Invocation or _rfi_ package implements all classes and tools to call functions on the community codes. The _rfi_ package implements support for communicating over MPI, raw sockets and Ibis. This code is mainly used internally to the framework and in most scripts you will not see it, but you will use it!

# <codecell>

from amuse.rfi.channel import is_mpd_running
print is_mpd_running()

# <markdowncell>

# ### couple ###

# <codecell>


# <markdowncell>

# ### io ###

# <codecell>


# <markdowncell>

# ### ext ###

# <codecell>


# <markdowncell>

# ### ic ###

# <codecell>


