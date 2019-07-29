"""
Initial Conditions

This module contains a number of packages with routines to generate initial conditions.

"""

from amuse.ic.plummer import new_plummer_model, new_plummer_sphere
from amuse.ic.salpeter import (
    new_salpeter_mass_distribution, new_salpeter_mass_distribution_nbody,
    new_powerlaw_mass_distribution, new_powerlaw_mass_distribution_nbody,
)
from amuse.ic.brokenimf import (
    new_broken_power_law_mass_distribution, new_scalo_mass_distribution,
    new_miller_scalo_mass_distribution, new_kroupa_mass_distribution,
)
from amuse.ic.flatimf import (
    new_flat_mass_distribution, new_flat_mass_distribution_nbody,
)
from amuse.ic.kingmodel import new_king_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.molecular_cloud import new_ism_cube
