"""
This file contains conversion functions to/from Astropy.

Requires astropy to be installed.
"""

from amuse.units import units
from amuse.units.si import m, kg, s, A, K, mol, cd

try:
    import astropy.units as apu
    import astropy.constants as apc
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "Astropy does not seem to be installed"
    )


def from_astropy(ap_quantity):
    "Convert a quantity from Astropy to AMUSE"

    # Find SI bases of the unit
    si_bases = ap_quantity.si.unit.bases
    si_powers = ap_quantity.si.unit.powers
    si_units = list(zip(si_powers, si_bases))

    # Find the quantity's value in base units
    si_value = ap_quantity.si.value

    # Reconstruct the quantity in AMUSE units
    amuse_quantity = si_value
    for base_unit in si_units:
        if base_unit[1].name == "m":
            amuse_quantity = amuse_quantity * (1 | units.m**base_unit[0])
        elif base_unit[1].name == "kg":
            amuse_quantity = amuse_quantity * (1 | units.kg**base_unit[0])
        elif base_unit[1].name == "s":
            amuse_quantity = amuse_quantity * (1 | units.s**base_unit[0])
        elif base_unit[1].name == "A":
            amuse_quantity = amuse_quantity * (1 | units.A**base_unit[0])
        elif base_unit[1].name == "K":
            amuse_quantity = amuse_quantity * (1 | units.K**base_unit[0])
        elif base_unit[1].name == "mol":
            amuse_quantity = amuse_quantity * (1 | units.mol**base_unit[0])
        elif base_unit[1].name == "cd":
            amuse_quantity = amuse_quantity * (1 | units.cd**base_unit[0])

    return amuse_quantity


def to_astropy(quantity):
    "Convert a quantity from AMUSE to Astropy"

    # Find the SI bases of the unit
    unit = quantity.unit
    unit_bases = unit.base

    # Find the quantity's value in base units
    value = quantity.value_in(unit.base_unit())
    
    # Reconstruct the quantity in Astropy units
    ap_quantity = value
    for base_unit in unit_bases:
        if base_unit[1] == m:
            ap_quantity = ap_quantity * apu.m**base_unit[0]
        elif base_unit[1] == kg:
            ap_quantity = ap_quantity * apu.kg**base_unit[0]
        elif base_unit[1] == s:
            ap_quantity = ap_quantity * apu.s**base_unit[0]
        elif base_unit[1] == A:
            ap_quantity = ap_quantity * apu.A**base_unit[0]
        elif base_unit[1] == K:
            ap_quantity = ap_quantity * apu.K**base_unit[0]
        elif base_unit[1] == mol:
            ap_quantity = ap_quantity * apu.mol**base_unit[0]
        elif base_unit[1] == cd:
            ap_quantity = ap_quantity * apu.cd**base_unit[0]

    return ap_quantity
