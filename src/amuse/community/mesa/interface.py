"""
MESA chooser
"""


def MESAInterface(version="15140", **options):
    if str(version) == "2208":
        from amuse.community.mesa_r2208.interface import MESAInterface as MESAInterface_2208
        return MESAInterface_2208(**options)
    if str(version) == "15140":
        from amuse.community.mesa_r15140.interface import MESAInterface as MESAInterface_15140
        return MESAInterface_15140(**options)
    raise AttributeError(
        "This version of MESA is not (yet) supported by AMUSE"
    )


def MESA(version="15140", **options):
    if str(version) == "2208":
        from amuse.community.mesa_r2208.interface import MESA as MESA_2208
        return MESA_2208(**options)
    if str(version) == "15140":
        from amuse.community.mesa_r15140.interface import MESA as MESA_15140
        return MESA_15140(**options)
    raise AttributeError(
        "This version of MESA is not (yet) supported by AMUSE"
    )


Mesa = MESA
