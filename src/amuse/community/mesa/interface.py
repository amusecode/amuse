"""
MESA chooser
"""

MESA_VERSIONS_AVAILABLE = []
try:
    from amuse.community.mesa_r2208 import \
        MesaInterface as MesaInterface_2208
    from amuse.community.mesa_r2208 import Mesa as Mesa_2208
    MESA_VERSIONS_AVAILABLE.append("2208")
except ImportError:
    pass

try:
    from amuse.community.mesa_r15140 import \
        MesaInterface as MesaInterface_15140
    from amuse.community.mesa_r15140 import Mesa as Mesa_15140
    MESA_VERSIONS_AVAILABLE.append("15140")
except ImportError:
    pass

if len(MESA_VERSIONS_AVAILABLE) == 0:
    raise ImportError("No MESA versions are available")

def MesaInterface(version="15140", **options):
    if str(version) == "2208":
        return MesaInterface_2208(**options)
    if str(version) == "15140":
        return MesaInterface_15140(**options)
    raise AttributeError(
        "This version of MESA is not (yet) supported by AMUSE"
    )

def Mesa(version="15140", **options):
    if str(version) == "2208":
        return Mesa_2208(**options)
    if str(version) == "15140":
        return Mesa_15140(**options)
    raise AttributeError(
        "This version of MESA is not (yet) supported by AMUSE"
    )

MESA = Mesa
