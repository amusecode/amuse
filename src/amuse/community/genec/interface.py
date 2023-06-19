"""
Interface for GENEC.
"""

import numpy as np
from amuse.datamodel import ParticlesWithFilteredAttributes
from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function, remote_function
from amuse.community import LiteratureReferencesMixIn
from amuse.community import CodeWithDataDirectories
from amuse.community import InCodeComponentImplementation
from amuse.community import NO_UNIT, ERROR_CODE
from amuse.community.interface import common
from amuse.community.interface.se import StellarEvolution
from amuse.community.interface.se import StellarEvolutionInterface
from amuse.community.interface.se import InternalStellarStructure
from amuse.community.interface.se import InternalStellarStructureInterface
from amuse.community.interface.stopping_conditions import (
    StoppingConditionInterface, StoppingConditions
)
from amuse.units import units

SPECIES_NAMES = {
    'h': 1,
    'he3': 2,
    'he': 3,
    'c12': 4,
    'c13': 5,
    'n14': 6,
    'n15': 7,
    'o16': 8,
    'o17': 9,
    'o18': 10,
    'ne20': 11,
    'ne22': 12,
    'mg24': 13,
    'mg25': 14,
    'mg26': 15,
    'c14': 16,
    'f18': 17,
    'f19': 18,
    'ne21': 19,
    'na23': 20,
    'al26': 21,
    'al27': 22,
    'si28': 23,
    'neut': 24,
    'prot': 25,
    'bid': 26,
    'bid1': 27,
}

GENEC_AMUSE_SPECIFIC = {
    'modell': ['int32', '', ''],
    'veryFirst': ['bool', '', ''],
}

# Parameters (but individual to each star)
GENEC_STAR_CHARACTERISTICS = {
    # 'GENEC name': [dtype, unit, description, AMUSE name (optional)]
    'initialised': ['bool', '', "True if the star is an intialised model"],
    'star_name': ['string', '', "Name of the star"],
    'nwmd': ['int32', '', "model number", "step"],
    'nwseq': [
        'int32', '', "number of the first model in the time-step series"
    ],
    'modanf': ['int32', '', ".b file number"],
    'nzmod': ['int32', '', "number of models in a run"],
    'end_at_phase': ['int32', '', "Stop if this phase is reached"],
    'end_at_model': ['int32', '', "Stop if this model number is reached"],
}

GENEC_STAR_PHYSICS = {
    'irot': ['int32', '', "rotating if set to 1"],
    'isol': ['int32', '', "solid rotation if set to 1"],
    'imagn': ['int32', '', "internal magnetic fields (0=none, 1=included)"],
    'ialflu': ['int32', '', "Ne-Na and Mg-Al networks if set to 1"],
    'ianiso': ['int32', '', "wind anisotropy if set to 1", "anisotropic_wind"],
    'ipop3': ['int32', '', "Z=0 models if set to 1"],
    'ibasnet': ['int32', '', "extended nuclear network if set to 1"],
    'phase': ['int32', '', "fusion phases"],
    'var_rates': [
        'bool', '',
        "allows to use different reaction rate files if set to True"
    ],
    'bintide': [
        'bool',
        '',
        "tidal interaction in binaries if set to True",
        "binary_enable_tides",
    ],
    'binm2': [
        'float64',
        'MSun',
        "mass of the companion",
        "binary_companion_mass",
    ],
    'periodini': [
        'float64',
        'day',
        "initial period of the binary",
        "binary_initial_period",
    ],
    'const_per': ['bool', '', "keep constant period if True"],
    'iprezams': ['int32', '', ""],
}

GENEC_STAR_COMPOSITION = {
    'initial_metallicity': ['float64', '', "initial metallicity of the model"],
    'zsol': ['float64', '', "reference solar metallicity"],
    'z': ['float64', '', "abundance of the neglected isotopes"],
    'iopac': ['int32', '', "choice of the opacity table if ikappa = 5"],
    'ikappa': ['int32', '', "opacity choice"],
}

GENEC_STAR_ROTATION = {
    'idiff': [
        'int32', '',
        "computation of the diffusion of Omega and chemicals if set to 1"
    ],
    'iadvec': [
        'int32', '',
        "advecto-diffusive version of the transport of Omega if set to 1"
    ],
    'istati': [
        'int32', '',
        "only local conservation of angular momentum if set to 1"
    ],
    'icoeff': [
        'int32', '', "prescription to be used for the diffusion coefficients"],
    'fenerg': [
        'float64', '', "fraction of the shear energy used for the mixing"],
    'richac': ['float64', '', "critical Richardson number value"],
    'igamma': [
        'int32', '', "treatment of the shear according to M&M96 if set to 1"
    ],
    'frein': ['float64', '', "magnetic braking if ≠ 0"],
    'K_Kawaler': [
        'float64', '', "K parameter in Kawaler1998 magnetic braking law"
    ],
    'Omega_saturation': [
        'float64', '', "Ωsat in Kawaler magnetic braking law"
    ],
    'rapcrilim': [
        'float64', '',
        "maximum Ωcrit ratio before the onset of mechanical mass loss"
    ],
    'zams_velocity': [
        'float64',
        '',
        "chosen velocity on the ZAMS (Veq if > 1.0, V/Vcrit if 0<zams_velocity<1, "
        "Ω/Ωcrit if -1<zams_velocity<0",
        'zams_velocity'
    ],
    'xfom': ['float64', '', "multiplying factor for surface Ω"],
    'omega': ['float64', '', "surface Ω"],
    'xdial': ['float64', '', ""],
    'idialo': ['int32', '', ""],
    'idialu': ['int32', '', ""],
    'Add_Flux': [
        'bool', '',
        "improves the angular momentum conservation treatment if set to True"
    ],
    'diff_only': [
        'bool', '',
        "applies the tidal mixing only in timesteps where we use diffusion "
        "(and not advection) if set to True"
    ],
    'B_initial': [
        'float64', '',
        "if >0, switches on the wind quenching, and evolves the magnetic "
        "field strength according to the conservation of magnetic flux"
    ],
    'add_diff': [
        'float64', '',
        "additional viscosity value, used to increase the core-envelope "
        "coupling"
    ],
    'n_mag': ['int32', '', "type of treatment for magnetic fields"],
    'alpha_F': [
        'float64', '',
        "α parameter in Fuller+ 2019. Values higher than 1 lower the "
        "threshold to trigger the instability (Qmin) and increase the "
        "magnetic viscosity (increased transport)"
    ],
    'nsmooth': [
        'int32', '',
        "number of layers used for smoothing the Ω gradient (used by "
        "Mag_diff_general). Default value is nsmooth=1, recommended "
        "value for Fuller+ 2019 is nsmooth=5"
    ],
    'qminsmooth': [
        'bool', '', "mag. instability always considered if set to True"
    ],
}

GENEC_STAR_SURFACE = {
    'imloss': ['int32', '', "choice of the mass-loss prescription"],
    'fmlos': [
        'float64', '',
        "except for IMLOSS=2 or 3, multiplying factor applied to the mass loss"
    ],
    'ifitm': [
        'int32', '',
        "management of the changes of fitm during redward evolution (and/or "
        "blueward evolution after a RSG phase)"
    ],
    'fitm': [
        'float64',
        '',
        "mass included in the interior",
    ],
    'fitmi': [
        'float64', '',
        "max value of FITM to which the star will come back when going back "
        "to the blue",
    ],
    'fitmi_default': [
        'float64', '',
        '',
        '',
    ],
    'deltal': ['float64', '', "triangle size for L at the surface"],
    'deltat': ['float64', '', "triangle size for T at the surface"],
    'nndr': [
        'int32', '',
        "management of the behaviour of (L,Teff) with respect to the triangle"
    ],
    'RSG_Mdot': ['int32', '', "mass-loss recipe used for RSG"],
    'SupraEddMdot': [
        'bool', '',
        "x3 multiplication factor to the wind in case of supra-Eddington "
        "layers if set to True"
    ],
    'Be_mdotfrac': ['float64', '', ""],
    'start_mdot': [
        'float64', '',
        "value of Ω/Ωcrit at which the Be mass loss starts to apply"
    ],
}

GENEC_STAR_CONVECTION = {
    'iledou': ['int32', '', "Ledoux criterion for convection if set to 1"],
    'idifcon': [
        'int32', '',
        "convection treated as a diffusion if set to 1 (used during O-b "
        "and Si-b)"
    ],
    'iover': ['int32', '', "overshooting taken into account if set to 1"],
    'elph': ['float64', '', "mixing length for the external convective zone"],
    'my': [
        'int32', '',
        "integration of the envelope on ρ rather than P if set to 1"
    ],
    'dovhp': ['float64', '', "value of the α overshooting parameter Λ = α HP"],
    'iunder': ['int32', '', "undershooting taken into account if set to 1"],
    'dunder': ['float64', '', "value of the undershooting parameter"],
}

GENEC_STAR_CONVERGENCE = {
    'gkorm': ['float64', '', "accepted deviation on the structure variables"],
    'alph': [
        'float64', '',
        "fraction of the correction applied for the next iteration"
    ],
    'agdr': [
        'float64', '',
        "absolute value of the maximum correction accepted on r, s, P and T"
    ],
    'faktor': [
        'float64', '',
        "parameter allowing to 'compensate' (during the computation only) for "
        "layers with Lr > Ltot, and hence avoid dL/dr < 0"
    ],
    'dgrp': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on P"
    ],
    'dgrl': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on L"
    ],
    'dgry': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on 4He"
    ],
    'dgrc': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on 12C"
    ],
    'dgro': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on 16O"
    ],
    'dgr20': [
        'float64', '',
        "relative variations (from one layer to the other) accepted on ?"
    ],
    'nbchx': [
        'int32', '',
        "iteration number for the computation of the change in chemical "
        "composition"
    ],
    'nrband': [
        'int32', '',
        "iteration number for the chemicals between model n and n+1"],
}

GENEC_STAR_TIME = {
    'xcn': [
        'float64', '',
        "multiplying factor applied on the time step for the next run"
    ],
    'islow': [
        'int32', '',
        "slow version of the program if not 0 by modification of the ideal "
        "nuclear time step ratxcn"
    ],
    'icncst': ['int32', '', "constant time step (equivalent to xcn=1.0)"],
    'tauH_fit': [
        'int32', '',
        "used to set the maximal timestep in case of critical velocity, as a "
        "fraction of the MS lifetime"
    ],
}

GENEC_STAR_VARIOUS = {
    'display_plot': [
        'bool', '', "display of the pgplot window if set to True"
    ],
    'iauto': [
        'int32', '',
        "management of the parameters change through the different phases of "
        "the evolution"
    ],
    'iprn': [
        'int32', '',
        "the full structure is printed every iprn models. If iprn > nzmod, "
        "only one structure will be printed, the one corresponding to model "
        "nwseq"
    ],
    'iout': [
        'int32', '',
        "number of layers used for the smoothing of the diffusion gradient "
        "at the border of the convective zones"
    ],
    'itmin': ['int32', '', "minimal number of iterations in henyey"],
    'xyfiles': ['bool', '', ""],
    'idebug': [
        'int32', '',
        "addition of some terminal printings useful for the debugging, with "
        "different possible levels. All values > 0 set verbose=True"
    ],
    'itests': ['int32', '', "use of a flag to test some pieces of code"],
    'verbose': [
        'bool', '',
        "increases the level of printings on the terminal and the .l file"
    ],
    'stop_deg': [
        'bool', '',
        "automatically stops a computation when Tc becomes degenerate"
    ],
    'n_snap': ['int32', '', "number of steps between snapshots [0]"],
}

# Stellar properties (but global for the star)
# i.e.: b file (not shells)
GENEC_STAR_PROPERTIES = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (empty = not used)]
    'gms': ['float64', 'MSun', "total mass", "mass"],
    'alter': ['float64', 'yr', "stellar age", "age"],
    'gls': ['float64', 'LSun', "stellar luminosity", "luminosity"],
    'teff': ['float64', 'K', "effective temperature", "temperature"],
    'glsv': ['float64', 'LSun', "previous luminosity"],
    'teffv': ['float64', 'K', "previous temperature"],
    'dzeitj': ['float64', 'yr', "time step (yr)", ""],
    'dzeit': ['float64', 's', "time step (s)", ""],
    'dzeitv': ['float64', 's', "previous time step", ""],
    'xmini': ['float64', 'MSun', "initial mass", "initial_mass"],
    'ab': ['float64', '', "binary separation", ""],
    'dm_lost': ['float64', 'MSun', "total mass lost", ""],
    'm': ['int32', '', "number of zones", "n_zones"],
    'summas': ['float64', 'MSun', "total mass", ""],  # = xmini
    # then a bunch of layered stuff, then
    "dk": ['float64', '', '', ''],
    "rlp": ['float64', '', '', ''],
    "rlt": ['float64', '', '', ''],
    "rlc": ['float64', '', '', ''],
    "rrp": ['float64', '', '', ''],
    "rrt": ['float64', '', '', ''],
    "rrc": ['float64', '', '', ''],
    "rtp": ['float64', '', '', ''],
    "rtt": ['float64', '', '', ''],
    "rtc": ['float64', '', '', ''],
    "tdiff": ['float64', '', '', ''],
    "suminenv": ['float64', '', '', ''],
    "xltotbeg": ['float64', '', '', ''],
    "dlelexprev": ['float64', '', '', ''],
    "radius": ['float64', 'RSun', '', ''],
    "zams_radius": ['float64', '', '', ''],
    'mbelx': ['int32', '', "number of extra elements"],
    'xtefflast': ['float64', '', ""],
    'xllast': ['float64', '', ""],
    'xrholast': ['float64', '', ""],
    'xclast': ['float64', '', ""],
    'xtclast': ['float64', '', ""],
    'inum': ['int32', '', ""],
    'nsugi': ['int32', '', ""],
    'period': ['float64', '', ""],
    'r_core': ['float64', '', ""],
    'vna': ['float64', '', "adiabatic gradient in the envelope"],
    'vnr': ['float64', '', "radiative gradient in the envelope"],
}

GENEC_ARRAY_3 = {
    "drl": ['float64', '', '', ''],  # length 3
    "drte": ['float64', '', '', ''],  # length 3
    "drp": ['float64', '', '', ''],  # length 3
    "drt": ['float64', '', '', ''],  # length 3
    "drr": ['float64', '', '', ''],  # length 3
}

GENEC_ARRAY_NPONDCOUCHE = {
    "CorrOmega": ['float64', '', '', ''],
}

# Structural properties (m layers)
GENEC_ZONE = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (optional)
    'q': ['float64', '', "ln(1-Mr/M)"],
    'p': ['float64', '', "ln(pressure)"],
    't': ['float64', '', "ln(temperature)"],
    'r': ['float64', '', "ln(radius)"],
    's': ['float64', '', "ln(internal luminosity + 1)"],
    'x': ['float64', '', "H abundance fraction"],
    'y3': ['float64', '', "He3 abundance fraction"],
    'y': ['float64', '', "He abundance fraction"],
    'xc12': ['float64', '', "C12 abundance fraction"],
    'xc13': ['float64', '', "C13 abundance fraction"],
    'xn14': ['float64', '', "N14 abundance fraction"],
    'xn15': ['float64', '', "N15 abundance fraction"],
    'xo16': ['float64', '', "O16 abundance fraction"],
    'xo17': ['float64', '', "O17 abundance fraction"],
    'xo18': ['float64', '', "O18 abundance fraction"],
    'xne20': ['float64', '', "Ne20 abundance fraction"],
    'xne22': ['float64', '', "Ne22 abundance fraction"],
    'xmg24': ['float64', '', "Mg24 abundance fraction"],
    'xmg25': ['float64', '', "Mg25 abundance fraction"],
    'xmg26': ['float64', '', "Mg26 abundance fraction"],
    'xf19': ['float64', '', "F19 abundance fraction"],
    'xne21': ['float64', '', "Ne21 abundance fraction"],
    'xna23': ['float64', '', "Na23 abundance fraction"],
    'xal27': ['float64', '', "Al27 abundance fraction"],
    'xsi28': ['float64', '', "Si28 abundance fraction"],
    'xc14': ['float64', '', "C14 abundance fraction"],
    'xf18': ['float64', '', "F18 abundance fraction"],
    'xal26': ['float64', '', "Al26 abundance fraction"],
    'xneut': ['float64', '', "Neutron abundance fraction"],
    'xprot': ['float64', '', "Proton abundance fraction"],
    'omegi': ['float64', '', "Rotation"],
    'xbid': ['float64', '', "Pseudo element abundance fraction"],
    'xbid1': ['float64', '', "Pseudo element 2 abundance fraction"],
    'vp': ['float64', '', "Previous ln(pressure)"],
    'vt': ['float64', '', "Previous ln(temperature)"],
    'vr': ['float64', '', "Previous ln(radius)"],
    'vs': ['float64', '', "Previous ln(internal luminosity + 1)"],
    'vx': ['float64', '', "Previous H abundance fraction"],
    'vy': ['float64', '', "Previous He abundance fraction"],
    'vy3': ['float64', '', "Previous He3 abundance fraction"],
    'vxc12': ['float64', '', "Previous C12 abundance fraction"],
    'vxc13': ['float64', '', "Previous C13 abundance fraction"],
    'vxn14': ['float64', '', "Previous N14 abundance fraction"],
    'vxn15': ['float64', '', "Previous N15 abundance fraction"],
    'vxo16': ['float64', '', "Previous O16 abundance fraction"],
    'vxo17': ['float64', '', "Previous O17 abundance fraction"],
    'vxo18': ['float64', '', "Previous O18 abundance fraction"],
    'vxne20': ['float64', '', "Previous Ne20 abundance fraction"],
    'vxne22': ['float64', '', "Previous Ne22 abundance fraction"],
    'vxmg24': ['float64', '', "Previous Mg24 abundance fraction"],
    'vxmg25': ['float64', '', "Previous Mg25 abundance fraction"],
    'vxmg26': ['float64', '', "Previous Mg26 abundance fraction"],
    'vxf19': ['float64', '', "Previous F19 abundance fraction"],
    'vxne21': ['float64', '', "Previous Ne21 abundance fraction"],
    'vxna23': ['float64', '', "Previous Na23 abundance fraction"],
    'vxal27': ['float64', '', "Previous Al27 abundance fraction"],
    'vxsi28': ['float64', '', "Previous Si28 abundance fraction"],
    'vxc14': ['float64', '', "Previous C14 abundance fraction"],
    'vxf18': ['float64', '', "Previous F18 abundance fraction"],
    'vxal26': ['float64', '', "Previous Al26 abundance fraction"],
    'vxneut': ['float64', '', "Previous Neutron abundance fraction"],
    'vxprot': ['float64', '', "Previous Proton abundance fraction"],
    'vomegi': ['float64', '', "Previous rotation"],
    'vxbid': ['float64', '', "Previous Pseudo element abundance fraction"],
    'vxbid1': ['float64', '', "Previous Pseudo element 2 abundance fraction"],
}

# Structural properties (mbelx elements, m layers)
GENEC_ARRAY_MBELX_M = {
    'abelx': ['float64', '', "extra elements"],
    'vabelx': ['float64', '', "previous extra elements"],
}

GENEC_NETDEF_SCALARS = {
    'xlostneu': ['float64', '', "lost neutrons"],
}
GENEC_NETDEF_ARRAYS = {
    'nbzel': ['int32', '', 'length 8'],
    'nbael': ['int32', '', 'length 8'],
    'abels': ['float64', '', 'length 8'],
}
GENEC_NETALU_ARRAYS = {
    'xnetalu': ['float64', '', 'length 5'],
}

# These are not (re)stored / set, but they are calculated e.g. for plotting
GENEC_ZONE_DERIVED = {
    'eps': ['float64', '', ""],
    'epsy': ['float64', '', ""],
    'eps_c_adv': ['float64', '', ""],
    'eps_ne_adv': ['float64', '', ""],
    'eps_o_adv': ['float64', '', ""],
    'eps_si_adv': ['float64', '', ""],
    'eps_grav': ['float64', '', ""],
    'eps_nu': ['float64', '', ""],
    'nabla_rad': ['float64', '', ""],
    'nabla_ad': ['float64', '', ""],
    'nabla_mu': ['float64', '', ""],
}

ALL_SETTERS = {
    **GENEC_STAR_CHARACTERISTICS,
    **GENEC_STAR_PHYSICS,
    **GENEC_STAR_COMPOSITION,
    **GENEC_STAR_ROTATION,
    **GENEC_STAR_SURFACE,
    **GENEC_STAR_CONVECTION,
    **GENEC_STAR_CONVERGENCE,
    **GENEC_STAR_TIME,
    **GENEC_STAR_VARIOUS,
    **GENEC_STAR_PROPERTIES,
    **GENEC_NETDEF_SCALARS,
    **GENEC_ARRAY_3,
    **GENEC_ARRAY_NPONDCOUCHE,
    **GENEC_ZONE,
    **GENEC_ARRAY_MBELX_M,
    **GENEC_NETDEF_ARRAYS,
    **GENEC_NETALU_ARRAYS,
}

SCALAR_SETTERS = {
    **GENEC_AMUSE_SPECIFIC,
    **GENEC_STAR_CHARACTERISTICS,
    **GENEC_STAR_PHYSICS,
    **GENEC_STAR_COMPOSITION,
    **GENEC_STAR_ROTATION,
    **GENEC_STAR_SURFACE,
    **GENEC_STAR_CONVECTION,
    **GENEC_STAR_CONVERGENCE,
    **GENEC_STAR_TIME,
    **GENEC_STAR_VARIOUS,
    **GENEC_STAR_PROPERTIES,
    **GENEC_NETDEF_SCALARS,
}

VECTOR_SETTERS = {
    **GENEC_ARRAY_3,
    **GENEC_ARRAY_NPONDCOUCHE,
    **GENEC_ZONE,
    **GENEC_ARRAY_MBELX_M,
    **GENEC_NETDEF_ARRAYS,
    **GENEC_NETALU_ARRAYS,
}

ALL_GETTERS = {
    **ALL_SETTERS,
    **GENEC_ZONE_DERIVED,
}


class GenecInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    StellarEvolutionInterface,
    InternalStellarStructureInterface,
    CodeWithDataDirectories,
    # StoppingConditionInterface,
):
    """
    GENEC is the Geneva Stellar Evolution Code

    References:
        .. [#] The Geneva Stellar Evolution Group
    """

    use_modules = ['AmuseInterface', ]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="genec_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)

    @remote_function
    def finalize_stellar_model():
        returns ()

    @remote_function
    def set_genec_path(genec_path='s'):
        returns ()

    @remote_function
    def commit_parameters():
        returns ()

    @legacy_function
    def new_particle():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The initial mass of the star")
        function.addParameter(
            'metallicity', dtype='float64', direction=function.IN,
            default=0.014,
            description="The initial metallicity of the star (default: 0.014)")
        function.addParameter(
            'zams_velocity', dtype='float64', direction=function.IN,
            default=0.0,
            description="The desired ZAMS velocity of the star (default: 0.0)")
        function.addParameter(
            'star_name', dtype='string', direction=function.IN,
            default='AmuseStar', description="The star's name")
        # function.addParameter(
        #     'age_tag', dtype='float64', direction=function.IN,
        #     description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_particle parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function

    @legacy_function
    def new_stellar_model():
        """
        New star from an existing model
        """
        function = LegacyFunctionSpecification()
        # function.can_handle_array = True
        function.must_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        for parameter in SCALAR_SETTERS.items():
            if parameter[1][1] == "":
                unit = NO_UNIT
            else:
                unit = getattr(units, parameter[1][1])
            function.addParameter(
                parameter[0],
                dtype=parameter[1][0],
                unit=unit,
                description=parameter[1][2],
                direction=function.IN,
            )
        for parameter in GENEC_ZONE.items():
            if parameter[1][1] == "":
                unit = NO_UNIT
            else:
                unit = getattr(units, parameter[1][1])
            function.addParameter(
                parameter[0],
                dtype=parameter[1][0],
                unit=unit,
                description=parameter[1][2],
                direction=function.IN,
            )
        # for parameter in GENEC_STAR_STRUCTURE_EXTRA.items():
        #     if parameter[1][1] == "":
        #         unit = NO_UNIT
        #     else:
        #         unit = getattr(units, parameter[1][1])
        #     function.addParameter(
        #         parameter[0],
        #         dtype=parameter[1][0],
        #         unit=unit,
        #         description=parameter[1][2],
        #         direction=function.IN,
        #     )
        function.addParameter('n', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function

    @remote_function(can_handle_array=True)
    def get_firstlast_zone(index_of_the_particle='i'):
        returns (first='i', last='i')

    @remote_function(can_handle_array=True)
    def get_mass_fraction_at_zone(index_of_the_particle='i', zone='i'):
        returns (dq_i='d')

    @remote_function(can_handle_array=True)
    def get_surface_velocity(index_of_the_particle='i'):
        returns (surface_velocity='d')

    @remote_function(can_handle_array=True)
    def get_luminosity_at_zone(index_of_the_particle='i', zone='i'):
        returns (lum_i='d' | units.LSun)

    @remote_function(can_handle_array=True)
    def get_mass_of_species(index_of_the_particle='i', species='i'):
        returns (species_mass='d')

    @remote_function(can_handle_array=True)
    def get_modell(index_of_the_particle='i'):
        returns (modell='int32')

    @remote_function(can_handle_array=True)
    def set_modell(
        index_of_the_particle='i',
        modell='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_veryFirst(index_of_the_particle='i'):
        returns (veryFirst='bool')

    @remote_function(can_handle_array=True)
    def set_veryFirst(
        index_of_the_particle='i',
        veryFirst='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_initialised(index_of_the_particle='i'):
        returns (initialised='bool')

    @remote_function(can_handle_array=True)
    def set_initialised(
        index_of_the_particle='i',
        initialised='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_star_name(index_of_the_particle='i'):
        returns (star_name='string')

    @remote_function(can_handle_array=True)
    def set_star_name(
        index_of_the_particle='i',
        star_name='string',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nwmd(index_of_the_particle='i'):
        returns (nwmd='int32')

    @remote_function(can_handle_array=True)
    def set_nwmd(
        index_of_the_particle='i',
        nwmd='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nwseq(index_of_the_particle='i'):
        returns (nwseq='int32')

    @remote_function(can_handle_array=True)
    def set_nwseq(
        index_of_the_particle='i',
        nwseq='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_modanf(index_of_the_particle='i'):
        returns (modanf='int32')

    @remote_function(can_handle_array=True)
    def set_modanf(
        index_of_the_particle='i',
        modanf='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nzmod(index_of_the_particle='i'):
        returns (nzmod='int32')

    @remote_function(can_handle_array=True)
    def set_nzmod(
        index_of_the_particle='i',
        nzmod='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_end_at_phase(index_of_the_particle='i'):
        returns (end_at_phase='int32')

    @remote_function(can_handle_array=True)
    def set_end_at_phase(
        index_of_the_particle='i',
        end_at_phase='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_end_at_model(index_of_the_particle='i'):
        returns (end_at_model='int32')

    @remote_function(can_handle_array=True)
    def set_end_at_model(
        index_of_the_particle='i',
        end_at_model='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_irot(index_of_the_particle='i'):
        returns (irot='int32')

    @remote_function(can_handle_array=True)
    def set_irot(
        index_of_the_particle='i',
        irot='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_isol(index_of_the_particle='i'):
        returns (isol='int32')

    @remote_function(can_handle_array=True)
    def set_isol(
        index_of_the_particle='i',
        isol='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_imagn(index_of_the_particle='i'):
        returns (imagn='int32')

    @remote_function(can_handle_array=True)
    def set_imagn(
        index_of_the_particle='i',
        imagn='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ialflu(index_of_the_particle='i'):
        returns (ialflu='int32')

    @remote_function(can_handle_array=True)
    def set_ialflu(
        index_of_the_particle='i',
        ialflu='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ianiso(index_of_the_particle='i'):
        returns (ianiso='int32')

    @remote_function(can_handle_array=True)
    def set_ianiso(
        index_of_the_particle='i',
        ianiso='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ipop3(index_of_the_particle='i'):
        returns (ipop3='int32')

    @remote_function(can_handle_array=True)
    def set_ipop3(
        index_of_the_particle='i',
        ipop3='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ibasnet(index_of_the_particle='i'):
        returns (ibasnet='int32')

    @remote_function(can_handle_array=True)
    def set_ibasnet(
        index_of_the_particle='i',
        ibasnet='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_phase(index_of_the_particle='i'):
        returns (phase='int32')

    @remote_function(can_handle_array=True)
    def set_phase(
        index_of_the_particle='i',
        phase='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_var_rates(index_of_the_particle='i'):
        returns (var_rates='bool')

    @remote_function(can_handle_array=True)
    def set_var_rates(
        index_of_the_particle='i',
        var_rates='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_bintide(index_of_the_particle='i'):
        returns (bintide='bool')

    @remote_function(can_handle_array=True)
    def set_bintide(
        index_of_the_particle='i',
        bintide='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_binm2(index_of_the_particle='i'):
        returns (binm2='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_binm2(
        index_of_the_particle='i',
        binm2='float64' | units.MSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_periodini(index_of_the_particle='i'):
        returns (periodini='float64' | units.day)

    @remote_function(can_handle_array=True)
    def set_periodini(
        index_of_the_particle='i',
        periodini='float64' | units.day,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_const_per(index_of_the_particle='i'):
        returns (const_per='bool')

    @remote_function(can_handle_array=True)
    def set_const_per(
        index_of_the_particle='i',
        const_per='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iprezams(index_of_the_particle='i'):
        returns (iprezams='int32')

    @remote_function(can_handle_array=True)
    def set_iprezams(
        index_of_the_particle='i',
        iprezams='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_initial_metallicity(index_of_the_particle='i'):
        returns (initial_metallicity='float64')

    @remote_function(can_handle_array=True)
    def set_initial_metallicity(
        index_of_the_particle='i',
        initial_metallicity='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_zsol(index_of_the_particle='i'):
        returns (zsol='float64')

    @remote_function(can_handle_array=True)
    def set_zsol(
        index_of_the_particle='i',
        zsol='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_z(index_of_the_particle='i'):
        returns (z='float64')

    @remote_function(can_handle_array=True)
    def set_z(
        index_of_the_particle='i',
        z='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iopac(index_of_the_particle='i'):
        returns (iopac='int32')

    @remote_function(can_handle_array=True)
    def set_iopac(
        index_of_the_particle='i',
        iopac='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ikappa(index_of_the_particle='i'):
        returns (ikappa='int32')

    @remote_function(can_handle_array=True)
    def set_ikappa(
        index_of_the_particle='i',
        ikappa='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idiff(index_of_the_particle='i'):
        returns (idiff='int32')

    @remote_function(can_handle_array=True)
    def set_idiff(
        index_of_the_particle='i',
        idiff='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iadvec(index_of_the_particle='i'):
        returns (iadvec='int32')

    @remote_function(can_handle_array=True)
    def set_iadvec(
        index_of_the_particle='i',
        iadvec='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_istati(index_of_the_particle='i'):
        returns (istati='int32')

    @remote_function(can_handle_array=True)
    def set_istati(
        index_of_the_particle='i',
        istati='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_icoeff(index_of_the_particle='i'):
        returns (icoeff='int32')

    @remote_function(can_handle_array=True)
    def set_icoeff(
        index_of_the_particle='i',
        icoeff='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fenerg(index_of_the_particle='i'):
        returns (fenerg='float64')

    @remote_function(can_handle_array=True)
    def set_fenerg(
        index_of_the_particle='i',
        fenerg='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_richac(index_of_the_particle='i'):
        returns (richac='float64')

    @remote_function(can_handle_array=True)
    def set_richac(
        index_of_the_particle='i',
        richac='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_igamma(index_of_the_particle='i'):
        returns (igamma='int32')

    @remote_function(can_handle_array=True)
    def set_igamma(
        index_of_the_particle='i',
        igamma='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_frein(index_of_the_particle='i'):
        returns (frein='float64')

    @remote_function(can_handle_array=True)
    def set_frein(
        index_of_the_particle='i',
        frein='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_K_Kawaler(index_of_the_particle='i'):
        returns (K_Kawaler='float64')

    @remote_function(can_handle_array=True)
    def set_K_Kawaler(
        index_of_the_particle='i',
        K_Kawaler='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Omega_saturation(index_of_the_particle='i'):
        returns (Omega_saturation='float64')

    @remote_function(can_handle_array=True)
    def set_Omega_saturation(
        index_of_the_particle='i',
        Omega_saturation='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rapcrilim(index_of_the_particle='i'):
        returns (rapcrilim='float64')

    @remote_function(can_handle_array=True)
    def set_rapcrilim(
        index_of_the_particle='i',
        rapcrilim='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_zams_velocity(index_of_the_particle='i'):
        returns (zams_velocity='float64')

    @remote_function(can_handle_array=True)
    def set_zams_velocity(
        index_of_the_particle='i',
        zams_velocity='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xfom(index_of_the_particle='i'):
        returns (xfom='float64')

    @remote_function(can_handle_array=True)
    def set_xfom(
        index_of_the_particle='i',
        xfom='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_omega(index_of_the_particle='i'):
        returns (omega='float64')

    @remote_function(can_handle_array=True)
    def set_omega(
        index_of_the_particle='i',
        omega='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xdial(index_of_the_particle='i'):
        returns (xdial='float64')

    @remote_function(can_handle_array=True)
    def set_xdial(
        index_of_the_particle='i',
        xdial='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idialo(index_of_the_particle='i'):
        returns (idialo='int32')

    @remote_function(can_handle_array=True)
    def set_idialo(
        index_of_the_particle='i',
        idialo='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idialu(index_of_the_particle='i'):
        returns (idialu='int32')

    @remote_function(can_handle_array=True)
    def set_idialu(
        index_of_the_particle='i',
        idialu='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Add_Flux(index_of_the_particle='i'):
        returns (Add_Flux='bool')

    @remote_function(can_handle_array=True)
    def set_Add_Flux(
        index_of_the_particle='i',
        Add_Flux='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_diff_only(index_of_the_particle='i'):
        returns (diff_only='bool')

    @remote_function(can_handle_array=True)
    def set_diff_only(
        index_of_the_particle='i',
        diff_only='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_B_initial(index_of_the_particle='i'):
        returns (B_initial='float64')

    @remote_function(can_handle_array=True)
    def set_B_initial(
        index_of_the_particle='i',
        B_initial='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_add_diff(index_of_the_particle='i'):
        returns (add_diff='float64')

    @remote_function(can_handle_array=True)
    def set_add_diff(
        index_of_the_particle='i',
        add_diff='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_n_mag(index_of_the_particle='i'):
        returns (n_mag='int32')

    @remote_function(can_handle_array=True)
    def set_n_mag(
        index_of_the_particle='i',
        n_mag='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alpha_F(index_of_the_particle='i'):
        returns (alpha_F='float64')

    @remote_function(can_handle_array=True)
    def set_alpha_F(
        index_of_the_particle='i',
        alpha_F='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nsmooth(index_of_the_particle='i'):
        returns (nsmooth='int32')

    @remote_function(can_handle_array=True)
    def set_nsmooth(
        index_of_the_particle='i',
        nsmooth='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_qminsmooth(index_of_the_particle='i'):
        returns (qminsmooth='bool')

    @remote_function(can_handle_array=True)
    def set_qminsmooth(
        index_of_the_particle='i',
        qminsmooth='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_imloss(index_of_the_particle='i'):
        returns (imloss='int32')

    @remote_function(can_handle_array=True)
    def set_imloss(
        index_of_the_particle='i',
        imloss='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fmlos(index_of_the_particle='i'):
        returns (fmlos='float64')

    @remote_function(can_handle_array=True)
    def set_fmlos(
        index_of_the_particle='i',
        fmlos='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ifitm(index_of_the_particle='i'):
        returns (ifitm='int32')

    @remote_function(can_handle_array=True)
    def set_ifitm(
        index_of_the_particle='i',
        ifitm='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fitm(index_of_the_particle='i'):
        returns (fitm='float64')

    @remote_function(can_handle_array=True)
    def set_fitm(
        index_of_the_particle='i',
        fitm='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fitmi(index_of_the_particle='i'):
        returns (fitmi='float64')

    @remote_function(can_handle_array=True)
    def set_fitmi(
        index_of_the_particle='i',
        fitmi='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fitmi_default(index_of_the_particle='i'):
        returns (fitmi_default='float64')

    @remote_function(can_handle_array=True)
    def set_fitmi_default(
        index_of_the_particle='i',
        fitmi_default='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_deltal(index_of_the_particle='i'):
        returns (deltal='float64')

    @remote_function(can_handle_array=True)
    def set_deltal(
        index_of_the_particle='i',
        deltal='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_deltat(index_of_the_particle='i'):
        returns (deltat='float64')

    @remote_function(can_handle_array=True)
    def set_deltat(
        index_of_the_particle='i',
        deltat='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nndr(index_of_the_particle='i'):
        returns (nndr='int32')

    @remote_function(can_handle_array=True)
    def set_nndr(
        index_of_the_particle='i',
        nndr='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_RSG_Mdot(index_of_the_particle='i'):
        returns (RSG_Mdot='int32')

    @remote_function(can_handle_array=True)
    def set_RSG_Mdot(
        index_of_the_particle='i',
        RSG_Mdot='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_SupraEddMdot(index_of_the_particle='i'):
        returns (SupraEddMdot='bool')

    @remote_function(can_handle_array=True)
    def set_SupraEddMdot(
        index_of_the_particle='i',
        SupraEddMdot='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Be_mdotfrac(index_of_the_particle='i'):
        returns (Be_mdotfrac='float64')

    @remote_function(can_handle_array=True)
    def set_Be_mdotfrac(
        index_of_the_particle='i',
        Be_mdotfrac='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_start_mdot(index_of_the_particle='i'):
        returns (start_mdot='float64')

    @remote_function(can_handle_array=True)
    def set_start_mdot(
        index_of_the_particle='i',
        start_mdot='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iledou(index_of_the_particle='i'):
        returns (iledou='int32')

    @remote_function(can_handle_array=True)
    def set_iledou(
        index_of_the_particle='i',
        iledou='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idifcon(index_of_the_particle='i'):
        returns (idifcon='int32')

    @remote_function(can_handle_array=True)
    def set_idifcon(
        index_of_the_particle='i',
        idifcon='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iover(index_of_the_particle='i'):
        returns (iover='int32')

    @remote_function(can_handle_array=True)
    def set_iover(
        index_of_the_particle='i',
        iover='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_elph(index_of_the_particle='i'):
        returns (elph='float64')

    @remote_function(can_handle_array=True)
    def set_elph(
        index_of_the_particle='i',
        elph='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_my(index_of_the_particle='i'):
        returns (my='int32')

    @remote_function(can_handle_array=True)
    def set_my(
        index_of_the_particle='i',
        my='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dovhp(index_of_the_particle='i'):
        returns (dovhp='float64')

    @remote_function(can_handle_array=True)
    def set_dovhp(
        index_of_the_particle='i',
        dovhp='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iunder(index_of_the_particle='i'):
        returns (iunder='int32')

    @remote_function(can_handle_array=True)
    def set_iunder(
        index_of_the_particle='i',
        iunder='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dunder(index_of_the_particle='i'):
        returns (dunder='float64')

    @remote_function(can_handle_array=True)
    def set_dunder(
        index_of_the_particle='i',
        dunder='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gkorm(index_of_the_particle='i'):
        returns (gkorm='float64')

    @remote_function(can_handle_array=True)
    def set_gkorm(
        index_of_the_particle='i',
        gkorm='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alph(index_of_the_particle='i'):
        returns (alph='float64')

    @remote_function(can_handle_array=True)
    def set_alph(
        index_of_the_particle='i',
        alph='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_agdr(index_of_the_particle='i'):
        returns (agdr='float64')

    @remote_function(can_handle_array=True)
    def set_agdr(
        index_of_the_particle='i',
        agdr='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_faktor(index_of_the_particle='i'):
        returns (faktor='float64')

    @remote_function(can_handle_array=True)
    def set_faktor(
        index_of_the_particle='i',
        faktor='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrp(index_of_the_particle='i'):
        returns (dgrp='float64')

    @remote_function(can_handle_array=True)
    def set_dgrp(
        index_of_the_particle='i',
        dgrp='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrl(index_of_the_particle='i'):
        returns (dgrl='float64')

    @remote_function(can_handle_array=True)
    def set_dgrl(
        index_of_the_particle='i',
        dgrl='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgry(index_of_the_particle='i'):
        returns (dgry='float64')

    @remote_function(can_handle_array=True)
    def set_dgry(
        index_of_the_particle='i',
        dgry='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrc(index_of_the_particle='i'):
        returns (dgrc='float64')

    @remote_function(can_handle_array=True)
    def set_dgrc(
        index_of_the_particle='i',
        dgrc='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgro(index_of_the_particle='i'):
        returns (dgro='float64')

    @remote_function(can_handle_array=True)
    def set_dgro(
        index_of_the_particle='i',
        dgro='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgr20(index_of_the_particle='i'):
        returns (dgr20='float64')

    @remote_function(can_handle_array=True)
    def set_dgr20(
        index_of_the_particle='i',
        dgr20='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nbchx(index_of_the_particle='i'):
        returns (nbchx='int32')

    @remote_function(can_handle_array=True)
    def set_nbchx(
        index_of_the_particle='i',
        nbchx='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nrband(index_of_the_particle='i'):
        returns (nrband='int32')

    @remote_function(can_handle_array=True)
    def set_nrband(
        index_of_the_particle='i',
        nrband='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xcn(index_of_the_particle='i'):
        returns (xcn='float64')

    @remote_function(can_handle_array=True)
    def set_xcn(
        index_of_the_particle='i',
        xcn='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_islow(index_of_the_particle='i'):
        returns (islow='int32')

    @remote_function(can_handle_array=True)
    def set_islow(
        index_of_the_particle='i',
        islow='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_icncst(index_of_the_particle='i'):
        returns (icncst='int32')

    @remote_function(can_handle_array=True)
    def set_icncst(
        index_of_the_particle='i',
        icncst='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_tauH_fit(index_of_the_particle='i'):
        returns (tauH_fit='int32')

    @remote_function(can_handle_array=True)
    def set_tauH_fit(
        index_of_the_particle='i',
        tauH_fit='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_display_plot(index_of_the_particle='i'):
        returns (display_plot='bool')

    @remote_function(can_handle_array=True)
    def set_display_plot(
        index_of_the_particle='i',
        display_plot='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iauto(index_of_the_particle='i'):
        returns (iauto='int32')

    @remote_function(can_handle_array=True)
    def set_iauto(
        index_of_the_particle='i',
        iauto='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iprn(index_of_the_particle='i'):
        returns (iprn='int32')

    @remote_function(can_handle_array=True)
    def set_iprn(
        index_of_the_particle='i',
        iprn='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iout(index_of_the_particle='i'):
        returns (iout='int32')

    @remote_function(can_handle_array=True)
    def set_iout(
        index_of_the_particle='i',
        iout='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_itmin(index_of_the_particle='i'):
        returns (itmin='int32')

    @remote_function(can_handle_array=True)
    def set_itmin(
        index_of_the_particle='i',
        itmin='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xyfiles(index_of_the_particle='i'):
        returns (xyfiles='bool')

    @remote_function(can_handle_array=True)
    def set_xyfiles(
        index_of_the_particle='i',
        xyfiles='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idebug(index_of_the_particle='i'):
        returns (idebug='int32')

    @remote_function(can_handle_array=True)
    def set_idebug(
        index_of_the_particle='i',
        idebug='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_itests(index_of_the_particle='i'):
        returns (itests='int32')

    @remote_function(can_handle_array=True)
    def set_itests(
        index_of_the_particle='i',
        itests='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_verbose(index_of_the_particle='i'):
        returns (verbose='bool')

    @remote_function(can_handle_array=True)
    def set_verbose(
        index_of_the_particle='i',
        verbose='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_stop_deg(index_of_the_particle='i'):
        returns (stop_deg='bool')

    @remote_function(can_handle_array=True)
    def set_stop_deg(
        index_of_the_particle='i',
        stop_deg='bool',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_n_snap(index_of_the_particle='i'):
        returns (n_snap='int32')

    @remote_function(can_handle_array=True)
    def set_n_snap(
        index_of_the_particle='i',
        n_snap='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gms(index_of_the_particle='i'):
        returns (gms='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_gms(
        index_of_the_particle='i',
        gms='float64' | units.MSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alter(index_of_the_particle='i'):
        returns (alter='float64' | units.julianyr)

    @remote_function(can_handle_array=True)
    def set_alter(
        index_of_the_particle='i',
        alter='float64' | units.julianyr,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gls(index_of_the_particle='i'):
        returns (gls='float64' | units.LSun)

    @remote_function(can_handle_array=True)
    def set_gls(
        index_of_the_particle='i',
        gls='float64' | units.LSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_teff(index_of_the_particle='i'):
        returns (teff='float64' | units.K)

    @remote_function(can_handle_array=True)
    def set_teff(
        index_of_the_particle='i',
        teff='float64' | units.K,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_glsv(index_of_the_particle='i'):
        returns (glsv='float64' | units.LSun)

    @remote_function(can_handle_array=True)
    def set_glsv(
        index_of_the_particle='i',
        glsv='float64' | units.LSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_teffv(index_of_the_particle='i'):
        returns (teffv='float64' | units.K)

    @remote_function(can_handle_array=True)
    def set_teffv(
        index_of_the_particle='i',
        teffv='float64' | units.K,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeitj(index_of_the_particle='i'):
        returns (dzeitj='float64' | units.julianyr)

    @remote_function(can_handle_array=True)
    def set_dzeitj(
        index_of_the_particle='i',
        dzeitj='float64' | units.julianyr,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeit(index_of_the_particle='i'):
        returns (dzeit='float64' | units.s)

    @remote_function(can_handle_array=True)
    def set_dzeit(
        index_of_the_particle='i',
        dzeit='float64' | units.s,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_time_step(
        index_of_the_particle='i', 
    ):
        returns (time_step='float64' | units.s)

    @remote_function(can_handle_array=True)
    def set_time_step(
        index_of_the_particle='i',
        time_step='float64' | units.s,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeitv(index_of_the_particle='i'):
        returns (dzeitv='float64' | units.s)

    @remote_function(can_handle_array=True)
    def set_dzeitv(
        index_of_the_particle='i',
        dzeitv='float64' | units.s,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xmini(index_of_the_particle='i'):
        returns (xmini='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_xmini(
        index_of_the_particle='i',
        xmini='float64' | units.MSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ab(index_of_the_particle='i'):
        returns (ab='float64')

    @remote_function(can_handle_array=True)
    def set_ab(
        index_of_the_particle='i',
        ab='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dm_lost(index_of_the_particle='i'):
        returns (dm_lost='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_dm_lost(
        index_of_the_particle='i',
        dm_lost='float64' | units.MSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_m(index_of_the_particle='i'):
        returns (m='int32')

    @remote_function(can_handle_array=True)
    def set_m(
        index_of_the_particle='i',
        m='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_summas(index_of_the_particle='i'):
        returns (summas='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_summas(
        index_of_the_particle='i',
        summas='float64' | units.MSun,
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dk(index_of_the_particle='i'):
        returns (dk='float64')

    @remote_function(can_handle_array=True)
    def set_dk(
        index_of_the_particle='i',
        dk='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rlp(index_of_the_particle='i'):
        returns (rlp='float64')

    @remote_function(can_handle_array=True)
    def set_rlp(
        index_of_the_particle='i',
        rlp='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rlt(index_of_the_particle='i'):
        returns (rlt='float64')

    @remote_function(can_handle_array=True)
    def set_rlt(
        index_of_the_particle='i',
        rlt='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rlc(index_of_the_particle='i'):
        returns (rlc='float64')

    @remote_function(can_handle_array=True)
    def set_rlc(
        index_of_the_particle='i',
        rlc='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rrp(index_of_the_particle='i'):
        returns (rrp='float64')

    @remote_function(can_handle_array=True)
    def set_rrp(
        index_of_the_particle='i',
        rrp='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rrt(index_of_the_particle='i'):
        returns (rrt='float64')

    @remote_function(can_handle_array=True)
    def set_rrt(
        index_of_the_particle='i',
        rrt='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rrc(index_of_the_particle='i'):
        returns (rrc='float64')

    @remote_function(can_handle_array=True)
    def set_rrc(
        index_of_the_particle='i',
        rrc='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rtp(index_of_the_particle='i'):
        returns (rtp='float64')

    @remote_function(can_handle_array=True)
    def set_rtp(
        index_of_the_particle='i',
        rtp='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rtt(index_of_the_particle='i'):
        returns (rtt='float64')

    @remote_function(can_handle_array=True)
    def set_rtt(
        index_of_the_particle='i',
        rtt='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rtc(index_of_the_particle='i'):
        returns (rtc='float64')

    @remote_function(can_handle_array=True)
    def set_rtc(
        index_of_the_particle='i',
        rtc='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_tdiff(index_of_the_particle='i'):
        returns (tdiff='float64')

    @remote_function(can_handle_array=True)
    def set_tdiff(
        index_of_the_particle='i',
        tdiff='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_suminenv(index_of_the_particle='i'):
        returns (suminenv='float64')

    @remote_function(can_handle_array=True)
    def set_suminenv(
        index_of_the_particle='i',
        suminenv='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xltotbeg(index_of_the_particle='i'):
        returns (xltotbeg='float64')

    @remote_function(can_handle_array=True)
    def set_xltotbeg(
        index_of_the_particle='i',
        xltotbeg='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dlelexprev(index_of_the_particle='i'):
        returns (dlelexprev='float64')

    @remote_function(can_handle_array=True)
    def set_dlelexprev(
        index_of_the_particle='i',
        dlelexprev='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_zams_radius(index_of_the_particle='i'):
        returns (zams_radius='float64')

    @remote_function(can_handle_array=True)
    def set_zams_radius(
        index_of_the_particle='i',
        zams_radius='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_mbelx(index_of_the_particle='i'):
        returns (mbelx='int32')

    @remote_function(can_handle_array=True)
    def set_mbelx(
        index_of_the_particle='i',
        mbelx='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xtefflast(index_of_the_particle='i'):
        returns (xtefflast='float64')

    @remote_function(can_handle_array=True)
    def set_xtefflast(
        index_of_the_particle='i',
        xtefflast='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xllast(index_of_the_particle='i'):
        returns (xllast='float64')

    @remote_function(can_handle_array=True)
    def set_xllast(
        index_of_the_particle='i',
        xllast='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xrholast(index_of_the_particle='i'):
        returns (xrholast='float64')

    @remote_function(can_handle_array=True)
    def set_xrholast(
        index_of_the_particle='i',
        xrholast='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xclast(index_of_the_particle='i'):
        returns (xclast='float64')

    @remote_function(can_handle_array=True)
    def set_xclast(
        index_of_the_particle='i',
        xclast='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xtclast(index_of_the_particle='i'):
        returns (xtclast='float64')

    @remote_function(can_handle_array=True)
    def set_xtclast(
        index_of_the_particle='i',
        xtclast='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_inum(index_of_the_particle='i'):
        returns (inum='int32')

    @remote_function(can_handle_array=True)
    def set_inum(
        index_of_the_particle='i',
        inum='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nsugi(index_of_the_particle='i'):
        returns (nsugi='int32')

    @remote_function(can_handle_array=True)
    def set_nsugi(
        index_of_the_particle='i',
        nsugi='int32',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_period(index_of_the_particle='i'):
        returns (period='float64')

    @remote_function(can_handle_array=True)
    def set_period(
        index_of_the_particle='i',
        period='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_r_core(index_of_the_particle='i'):
        returns (r_core='float64')

    @remote_function(can_handle_array=True)
    def set_r_core(
        index_of_the_particle='i',
        r_core='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vna(index_of_the_particle='i'):
        returns (vna='float64')

    @remote_function(can_handle_array=True)
    def set_vna(
        index_of_the_particle='i',
        vna='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vnr(index_of_the_particle='i'):
        returns (vnr='float64')

    @remote_function(can_handle_array=True)
    def set_vnr(
        index_of_the_particle='i',
        vnr='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xlostneu(index_of_the_particle='i'):
        returns (xlostneu='float64')

    @remote_function(can_handle_array=True)
    def set_xlostneu(
        index_of_the_particle='i',
        xlostneu='float64',
    ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_q_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (q='float64')

    @remote_function(can_handle_array=True)
    def set_q_at_zone(
        index_of_the_particle='i', zone='i',
        q='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_p_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (p='float64')

    @remote_function(can_handle_array=True)
    def set_p_at_zone(
        index_of_the_particle='i', zone='i',
        p='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_t_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (t='float64')

    @remote_function(can_handle_array=True)
    def set_t_at_zone(
        index_of_the_particle='i', zone='i',
        t='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_r_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (r='float64')

    @remote_function(can_handle_array=True)
    def set_r_at_zone(
        index_of_the_particle='i', zone='i',
        r='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_s_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (s='float64')

    @remote_function(can_handle_array=True)
    def set_s_at_zone(
        index_of_the_particle='i', zone='i',
        s='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_x_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (x='float64')

    @remote_function(can_handle_array=True)
    def set_x_at_zone(
        index_of_the_particle='i', zone='i',
        x='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_y3_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (y3='float64')

    @remote_function(can_handle_array=True)
    def set_y3_at_zone(
        index_of_the_particle='i', zone='i',
        y3='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_y_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (y='float64')

    @remote_function(can_handle_array=True)
    def set_y_at_zone(
        index_of_the_particle='i', zone='i',
        y='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xc12_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xc12='float64')

    @remote_function(can_handle_array=True)
    def set_xc12_at_zone(
        index_of_the_particle='i', zone='i',
        xc12='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xc13_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xc13='float64')

    @remote_function(can_handle_array=True)
    def set_xc13_at_zone(
        index_of_the_particle='i', zone='i',
        xc13='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xn14_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xn14='float64')

    @remote_function(can_handle_array=True)
    def set_xn14_at_zone(
        index_of_the_particle='i', zone='i',
        xn14='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xn15_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xn15='float64')

    @remote_function(can_handle_array=True)
    def set_xn15_at_zone(
        index_of_the_particle='i', zone='i',
        xn15='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xo16_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xo16='float64')

    @remote_function(can_handle_array=True)
    def set_xo16_at_zone(
        index_of_the_particle='i', zone='i',
        xo16='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xo17_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xo17='float64')

    @remote_function(can_handle_array=True)
    def set_xo17_at_zone(
        index_of_the_particle='i', zone='i',
        xo17='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xo18_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xo18='float64')

    @remote_function(can_handle_array=True)
    def set_xo18_at_zone(
        index_of_the_particle='i', zone='i',
        xo18='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xne20_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xne20='float64')

    @remote_function(can_handle_array=True)
    def set_xne20_at_zone(
        index_of_the_particle='i', zone='i',
        xne20='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xne22_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xne22='float64')

    @remote_function(can_handle_array=True)
    def set_xne22_at_zone(
        index_of_the_particle='i', zone='i',
        xne22='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xmg24_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xmg24='float64')

    @remote_function(can_handle_array=True)
    def set_xmg24_at_zone(
        index_of_the_particle='i', zone='i',
        xmg24='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xmg25_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xmg25='float64')

    @remote_function(can_handle_array=True)
    def set_xmg25_at_zone(
        index_of_the_particle='i', zone='i',
        xmg25='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xmg26_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xmg26='float64')

    @remote_function(can_handle_array=True)
    def set_xmg26_at_zone(
        index_of_the_particle='i', zone='i',
        xmg26='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xf19_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xf19='float64')

    @remote_function(can_handle_array=True)
    def set_xf19_at_zone(
        index_of_the_particle='i', zone='i',
        xf19='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xne21_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xne21='float64')

    @remote_function(can_handle_array=True)
    def set_xne21_at_zone(
        index_of_the_particle='i', zone='i',
        xne21='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xna23_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xna23='float64')

    @remote_function(can_handle_array=True)
    def set_xna23_at_zone(
        index_of_the_particle='i', zone='i',
        xna23='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xal27_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xal27='float64')

    @remote_function(can_handle_array=True)
    def set_xal27_at_zone(
        index_of_the_particle='i', zone='i',
        xal27='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xsi28_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xsi28='float64')

    @remote_function(can_handle_array=True)
    def set_xsi28_at_zone(
        index_of_the_particle='i', zone='i',
        xsi28='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xc14_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xc14='float64')

    @remote_function(can_handle_array=True)
    def set_xc14_at_zone(
        index_of_the_particle='i', zone='i',
        xc14='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xf18_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xf18='float64')

    @remote_function(can_handle_array=True)
    def set_xf18_at_zone(
        index_of_the_particle='i', zone='i',
        xf18='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xal26_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xal26='float64')

    @remote_function(can_handle_array=True)
    def set_xal26_at_zone(
        index_of_the_particle='i', zone='i',
        xal26='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xneut_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xneut='float64')

    @remote_function(can_handle_array=True)
    def set_xneut_at_zone(
        index_of_the_particle='i', zone='i',
        xneut='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xprot_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xprot='float64')

    @remote_function(can_handle_array=True)
    def set_xprot_at_zone(
        index_of_the_particle='i', zone='i',
        xprot='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_omegi_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (omegi='float64')

    @remote_function(can_handle_array=True)
    def set_omegi_at_zone(
        index_of_the_particle='i', zone='i',
        omegi='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xbid_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xbid='float64')

    @remote_function(can_handle_array=True)
    def set_xbid_at_zone(
        index_of_the_particle='i', zone='i',
        xbid='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xbid1_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (xbid1='float64')

    @remote_function(can_handle_array=True)
    def set_xbid1_at_zone(
        index_of_the_particle='i', zone='i',
        xbid1='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vp_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vp='float64')

    @remote_function(can_handle_array=True)
    def set_vp_at_zone(
        index_of_the_particle='i', zone='i',
        vp='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vt_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vt='float64')

    @remote_function(can_handle_array=True)
    def set_vt_at_zone(
        index_of_the_particle='i', zone='i',
        vt='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vr_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vr='float64')

    @remote_function(can_handle_array=True)
    def set_vr_at_zone(
        index_of_the_particle='i', zone='i',
        vr='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vs_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vs='float64')

    @remote_function(can_handle_array=True)
    def set_vs_at_zone(
        index_of_the_particle='i', zone='i',
        vs='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vx_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vx='float64')

    @remote_function(can_handle_array=True)
    def set_vx_at_zone(
        index_of_the_particle='i', zone='i',
        vx='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vy_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vy='float64')

    @remote_function(can_handle_array=True)
    def set_vy_at_zone(
        index_of_the_particle='i', zone='i',
        vy='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vy3_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vy3='float64')

    @remote_function(can_handle_array=True)
    def set_vy3_at_zone(
        index_of_the_particle='i', zone='i',
        vy3='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxc12_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxc12='float64')

    @remote_function(can_handle_array=True)
    def set_vxc12_at_zone(
        index_of_the_particle='i', zone='i',
        vxc12='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxc13_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxc13='float64')

    @remote_function(can_handle_array=True)
    def set_vxc13_at_zone(
        index_of_the_particle='i', zone='i',
        vxc13='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxn14_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxn14='float64')

    @remote_function(can_handle_array=True)
    def set_vxn14_at_zone(
        index_of_the_particle='i', zone='i',
        vxn14='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxn15_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxn15='float64')

    @remote_function(can_handle_array=True)
    def set_vxn15_at_zone(
        index_of_the_particle='i', zone='i',
        vxn15='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxo16_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxo16='float64')

    @remote_function(can_handle_array=True)
    def set_vxo16_at_zone(
        index_of_the_particle='i', zone='i',
        vxo16='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxo17_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxo17='float64')

    @remote_function(can_handle_array=True)
    def set_vxo17_at_zone(
        index_of_the_particle='i', zone='i',
        vxo17='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxo18_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxo18='float64')

    @remote_function(can_handle_array=True)
    def set_vxo18_at_zone(
        index_of_the_particle='i', zone='i',
        vxo18='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxne20_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxne20='float64')

    @remote_function(can_handle_array=True)
    def set_vxne20_at_zone(
        index_of_the_particle='i', zone='i',
        vxne20='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxne22_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxne22='float64')

    @remote_function(can_handle_array=True)
    def set_vxne22_at_zone(
        index_of_the_particle='i', zone='i',
        vxne22='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxmg24_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxmg24='float64')

    @remote_function(can_handle_array=True)
    def set_vxmg24_at_zone(
        index_of_the_particle='i', zone='i',
        vxmg24='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxmg25_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxmg25='float64')

    @remote_function(can_handle_array=True)
    def set_vxmg25_at_zone(
        index_of_the_particle='i', zone='i',
        vxmg25='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxmg26_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxmg26='float64')

    @remote_function(can_handle_array=True)
    def set_vxmg26_at_zone(
        index_of_the_particle='i', zone='i',
        vxmg26='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxf19_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxf19='float64')

    @remote_function(can_handle_array=True)
    def set_vxf19_at_zone(
        index_of_the_particle='i', zone='i',
        vxf19='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxne21_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxne21='float64')

    @remote_function(can_handle_array=True)
    def set_vxne21_at_zone(
        index_of_the_particle='i', zone='i',
        vxne21='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxna23_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxna23='float64')

    @remote_function(can_handle_array=True)
    def set_vxna23_at_zone(
        index_of_the_particle='i', zone='i',
        vxna23='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxal27_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxal27='float64')

    @remote_function(can_handle_array=True)
    def set_vxal27_at_zone(
        index_of_the_particle='i', zone='i',
        vxal27='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxsi28_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxsi28='float64')

    @remote_function(can_handle_array=True)
    def set_vxsi28_at_zone(
        index_of_the_particle='i', zone='i',
        vxsi28='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxc14_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxc14='float64')

    @remote_function(can_handle_array=True)
    def set_vxc14_at_zone(
        index_of_the_particle='i', zone='i',
        vxc14='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxf18_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxf18='float64')

    @remote_function(can_handle_array=True)
    def set_vxf18_at_zone(
        index_of_the_particle='i', zone='i',
        vxf18='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxal26_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxal26='float64')

    @remote_function(can_handle_array=True)
    def set_vxal26_at_zone(
        index_of_the_particle='i', zone='i',
        vxal26='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxneut_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxneut='float64')

    @remote_function(can_handle_array=True)
    def set_vxneut_at_zone(
        index_of_the_particle='i', zone='i',
        vxneut='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxprot_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxprot='float64')

    @remote_function(can_handle_array=True)
    def set_vxprot_at_zone(
        index_of_the_particle='i', zone='i',
        vxprot='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vomegi_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vomegi='float64')

    @remote_function(can_handle_array=True)
    def set_vomegi_at_zone(
        index_of_the_particle='i', zone='i',
        vomegi='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxbid_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxbid='float64')

    @remote_function(can_handle_array=True)
    def set_vxbid_at_zone(
        index_of_the_particle='i', zone='i',
        vxbid='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vxbid1_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (vxbid1='float64')

    @remote_function(can_handle_array=True)
    def set_vxbid1_at_zone(
        index_of_the_particle='i', zone='i',
        vxbid1='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_abelx_at_zone(
        index_of_the_particle='i', zone='i',
        element='i',
    ):
        returns (abelx='float64')

    @remote_function(can_handle_array=True)
    def set_abelx_at_zone(
        index_of_the_particle='i', zone='i',
        element='i',
        abelx='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vabelx_at_zone(
        index_of_the_particle='i', zone='i',
        element='i',
    ):
        returns (vabelx='float64')

    @remote_function(can_handle_array=True)
    def set_vabelx_at_zone(
        index_of_the_particle='i', zone='i',
        element='i',
        vabelx='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_drl(
        index_of_the_particle='i', number='i'
    ):
        returns (drl='float64')

    @remote_function(can_handle_array=True)
    def set_drl(
        index_of_the_particle='i', number='i',
        drl='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_drte(
        index_of_the_particle='i', number='i'
    ):
        returns (drte='float64')

    @remote_function(can_handle_array=True)
    def set_drte(
        index_of_the_particle='i', number='i',
        drte='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_drp(
        index_of_the_particle='i', number='i'
    ):
        returns (drp='float64')

    @remote_function(can_handle_array=True)
    def set_drp(
        index_of_the_particle='i', number='i',
        drp='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_drt(
        index_of_the_particle='i', number='i'
    ):
        returns (drt='float64')

    @remote_function(can_handle_array=True)
    def set_drt(
        index_of_the_particle='i', number='i',
        drt='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_drr(
        index_of_the_particle='i', number='i'
    ):
        returns (drr='float64')

    @remote_function(can_handle_array=True)
    def set_drr(
        index_of_the_particle='i', number='i',
        drr='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_eps_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps='float64')

    @remote_function(can_handle_array=True)
    def get_epsy_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (epsy='float64')

    @remote_function(can_handle_array=True)
    def get_eps_c_adv_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_c_adv='float64')

    @remote_function(can_handle_array=True)
    def get_eps_ne_adv_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_ne_adv='float64')

    @remote_function(can_handle_array=True)
    def get_eps_o_adv_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_o_adv='float64')

    @remote_function(can_handle_array=True)
    def get_eps_si_adv_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_si_adv='float64')

    @remote_function(can_handle_array=True)
    def get_eps_grav_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_grav='float64')

    @remote_function(can_handle_array=True)
    def get_eps_nu_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (eps_nu='float64')

    @remote_function(can_handle_array=True)
    def get_nabla_rad_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (nabla_rad='float64')

    @remote_function(can_handle_array=True)
    def get_nabla_ad_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (nabla_ad='float64')

    @remote_function(can_handle_array=True)
    def get_nabla_mu_at_zone(
        index_of_the_particle='i', zone='i'
    ):
        returns (nabla_mu='float64')

    @remote_function(can_handle_array=True)
    def get_nbzel(
        index_of_the_particle='i', number='i'
    ):
        returns (nbzel='int32')

    @remote_function(can_handle_array=True)
    def set_nbzel(
        index_of_the_particle='i', number='i',
        nbzel='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nbael(
        index_of_the_particle='i', number='i'
    ):
        returns (nbael='int32')

    @remote_function(can_handle_array=True)
    def set_nbael(
        index_of_the_particle='i', number='i',
        nbael='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_abels(
        index_of_the_particle='i', number='i'
    ):
        returns (abels='float64')

    @remote_function(can_handle_array=True)
    def set_abels(
        index_of_the_particle='i', number='i',
        abels='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xnetalu(
        index_of_the_particle='i', number='i'
    ):
        returns (xnetalu='float64')

    @remote_function(can_handle_array=True)
    def set_xnetalu(
        index_of_the_particle='i', number='i',
        xnetalu='float64'):
        returns ()


class Genec(StellarEvolution, InternalStellarStructure):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(
            self,  GenecInterface(**options), **options)
        # self.stopping_conditions = StoppingConditions(self)
        self.model_time = 0.0 | units.yr

    def define_parameters(self, handler):
        handler.add_method_parameter(
            None,
            "set_genec_path",
            "path_to_data",
            "Path to the data directory",
            default_value=self.data_directory
        )

        # handler.add_method_parameter(
        #     "get_par_stopping_condition",
        #     "set_par_stopping_condition",
        #     "stopping_condition",
        #     "GENEC parameter stopping_condition",
        # )

        # handler.add_method_parameter(
        #     "get_min_timestep_stop_condition",
        #     "set_min_timestep_stop_condition",
        #     "min_timestep_stop_condition",
        #     "The minimum timestep stop condition",
        #     default_value=1e-4 | units.julianyr
        # )

    def define_particle_sets(self, handler):

        for set_name in ['particles', ]:
            handler.define_set(set_name, 'index_of_the_particle')
            InternalStellarStructure.define_particle_sets(
                self, handler, set_name=set_name
            )
            handler.set_new(set_name, 'new_particle')

            for parameter in SCALAR_SETTERS.items():
                if (
                    len(parameter[1]) == 4
                    and parameter[1][3] != ""
                ):
                    names = (
                        parameter[1][3],
                    )
                    handler.add_getter(
                        set_name,
                        f'get_{parameter[0]}',
                        names=names,
                    )
                    handler.add_setter(
                        set_name,
                        f'set_{parameter[0]}',
                        names=names,
                    )
                elif parameter[0] not in [
                    "radius",
                ]:
                    handler.add_method(
                        set_name,
                        f'get_{parameter[0]}',
                    )
                    handler.add_method(
                        set_name,
                        f'set_{parameter[0]}',
                    )

            handler.add_getter(set_name, 'get_radius')
            handler.add_getter(
                set_name, 'get_surface_velocity', names=('surface_velocity',)
            )
            # handler.add_getter(set_name, 'get_mass')
            # handler.add_getter(set_name, 'get_age')
            # handler.add_getter(set_name, 'get_luminosity')
            # handler.add_getter(set_name, 'get_temperature')
            handler.add_getter(set_name, 'get_time_step', names=('time_step',))
            handler.add_setter(set_name, 'set_time_step', names=('time_step',))
            # handler.add_getter(
            #     set_name, 'get_number_of_species', names=('n_species',)
            # )

            # handler.add_method(set_name, 'get_radius_profile')
            # handler.add_method(set_name, 'get_temperature_profile')
            # handler.add_method(set_name, 'get_luminosity_profile')
            handler.add_method(set_name, 'get_phase')
            handler.add_method(set_name, 'get_mass_profile')
            handler.add_method(set_name, 'get_eps_profile')
            handler.add_method(set_name, 'get_epsy_profile')
            handler.add_method(set_name, 'get_eps_c_adv_profile')
            handler.add_method(set_name, 'get_eps_ne_adv_profile')
            handler.add_method(set_name, 'get_eps_o_adv_profile')
            handler.add_method(set_name, 'get_eps_si_adv_profile')
            handler.add_method(set_name, 'get_eps_grav_profile')
            handler.add_method(set_name, 'get_eps_nu_profile')
            handler.add_method(set_name, 'get_nabla_rad_profile')
            handler.add_method(set_name, 'get_nabla_ad_profile')
            handler.add_method(set_name, 'get_nabla_mu_profile')
            handler.add_method(set_name, 'get_luminosity_profile')
            handler.add_method(set_name, 'get_cumulative_mass_profile')
            handler.add_method(set_name, 'get_mass_fraction_at_zone')

            handler.add_method(set_name, 'evolve_one_step')
            handler.add_method(set_name, 'evolve_for')
            handler.set_delete(set_name, 'delete_star')

            handler.add_method(set_name, 'get_star_name')
            handler.add_method(set_name, 'get_internal_structure')

            # for star_prop in [
            #         "nabla_rad", "nabla_ad", "nabla_mu", "eps", "epsy",
            #         "eps_c_adv", "eps_ne_adv", "eps_o_adv", "eps_si_adv",
            #         "eps_grav", "eps_nu", "radius", "temperature", "density",
            #         "luminosity", "pressure",
            # ]:
            #     handler.add_gridded_getter(
            #         set_name,
            #         f'get_{star_prop}_at_zone', 'get_firstlast_zone',
            #         names=(f'{star_prop}_profile',)
            #     )
            #     # handler.add_gridded_setter(
            #     #     set_name,
            #     #     'set_temperature_at_zone', 'get_firstlast_zone',
            #     #     names=('temperature_profile',)
            #     # )

            # for species in SPECIES_NAMES:
            #     handler.add_gridded_getter(
            #         set_name,
            #         f'get_mass_fraction_of_{species}_at_zone',
            #         'get_firstlast_zone',
            #         names=(f'abundance_{species}',)
            #     )

            #     handler.add_gridded_setter(
            #         set_name,
            #         f'set_mass_fraction_of_{species}_at_zone',
            #         'get_firstlast_zone',
            #         names=(f'abundance_{species}',)
            #     )

    # @property
    # def particles(self):
    #     basic_attributes = [
    #         "age", "mass", "radius", "temperature", "luminosity",
    #     ]
    #     return ParticlesWithFilteredAttributes(
    #         self.fullparticles,
    #         basic_attributes,
    #     )

    def define_state(self, handler):
        common.CommonCode.define_state(self, handler)
        # StellarEvolution.define_state(self, handler)
        # InternalStellarStructure.define_state(self, handler)

        # Only allow setting of star_name in EDIT or UPDATE states
        # I.e. must do initialize_code and commit_parameters FIRST!

        # Initialized (initialize_code)
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')

        # -> Edit (commit_parameters)
        handler.add_transition('INITIALIZED', 'EDIT', 'commit_parameters')
        # handler.add_method('EDIT', 'set_star_name')
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('UPDATE', 'new_particle')

        handler.add_transition(
            'RUN', 'UPDATE', 'finalize_stellar_model', False)

        # -> Run (commit_particles)
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('EDIT', 'UPDATE', 'commit_particles')
        handler.add_method('RUN', 'evolve_one_step')
        handler.add_method('!UPDATE', 'evolve_one_step')

        # -> Update
        handler.add_transition('RUN', 'UPDATE', 'finalize_stellar_model')

        for state in ["EDIT", "UPDATE"]:
            for method in (
                "new_particle_from_model", "new_stellar_model",
                "set_abelx_at_zone", "set_vabelx_at_zone", "set_nbzel",
                "set_nbael", "set_abels", "set_xnetalu"
            ):
                handler.add_method(state, method)

        for state in ["UPDATE", "RUN"]:
            for parameter in ALL_GETTERS:
                handler.add_method(state, f'get_{parameter[0]}')
            handler.add_method(state, 'get_internal_structure')
            handler.add_method(state, 'get_radius')
            handler.add_method(state, 'get_number_of_species')
            handler.add_method(state, 'get_temperature')
            handler.add_method(state, 'get_luminosity')
            handler.add_method(state, 'get_time_step')
            handler.add_method(state, 'get_mass')
            handler.add_method(state, 'get_age')
            handler.add_method(state, 'get_surface_velocity')

            handler.add_method(state, 'get_chemical_abundance_profiles')
            handler.add_method(state, 'get_mass_fraction_at_zone')
            handler.add_method(state, 'get_mass_fraction_of_species_at_zone')
            handler.add_method(state, 'get_mu_at_zone')
            handler.add_method(state, 'get_pressure_at_zone')
            handler.add_method(state, 'get_radius_at_zone')
            # handler.add_method(state, 'get_eps_at_zone')
            # handler.add_method(state, 'get_epsy_at_zone')
            # handler.add_method(state, 'get_eps_c_adv_at_zone')
            # handler.add_method(state, 'get_eps_ne_adv_at_zone')
            # handler.add_method(state, 'get_eps_o_adv_at_zone')
            # handler.add_method(state, 'get_eps_si_adv_at_zone')
            # handler.add_method(state, 'get_eps_grav_at_zone')
            # handler.add_method(state, 'get_eps_nu_at_zone')
            handler.add_method(state, 'get_temperature_at_zone')
            handler.add_method(state, 'get_density_at_zone')
            # handler.add_method(state, 'get_luminosity_at_zone')
            for species in SPECIES_NAMES:
                handler.add_method(
                    state, f'get_mass_fraction_of_{species}_at_zone'
                )
        for state in ["EDIT", "UPDATE", "!RUN"]:
            for parameter in ALL_SETTERS:
                handler.add_method(state, f'set_{parameter[0]}')
            handler.add_method(state, 'set_time_step')
        for state in ["!UPDATE"]:
            handler.add_method(state, "set_bintide")
            handler.add_method(state, "set_zams_velocity")
        for state in ["!EDIT"]:
            handler.add_method(state, "set_zams_velocity")

        handler.add_method('UPDATE', 'set_n_snap')
        # handler.add_method('UPDATE', 'set_ipoly')
        handler.add_method('UPDATE', 'set_chemical_abundance_profiles')
        handler.add_method('UPDATE', 'set_mass_fraction_of_species_at_zone')
        # handler.add_method('UPDATE', 'set_radius')
        # handler.add_method('UPDATE', 'set_pressure_at_zone')
        handler.add_method('UPDATE', 'set_radius_at_zone')
        handler.add_method('UPDATE', 'set_temperature_at_zone')
        handler.add_method('UPDATE', 'set_density_at_zone')
        # handler.add_method('UPDATE', 'set_luminosity_at_zone')
        handler.add_method('UPDATE', 'set_mass_fraction_at_zone')
        # handler.add_method('UPDATE', 'set_time_step')
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        for state in ["UPDATE"]:
            for species in SPECIES_NAMES:
                handler.add_method(
                    state, f'set_mass_fraction_of_{species}_at_zone'
                )
        # -> Run (recommit_particles)

        # handler.add_method('UPDATE', 'set_star_name')
        # handler.add_method('UPDATE', 'new_particle')

    def define_methods(self, handler):
        InternalStellarStructure.define_methods(self, handler)
        StellarEvolution.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (units.MSun, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
        )

        # handler.add_method(
        #     "read_genec_model",
        #     (handler.NO_UNIT),
        #     (handler.INDEX, handler.ERROR_CODE)
        # )

        # specifically add this here since the unit is different
        handler.add_method(  
            "get_time_step",
            (handler.INDEX,),
            (units.s, handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_surface_velocity",
            (handler.INDEX),
            (units.km/units.s, handler.ERROR_CODE,)
        )

    def get_eps_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_epsy_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying epsy profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_epsy_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_c_adv_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_c_adv_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_ne_adv_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_ne_adv_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_o_adv_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_o_adv_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_si_adv_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_si_adv_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_grav_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_grav_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_eps_nu_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying eps_c_adv profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_eps_nu_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_mass_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying mass profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_mass_fraction_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_cumulative_mass_profile(
            self, indices_of_the_stars, number_of_zones=None):
        frac_profile = self.get_mass_profile(
            indices_of_the_stars, number_of_zones=number_of_zones)
        return frac_profile.cumsum()

    def get_nabla_rad_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying nabla_rad profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_nabla_rad_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_nabla_ad_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying nabla_ad profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_nabla_ad_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_nabla_mu_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying nabla_mu profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_nabla_mu_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_luminosity_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying luminosity profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_m(indices_of_the_stars)
        return self.get_luminosity_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_internal_structure(self, index_of_the_star):
        internal_structure = {
            # 'veryFirst': self.get_veryFirst(index_of_the_star),
        }
        for star_property in SCALAR_SETTERS:
            func = getattr(self, f'get_{star_property}')
            internal_structure[star_property] = func(index_of_the_star)
        number_of_zones = self.get_m(index_of_the_star)
        for structure_property in {
            **GENEC_ZONE,
        }:
            func = getattr(self, f'get_{structure_property}_at_zone')
            # internal_structure.append(func(index_of_the_star))

            internal_structure[structure_property] = func(
                index_of_the_star,
                (list(range(number_of_zones)) | units.none)[::-1]
            )
        extra_elements = self.get_mbelx(index_of_the_star)
        for structure_property in {
            **GENEC_ARRAY_MBELX_M,
        }:
            func = getattr(self, f'get_{structure_property}_at_zone')
            x = np.array([])
            for element in range(extra_elements):
                x = np.append(
                    x, func(
                        index_of_the_star,
                        (list(range(number_of_zones)) | units.none)[::-1],
                        1+element
                    )
                )
            internal_structure[structure_property] = x
        for netdef_property in {
            **GENEC_NETDEF_ARRAYS,
        }:
            func = getattr(self, f'get_{netdef_property}')

            internal_structure[netdef_property] = func(
                index_of_the_star,
                (list(range(1, 9)) | units.none)
            )
        for netalu_property in {
            **GENEC_NETALU_ARRAYS,
        }:
            func = getattr(self, f'get_{netalu_property}')

            internal_structure[netalu_property] = func(
                index_of_the_star,
                (list(range(1, 6)) | units.none)
            )
        for array_3_property in {
            **GENEC_ARRAY_3,
        }:
            func = getattr(self, f'get_{array_3_property}')

            internal_structure[array_3_property] = func(
                index_of_the_star,
                (list(range(1, 4)) | units.none)
            )

        return internal_structure

    def new_particle_from_model(
        self, internal_structure, current_age=0 | units.yr, key=None
    ):
        index_of_the_particle = self.new_stellar_model(
            # Extra - for AMUSE
            internal_structure['modell'],
            internal_structure['veryFirst'],
            # Characteristics
            internal_structure['initialised'],
            internal_structure['star_name'],
            internal_structure['nwmd'],
            internal_structure['nwseq'],
            internal_structure['modanf'],
            internal_structure['nzmod'],
            internal_structure['end_at_phase'],
            internal_structure['end_at_model'],
            # Physics
            internal_structure['irot'],
            internal_structure['isol'],
            internal_structure['imagn'],
            internal_structure['ialflu'],
            internal_structure['ianiso'],
            internal_structure['ipop3'],
            internal_structure['ibasnet'],
            internal_structure['phase'],
            internal_structure['var_rates'],
            internal_structure['bintide'],
            internal_structure['binm2'],
            internal_structure['periodini'],
            internal_structure['const_per'],
            internal_structure['iprezams'],
            # Composition
            internal_structure['initial_metallicity'],
            internal_structure['zsol'],
            internal_structure['z'],
            internal_structure['iopac'],
            internal_structure['ikappa'],
            # Rotation
            internal_structure['idiff'],
            internal_structure['iadvec'],
            internal_structure['istati'],
            internal_structure['icoeff'],
            internal_structure['fenerg'],
            internal_structure['richac'],
            internal_structure['igamma'],
            internal_structure['frein'],
            internal_structure['K_Kawaler'],
            internal_structure['Omega_saturation'],
            internal_structure['rapcrilim'],
            internal_structure['zams_velocity'],
            internal_structure['xfom'],
            internal_structure['omega'],
            internal_structure['xdial'],
            internal_structure['idialo'],
            internal_structure['idialu'],
            internal_structure['Add_Flux'],
            internal_structure['diff_only'],
            internal_structure['B_initial'],
            internal_structure['add_diff'],
            internal_structure['n_mag'],
            internal_structure['alpha_F'],
            internal_structure['nsmooth'],
            internal_structure['qminsmooth'],
            # Surface
            internal_structure['imloss'],
            internal_structure['fmlos'],
            internal_structure['ifitm'],
            internal_structure['fitm'],
            internal_structure['fitmi'],
            internal_structure['fitmi_default'],
            internal_structure['deltal'],
            internal_structure['deltat'],
            internal_structure['nndr'],
            internal_structure['RSG_Mdot'],
            internal_structure['SupraEddMdot'],
            internal_structure['Be_mdotfrac'],
            internal_structure['start_mdot'],
            # Convection
            internal_structure['iledou'],
            internal_structure['idifcon'],
            internal_structure['iover'],
            internal_structure['elph'],
            internal_structure['my'],
            internal_structure['dovhp'],
            internal_structure['iunder'],
            internal_structure['dunder'],
            # Convergence
            internal_structure['gkorm'],
            internal_structure['alph'],
            internal_structure['agdr'],
            internal_structure['faktor'],
            internal_structure['dgrp'],
            internal_structure['dgrl'],
            internal_structure['dgry'],
            internal_structure['dgrc'],
            internal_structure['dgro'],
            internal_structure['dgr20'],
            internal_structure['nbchx'],
            internal_structure['nrband'],
            # Time
            internal_structure['xcn'],
            internal_structure['islow'],
            internal_structure['icncst'],
            internal_structure['tauH_fit'],
            # Various
            internal_structure['display_plot'],
            internal_structure['iauto'],
            internal_structure['iprn'],
            internal_structure['iout'],
            internal_structure['itmin'],
            internal_structure['xyfiles'],
            internal_structure['idebug'],
            internal_structure['itests'],
            internal_structure['verbose'],
            internal_structure['stop_deg'],
            internal_structure['n_snap'],
            # Properties
            internal_structure['gms'],
            internal_structure['alter'],
            internal_structure['gls'],
            internal_structure['teff'],
            internal_structure['glsv'],
            internal_structure['teffv'],
            internal_structure['dzeitj'],
            internal_structure['dzeit'],
            internal_structure['dzeitv'],
            internal_structure['xmini'],
            internal_structure['ab'],
            internal_structure['dm_lost'],
            internal_structure['m'],
            internal_structure['summas'],
            internal_structure['dk'],
            internal_structure['rlp'],
            internal_structure['rlt'],
            internal_structure['rlc'],
            internal_structure['rrp'],
            internal_structure['rrt'],
            internal_structure['rrc'],
            internal_structure['rtp'],
            internal_structure['rtt'],
            internal_structure['rtc'],
            internal_structure['tdiff'],
            internal_structure['suminenv'],
            internal_structure['xltotbeg'],
            internal_structure['dlelexprev'],
            internal_structure['radius'],
            internal_structure['zams_radius'],
            internal_structure['mbelx'],
            internal_structure['xtefflast'],
            internal_structure['xllast'],
            internal_structure['xrholast'],
            internal_structure['xclast'],
            internal_structure['xtclast'],
            internal_structure['inum'],
            internal_structure['nsugi'],
            internal_structure['period'],
            internal_structure['r_core'],
            internal_structure['vna'],
            internal_structure['vnr'],
            # Zone
            internal_structure['q'],
            internal_structure['p'],
            internal_structure['t'],
            internal_structure['r'],
            internal_structure['s'],
            internal_structure['x'],
            internal_structure['y3'],
            internal_structure['y'],
            internal_structure['xc12'],
            internal_structure['xc13'],
            internal_structure['xn14'],
            internal_structure['xn15'],
            internal_structure['xo16'],
            internal_structure['xo17'],
            internal_structure['xo18'],
            internal_structure['xne20'],
            internal_structure['xne22'],
            internal_structure['xmg24'],
            internal_structure['xmg25'],
            internal_structure['xmg26'],
            internal_structure['xf19'],
            internal_structure['xne21'],
            internal_structure['xna23'],
            internal_structure['xal27'],
            internal_structure['xsi28'],
            internal_structure['xc14'],
            internal_structure['xf18'],
            internal_structure['xal26'],
            internal_structure['xneut'],
            internal_structure['xprot'],
            internal_structure['omegi'],
            internal_structure['xbid'],
            internal_structure['xbid1'],
            internal_structure['vp'],
            internal_structure['vt'],
            internal_structure['vr'],
            internal_structure['vs'],
            internal_structure['vx'],
            internal_structure['vy'],
            internal_structure['vy3'],
            internal_structure['vxc12'],
            internal_structure['vxc13'],
            internal_structure['vxn14'],
            internal_structure['vxn15'],
            internal_structure['vxo16'],
            internal_structure['vxo17'],
            internal_structure['vxo18'],
            internal_structure['vxne20'],
            internal_structure['vxne22'],
            internal_structure['vxmg24'],
            internal_structure['vxmg25'],
            internal_structure['vxmg26'],
            internal_structure['vxf19'],
            internal_structure['vxne21'],
            internal_structure['vxna23'],
            internal_structure['vxal27'],
            internal_structure['vxsi28'],
            internal_structure['vxc14'],
            internal_structure['vxf18'],
            internal_structure['vxal26'],
            internal_structure['vxneut'],
            internal_structure['vxprot'],
            internal_structure['vomegi'],
            internal_structure['vxbid'],
            internal_structure['vxbid1'],
            # netdef scalars
            internal_structure['xlostneu'],
            # internal_structure['abelx'],
            # internal_structure['vabelx'],
        )
        extra_elements = internal_structure['mbelx']
        zones = internal_structure['m']
        for element in range(extra_elements):
            for zone in range(zones):
                self.set_abelx_at_zone(
                    index_of_the_particle, zone, 1+element,
                    internal_structure['abelx'][element*zones + zone]
                )
                self.set_vabelx_at_zone(
                    index_of_the_particle, zone, 1+element,
                    internal_structure['vabelx'][element*zones + zone]
                )
        for i in range(8):
            self.set_nbzel(
                index_of_the_particle, i+1,
                internal_structure['nbzel'][i]
            )
            self.set_nbael(
                index_of_the_particle, i+1,
                internal_structure['nbael'][i]
            )
            self.set_abels(
                index_of_the_particle, i+1,
                internal_structure['abels'][i]
            )
        for i in range(5):
            self.set_xnetalu(
                index_of_the_particle, i+1,
                internal_structure['xnetalu'][i]
            )
        for i in range(3):
            self.set_drl(
                index_of_the_particle, i+1,
                internal_structure['drl'][i]
            )
            self.set_drte(
                index_of_the_particle, i+1,
                internal_structure['drte'][i]
            )
            self.set_drp(
                index_of_the_particle, i+1,
                internal_structure['drp'][i]
            )
            self.set_drt(
                index_of_the_particle, i+1,
                internal_structure['drt'][i]
            )
            self.set_drr(
                index_of_the_particle, i+1,
                internal_structure['drr'][i]
            )
        return
