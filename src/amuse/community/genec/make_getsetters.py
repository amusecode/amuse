import sys
from amuse.community import NO_UNIT

PY_INTERFACE_HEADER = """
\"\"\"
Interface for GENEC.
\"\"\"

import numpy as np
from amuse.datamodel import ParticlesWithFilteredAttributes
from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function, remote_function
from amuse.community import LiteratureReferencesMixIn
from amuse.community import CodeWithDataDirectories
from amuse.community import InCodeComponentImplementation
from amuse.community import NO_UNIT, ERROR_CODE
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
"""

# Parameters (but individual to each star)
GENEC_STAR_CHARACTERISTICS = {
    # 'GENEC name': [dtype, unit, description, AMUSE name (optional)]
    'initialised': ['bool', '', "True if the star is an intialised model"],
    'starname': ['string', '', "Name of the star"],
    'nwmd': ['int32', '', "model number"],
    'nwseq': ['int32', '', "number of the first model in the time-step series"],
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
    'ianiso': ['int32', '', "wind anisotropy if set to 1"],
    'ipop3': ['int32', '', "Z=0 models if set to 1"],
    'ibasnet': ['int32', '', "extended nuclear network if set to 1"],
    'phase': ['int32', '', "fusion phases"],
    'var_rates': ['bool', '', "allows to use different reaction rate files if set to True"],
    'bintide': ['bool', '', "tidal interaction in binaries if set to True"],
    'binm2': ['float64', '', "mass of the companion"],
    'periodini': ['float64', '', "initial period of the binary"],
    'const_per': ['bool', '', "keep constant period if True"],
    'iprezams': ['int32', '', ""],
}

GENEC_STAR_COMPOSITION = {
    'zinit': ['float64', '', "initial metallicity of the model"],
    'zsol': ['float64', '', "reference solar metallicity"],
    'z': ['float64', '', "abundance of the neglected isotopes"],
    'iopac': ['int32', '', "choice of the opacity table if ikappa = 5"],
    'ikappa': ['int32', '', "opacity choice"],
}

GENEC_STAR_ROTATION = {
    'idiff': ['int32', '', "computation of the diffusion of Omega and chemicals if set to 1"],
    'iadvec': ['int32', '', "advecto-diffusive version of the transport of Omega if set to 1"],
    'istati': ['int32', '', "only local conservation of angular momentum if set to 1"],
    'icoeff': ['int32', '', "prescription to be used for the diffusion coefficients"],
    'fenerg': ['float64', '', "fraction of the shear energy used for the mixing"],
    'richac': ['float64', '', "critical Richardson number value"],
    'igamma': ['int32', '', "treatment of the shear according to M&M96 if set to 1"],
    'frein': ['float64', '', "magnetic braking if ≠ 0"],
    'K_Kawaler': ['float64', '', "K parameter in Kawaler1998 magnetic braking law"],
    'Omega_saturation': ['float64', '', "Ωsat in Kawaler magnetic braking law"],
    'rapcrilim': ['float64', '', "maximum Ωcrit ratio before the onset of mechanical mass loss"],
    'vwant': ['float64', '', "chosen velocity on the ZAMS (Veq if > 1.0, V/Vcrit if 0<vwant<1, Ω/Ωcrit if -1<vwant<0"],
    'xfom': ['float64', '', "multiplying factor for surface Ω"],
    'omega': ['float64', '', "surface Ω"],
    'xdial': ['float64', '', ""],
    'idialo': ['int32', '', ""],
    'idialu': ['int32', '', ""],
    'Add_Flux': ['bool', '', "improves the angular momentum conservation treatment if set to True"],
    'diff_only': ['bool', '', "applies the tidal mixing only in timesteps where we use diffusion (and not advection) if set to True"],
    'B_initial': ['float64', '', "if >0, switches on the wind quenching, and evolves the magnetic field strength according to the conservation of magnetic flux"],
    'add_diff': ['float64', '', "additional viscosity value, used to increase the core-envelope coupling"],
    'n_mag': ['int32', '', "type of treatment for magnetic fields"],
    'alpha_F': ['float64', '', "α parameter in Fuller+ 2019. Values higher than 1 lower the threshold to trigger the instability (Qmin) and increase the magnetic viscosity (increased transport)"],
    'nsmooth': ['int32', '', "number of layers used for smoothing the Ω gradient (used by Mag_diff_general). Default value is nsmooth=1, recommended value for Fuller+ 2019 is nsmooth=5"],
    'qminsmooth': ['bool', '', "mag. instability always considered if set to True"],
}

GENEC_STAR_SURFACE = {
    'imloss': ['int32', '', "choice of the mass-loss prescription"],
    'fmlos': ['float64', '', "except for IMLOSS=2 or 3, multiplying factor applied to the mass loss"],
    'ifitm': ['int32', '', "management of the changes of fitm during redward evolution (and/or blueward evolution after a RSG phase)"],
    'fitm': ['float64', '', "mass included in the interior"],
    'fitmi': ['float64', '', "max value of FITM to which the star will come back when going back to the blue"],
    'deltal': ['float64', '', "triangle size for L at the surface"],
    'deltat': ['float64', '', "triangle size for T at the surface"],
    'nndr': ['int32', '', "management of the behaviour of (L,Teff) with respect to the triangle"],
    'RSG_Mdot': ['int32', '', "mass-loss recipe used for RSG"],
    'SupraEddMdot': ['bool', '', "x3 multiplication factor to the wind in case of supra-Eddington layers if set to True"],
    'Be_mdotfrac': ['float64', '', ""],
    'start_mdot': ['float64', '', "value of Ω/Ωcrit at which the Be mass loss starts to apply"],
}

GENEC_STAR_CONVECTION = {
    'iledou': ['int32', '', "Ledoux criterion for convection if set to 1"],
    'idifcon': ['int32', '', "convection treated as a diffusion if set to 1 (used during O-b and Si-b)"],
    'iover': ['int32', '', "overshooting taken into account if set to 1"],
    'elph': ['float64', '', "mixing length for the external convective zone"],
    'my': ['int32', '', "integration of the envelope on ρ rather than P if set to 1"],
    'dovhp': ['float64', '', "value of the α overshooting parameter Λ = α HP"],
    'iunder': ['int32', '', "undershooting taken into account if set to 1"],
    'dunder': ['float64', '', "value of the undershooting parameter"],
}

GENEC_STAR_CONVERGENCE = {
    'gkorm': ['float64', '', "accepted deviation on the structure variables"],
    'alph': ['float64', '', "fraction of the correction applied for the next iteration"],
    'agdr': ['float64', '', "absolute value of the maximum correction accepted on r, s, P and T"],
    'faktor': ['float64', '', "parameter allowing to 'compensate' (during the computation only) for layers with Lr > Ltot, and hence avoid dL/dr < 0"],
    'dgrp': ['float64', '', "relative variations (from one layer to the other) accepted on P"],
    'dgrl': ['float64', '', "relative variations (from one layer to the other) accepted on L"],
    'dgry': ['float64', '', "relative variations (from one layer to the other) accepted on 4He"],
    'dgrc': ['float64', '', "relative variations (from one layer to the other) accepted on 12C"],
    'dgro': ['float64', '', "relative variations (from one layer to the other) accepted on 16O"],
    'dgr20': ['float64', '', "relative variations (from one layer to the other) accepted on ?"],
    'nbchx': ['int32', '', "iteration number for the computation of the change in chemical composition"],
    'nrband': ['int32', '', "iteration number for the chemicals between model n and n+1"],
}

GENEC_STAR_TIME = {
    'xcn': ['float64', '', "multiplying factor applied on the time step for the next run"],
    'islow': ['int32', '', "slow version of the program if not 0 by modification of the ideal nuclear time step ratxcn"],
    'icncst': ['int32', '', "constant time step (equivalent to xcn=1.0)"],
    'tauH_fit': ['int32', '', "used to set the maximal timestep in case of critical velocity, as a fraction of the MS lifetime"],
}

GENEC_STAR_VARIOUS = {
    'display_plot': ['bool', '', "display of the pgplot window if set to True"],
    'iauto': ['int32', '', "management of the parameters change through the different phases of the evolution"],
    'iprn': ['int32', '', "the full structure is printed every iprn models. If iprn > nzmod, only one structure will be printed, the one corresponding to model nwseq"],
    'iout': ['int32', '', "number of layers used for the smoothing of the diffusion gradient at the border of the convective zones"],
    'itmin': ['int32', '', "minimal number of iterations in henyey"],
    'xyfiles': ['bool', '', ""],
    'idebug': ['int32', '', "addition of some terminal printings useful for the debugging, with different possible levels. All values > 0 set verbose=True"],
    'itests': ['int32', '', "use of a flag to test some pieces of code"],
    'verbose': ['bool', '', "increases the level of printings on the terminal and the .l file"],
    'stop_deg': ['bool', '', "automatically stops a computation when Tc becomes degenerate"],
    'n_snap': ['int32', '', "number of steps between snapshots [0]"],
}

# Stellar properties (but global for the star)
GENEC_STAR_PROPERTIES = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (empty = not used)]
    'm': ['int32', '', "number of zones", "n_zones"],
    'gms': ['float64', 'MSun', "total mass", "mass"],
    'alter': ['float64', 'julianyr', "stellar age", "age"],
    'gls': ['float64', 'LSun', "stellar luminosity", "luminosity"],
    'teff': ['float64', 'K', "effective temperature", "temperature"],
    'glsv': ['float64', 'LSun', "previous luminosity"],
    'teffv': ['float64', 'K', "previous temperature"],
    'dzeitj': ['float64', 'julianyr', "time step (yr)", "time_step"],
    'dzeit': ['float64', 's', "time step (s)", ""],
    'dzeitv': ['float64', 's', "previous time step", ""],
    'xmini': ['float64', 'MSun', "initial mass", "initial_mass"],
    'summas': ['float64', 'MSun', "total mass", ""],
    'ab': ['float64', '', "binary separation", ""],
    'dm_lost': ['float64', 'MSun', "total mass lost", ""],
}

GENEC_STAR_EXTRA = {
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

# Structural properties (m layers)
GENEC_STAR_STRUCTURE = {
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
GENEC_STAR_STRUCTURE_EXTRA = {
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
GENEC_STAR_STRUCTURE_DERIVED = {
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
    **GENEC_STAR_EXTRA,
    **GENEC_STAR_STRUCTURE,
    **GENEC_STAR_STRUCTURE_EXTRA,
    **GENEC_NETDEF_SCALARS,
}
SCALAR_SETTERS = {
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
    **GENEC_STAR_EXTRA,
    **GENEC_NETDEF_SCALARS,
}

VECTOR_SETTERS = {
    **GENEC_STAR_STRUCTURE,
    **GENEC_STAR_STRUCTURE_EXTRA,
}


def star_parameters(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
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
        **GENEC_STAR_EXTRA,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}(index_of_the_particle='i'):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}(\n"
                f"        index_of_the_particle='i',\n"
                f"        {name}={dtype},\n"
                f"    ):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}(index_of_the_particle, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(out):: {name}\n"
                f"    {name} = GenecStar%{name}\n"
                f"    get_{name} = 0\n"
                f"end function get_{name}\n"
                f"\n"
                f"integer function set_{name}(index_of_the_particle, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    GenecStar%{name} = {name}\n"
                f"    set_{name} = 0\n"
                f"end function set_{name}\n"
                f"\n"
            )
    return

def star_structure(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_STAR_STRUCTURE,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}_at_zone(\n"
                f"        index_of_the_particle='i', zone='i'\n"
                f"    ):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}_at_zone(\n"
                f"        index_of_the_particle='i', zone='i',\n"
                f"        {name}={dtype}):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}_at_zone(index_of_the_particle, zone, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: zone\n"
                f"    real(kindreal), intent(out):: {name}\n"
                f"    integer:: i\n"
                f"    get_{name}_at_zone = -1\n"
                f"    i = GenecStar%m-zone\n"
                # f"    if (zone <= GenecStar%m) then\n"
                f"    if (.true.) then\n"
                f"        {name} = GenecStar%{name}(i)\n"
                f"        get_{name}_at_zone = 0\n"
                f"    end if\n"
                f"end function get_{name}_at_zone\n"
                f"\n"
                f"integer function set_{name}_at_zone(index_of_the_particle, zone, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: zone\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    integer:: i\n"
                f"    set_{name}_at_zone = -1\n"
                f"    i = GenecStar%m-zone\n"
                # f"    if (zone <= GenecStar%m) then\n"
                f"    if (.true.) then\n"
                f"        GenecStar%{name}(i) = {name}\n"
                f"        set_{name}_at_zone = 0\n"
                f"    end if\n"
                f"end function set_{name}_at_zone\n"
                f"\n"
            )

def star_structure_extra(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_STAR_STRUCTURE_EXTRA,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}_at_zone(\n"
                f"        index_of_the_particle='i', zone='i',\n"
                f"        element='i',\n"
                f"    ):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}_at_zone(\n"
                f"        index_of_the_particle='i', zone='i',\n"
                f"        element='i',\n"
                f"        {name}={dtype}):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}_at_zone(index_of_the_particle, zone, element, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: zone\n"
                f"    integer, intent(in):: element\n"
                f"    real(kindreal), intent(out):: {name}\n"
                f"    integer:: i\n"
                f"    get_{name}_at_zone = -1\n"
                f"    i = GenecStar%m-zone\n"
                f"    if (zone <= GenecStar%m) then\n"
                f"        {name} = GenecStar%{name}(element,i)\n"
                f"        get_{name}_at_zone = 0\n"
                f"    end if\n"
                f"end function get_{name}_at_zone\n"
                f"\n"
                f"integer function set_{name}_at_zone(index_of_the_particle, zone, element, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: zone\n"
                f"    integer, intent(in):: element\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    integer:: i\n"
                f"    set_{name}_at_zone = -1\n"
                f"    i = GenecStar%m-zone\n"
                # f"    if (zone <= GenecStar%m) then\n"
                f"    if (.true.) then\n"
                f"        GenecStar%{name}(element,i) = {name}\n"
                f"        set_{name}_at_zone = 0\n"
                f"    end if\n"
                f"end function set_{name}_at_zone\n"
                f"\n"
            )

def star_structure_derived(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_STAR_STRUCTURE_DERIVED,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}_at_zone(\n"
                f"        index_of_the_particle='i', zone='i'\n"
                f"    ):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}_at_zone(index_of_the_particle, zone, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: zone\n"
                f"    real(kindreal), intent(out):: {name}\n"
                f"    integer:: i\n"
                f"    get_{name}_at_zone = -1\n"
                f"    i = GenecStar%m-zone\n"
                # f"    if (zone <= GenecStar%m) then\n"
                f"    if (.true.) then\n"
                f"        {name} = GenecStar%{name}(i)\n"
                f"        get_{name}_at_zone = 0\n"
                f"    end if\n"
                f"end function get_{name}_at_zone\n"
                f"\n"
            )


def netdef_scalars(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_NETDEF_SCALARS,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}(index_of_the_particle='i'):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}(\n"
                f"        index_of_the_particle='i',\n"
                f"        {name}={dtype},\n"
                f"    ):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}(index_of_the_particle, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(out):: {name}\n"
                f"    {name} = GenecStar%{name}\n"
                f"    get_{name} = 0\n"
                f"end function get_{name}\n"
                f"\n"
                f"integer function set_{name}(index_of_the_particle, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    GenecStar%{name} = {name}\n"
                f"    set_{name} = 0\n"
                f"end function set_{name}\n"
                f"\n"
            )
    return


def netdef_arrays(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_NETDEF_ARRAYS,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}(\n"
                f"        index_of_the_particle='i', number='i'\n"
                f"    ):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}(\n"
                f"        index_of_the_particle='i', number='i',\n"
                f"        {name}={dtype}):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}(index_of_the_particle, i, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: i\n"
                f"    {f_dtype}, intent(out):: {name}\n"
                f"    get_{name} = -1\n"
                f"    if ((i > 0) .and. (i < 9)) then\n"
                f"        {name} = GenecStar%{name}(i)\n"
                f"        get_{name} = 0\n"
                f"    end if\n"
                f"end function get_{name}\n"
                f"\n"
                f"integer function set_{name}(index_of_the_particle, i, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: i\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    set_{name} = -1\n"
                f"    if ((i > 0) .and. (i < 9)) then\n"
                f"        GenecStar%{name}(i) = {name}\n"
                f"        set_{name} = 0\n"
                f"    end if\n"
                f"end function set_{name}\n"
                f"\n"
            )

def netalu_arrays(outfilename):
    py_outfile = f"{outfilename}.py"
    f_outfile = f"{outfilename}.f90"
    for parameter in {
        **GENEC_NETALU_ARRAYS,
    }.items():
        name = parameter[0]
        dtype = parameter[1][0]
        unit = parameter[1][1]
        description = parameter[1][2]
        if len(parameter[1]) > 3:
            amuse_name = parameter[1][3]
        else:
            amuse_name = None
        f_dtype = get_f_dtype(dtype)
        if unit != '':
            dtype = f"'{dtype}' | units.{unit}"
        else:
            dtype = f"'{dtype}'"
        with open(py_outfile, "a") as py_out:
            py_out.write(
                f"    @remote_function(can_handle_array=True)\n"
                f"    def get_{name}(\n"
                f"        index_of_the_particle='i', number='i'\n"
                f"    ):\n"
                f"        returns ({name}={dtype})\n"
                f"\n"
                f"    @remote_function(can_handle_array=True)\n"
                f"    def set_{name}(\n"
                f"        index_of_the_particle='i', number='i',\n"
                f"        {name}={dtype}):\n"
                f"        returns ()\n"
                f"\n"
            )
        with open(f_outfile, "a") as f_out:
            f_out.write(
                f"integer function get_{name}(index_of_the_particle, i, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: i\n"
                f"    real(kindreal), intent(out):: {name}\n"
                f"    get_{name} = -1\n"
                f"    if ((i > 0) .and. (i < 6)) then\n"
                f"        {name} = GenecStar%{name}(i)\n"
                f"        get_{name} = 0\n"
                f"    end if\n"
                f"end function get_{name}\n"
                f"\n"
                f"integer function set_{name}(index_of_the_particle, i, {name})\n"
                f"    implicit none\n"
                f"    integer, intent(in):: index_of_the_particle\n"
                f"    integer, intent(in):: i\n"
                f"    {f_dtype}, intent(in):: {name}\n"
                f"    set_{name} = -1\n"
                f"    if ((i > 0) .and. (i < 6)) then\n"
                f"        GenecStar%{name}(i) = {name}\n"
                f"        set_{name} = 0\n"
                f"    end if\n"
                f"end function set_{name}\n"
                f"\n"
            )

def f_new_stellar_model(outfilename):
    f_outfile = f"{outfilename}.f90"
    with open(f_outfile, "a") as f_out:
        f_out.write(
            "function new_stellar_model(&\n"
            "      index_of_the_star,&\n"
        )
        line_length = line_length_start = 6
        first = True
        max_linelength = 100 - 2  # 2 for ",&" buffer
        for parameter in {**SCALAR_SETTERS, **GENEC_STAR_STRUCTURE}.items():
            name = parameter[0]
            f_dtype = get_f_dtype(parameter[1][0])
            if line_length + 1 + len(name) <= max_linelength:
                if first:
                    name = f"      {name}"
                    first = False
                else:
                    name = f",{name}"
                line_length += len(name)
                f_out.write(f"{name}")
            else:
                line_length = line_length_start + len(name)
                f_out.write(
                    f",&\n"
                    f"      {name}"
                )
        f_out.write(
            "&\n"
            "      )\n"
            "    implicit none\n"
            "    integer:: index_of_the_star\n"
        )
        line_length = line_length_start = 13
        previous_f_dtype = None
        for parameter in {**SCALAR_SETTERS, **GENEC_STAR_STRUCTURE}.items():
            name = parameter[0]
            f_dtype = get_f_dtype(parameter[1][0])
            if f_dtype != previous_f_dtype:
                f_out.write(
                    f"\n"
                    f"    {f_dtype}, intent(in):: &\n"
                )
                line_length += len(name)
                name = line_length_start * " " + name
                previous_f_dtype = f_dtype
                f_out.write(f"{name}")
            else:
                if line_length + 1 + len(name) <= max_linelength:
                    name = f",{name}"
                    line_length += len(name)
                    f_out.write(f"{name}")
                else:
                    name = line_length_start * " " + name
                    line_length = len(name)
                    f_out.write(
                        f",&\n"
                        f"{name}"
                    )
        f_out.write(
            "\n"
            "    integer:: new_stellar_model\n"
            "\n"
        )
        for parameter in SCALAR_SETTERS.items():
            name = parameter[0]
            f_dtype = get_f_dtype(parameter[1][0])
            if f_dtype == "logical":
                f_out.write(
                    #f"    if (GenecStar%{name} .neqv. {name}) then\n"
                    #f'        write(*,*) "{name}:", GenecStar%{name}, {name}\n'
                    #f"    endif\n"
                    f"    GenecStar%{name} = {name}\n"
                )
            else:
                f_out.write(
                    #f"    if (GenecStar%{name} /= {name}) then\n"
                    #f'        write(*,*) "{name}:", GenecStar%{name}, {name}\n'
                    #f"    endif\n"
                    f"    GenecStar%{name} = {name}\n"
                )
        for parameter in GENEC_STAR_STRUCTURE.items():
            name = parameter[0]
            f_out.write(
                f"    GenecStar%{name} = 0.d0\n"
                f"    GenecStar%{name} = {name}\n"
            )
        f_out.write(
            "\n"
            "    new_stellar_model = 0\n"
            "end function new_stellar_model\n"
            "\n"
        )


def py_new_stellar_model(outfilename):
    py_outfile = f"{outfilename}.py"
    with open(py_outfile, "a") as py_out:
        py_out.write(
            "    def new_particle_from_model(\n"
            "        self, internal_structure, current_age=0 | units.julianyr, key=None\n"
            "    ):\n"
            "        self.new_stellar_model(\n"
        )
        for parameter in ALL_SETTERS.items():
            name = parameter[0]
            py_out.write(
                f"            internal_structure['{name}'],\n"
            )
        py_out.write(
            "        )\n"
            "        return\n"
        )


def py_header(outfilename):
    py_outfile = f"{outfilename}.py"
    with open(py_outfile, "a") as py_out:
        py_out.write(
            PY_INTERFACE_HEADER
        )
    return


def get_f_dtype(dtype):
    if dtype in ["i", "int32"]:
        f_dtype = "integer"
    elif dtype in ["d", "float64"]:
        f_dtype = "real(kindreal)"
    elif dtype in ["b", "bool"]:
        f_dtype = "logical"
    elif dtype in ["s", "string"]:
        f_dtype = "character(256)"
    else:
        print(f"wrong dype: {dtype}")
        sys.exit()
    return f_dtype

def main():
    outfilename = "new_interface"
    # py_header(outfilename)
    # star_parameters(outfilename)
    # star_structure(outfilename)
    # star_structure_extra(outfilename)
    # star_structure_derived(outfilename)
    # netdef_scalars(outfilename)
    netdef_arrays(outfilename)
    netalu_arrays(outfilename)
    # py_new_stellar_model(outfilename)
    # f_new_stellar_model(outfilename)


if __name__ == "__main__":
    # for i in ALL_SETTERS.items():
    #     i = i[0]
    #     print(f"        {i} = Star%{i}")
    main()
