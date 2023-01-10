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

# Parameters (but individual to each star)
GENEC_STAR_PARAMETERS = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (optional)]
    'initialised': ['bool', '', "True if the star is an intialised model"],
    'starname': ['string', '', "Name of the star"],
    'nwseq': ['int32', '', ""],
    'modanf': ['int32', '', ""],
    'nzmod': ['int32', '', ""],
    'end_at_phase': ['int32', '', "Stop if this phase is reached"],
    'end_at_model': ['int32', '', "Stop if this model number is reached"],
    'irot': ['int32', '', ""],
    'isol': ['int32', '', ""],
    'imagn': ['int32', '', ""],
    'ialflu': ['int32', '', ""],
    'ianiso': ['int32', '', ""],
    'ipop3': ['int32', '', ""],
    'ibasnet': ['int32', '', ""],
    'phase': ['int32', '', ""],
    'var_rates': ['bool', '', ""],
    'bintide': ['bool', '', ""],
    'binm2': ['float64', '', ""],
    'periodini': ['float64', '', ""],
    'const_per': ['bool', '', ""],
    'iprezams': ['int32', '', ""],
    'zinit': ['float64', '', ""],
    'zsol': ['float64', '', ""],
    'z': ['float64', '', ""],
    'iopac': ['int32', '', ""],
    'ikappa': ['int32', '', ""],
    'idiff': ['int32', '', ""],
    'iadvec': ['int32', '', ""],
    'istati': ['int32', '', ""],
    'icoeff': ['int32', '', ""],
    'fenerg': ['float64', '', ""],
    'richac': ['float64', '', ""],
    'igamma': ['int32', '', ""],
    'frein': ['float64', '', ""],
    'K_Kawaler': ['float64', '', ""],
    'Omega_saturation': ['float64', '', ""],
    'rapcrilim': ['float64', '', ""],
    'vwant': ['float64', '', ""],
    'xfom': ['float64', '', ""],
    'omega': ['float64', '', ""],
    'xdial': ['float64', '', ""],
    'idialo': ['int32', '', ""],
    'idialu': ['int32', '', ""],
    'Add_Flux': ['bool', '', ""],
    'diff_only': ['bool', '', ""],
    'B_initial': ['float64', '', ""],
    'add_diff': ['float64', '', ""],
    'n_mag': ['int32', '', ""],
    'alpha_F': ['float64', '', ""],
    'nsmooth': ['int32', '', ""],
    'qminsmooth': ['bool', '', ""],
    'imloss': ['int32', '', ""],
    'fmlos': ['float64', '', ""],
    'ifitm': ['int32', '', ""],
    'fitm': ['float64', '', ""],
    'fitmi': ['float64', '', ""],
    'deltal': ['float64', '', ""],
    'deltat': ['float64', '', ""],
    'nndr': ['int32', '', ""],
    'RSG_Mdot': ['int32', '', ""],
    'SupraEddMdot': ['bool', '', ""],
    'Be_mdotfrac': ['float64', '', ""],
    'start_mdot': ['float64', '', ""],
    'iledou': ['int32', '', ""],
    'idifcon': ['int32', '', ""],
    'iover': ['int32', '', ""],
    'elph': ['float64', '', ""],
    'my': ['int32', '', ""],
    'dovhp': ['float64', '', ""],
    'iunder': ['int32', '', ""],
    'dunder': ['float64', '', ""],
    'gkorm': ['float64', '', ""],
    'alph': ['float64', '', ""],
    'agdr': ['float64', '', ""],
    'faktor': ['float64', '', ""],
    'dgrp': ['float64', '', ""],
    'dgrl': ['float64', '', ""],
    'dgry': ['float64', '', ""],
    'dgrc': ['float64', '', ""],
    'dgro': ['float64', '', ""],
    'dgr20': ['float64', '', ""],
    'nbchx': ['int32', '', ""],
    'nrband': ['int32', '', ""],
    'xcn': ['float64', '', ""],
    'islow': ['int32', '', ""],
    'icncst': ['int32', '', ""],
    'tauH_fit': ['int32', '', ""],
    'display_plot': ['bool', '', ""],
    'iauto': ['int32', '', ""],
    'iprn': ['int32', '', ""],
    'iout': ['int32', '', ""],
    'itmin': ['int32', '', ""],
    'xyfiles': ['bool', '', ""],
    'idebug': ['int32', '', ""],
    'itests': ['int32', '', ""],
    'verbose': ['bool', '', ""],
    'stop_deg': ['bool', '', ""],
    'n_snap': ['int32', '', "number of steps between snapshots [0]"],
}

# Stellar properties (but global for the star)
GENEC_STAR_PROPERTIES = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (optional)]
    'm': ['int32', '', "number of zones", "n_zones"],
    'gms': ['float64', 'MSun', "total mass", 'mass'],
    'alter': ['float64', 'julianyr', "stellar age", 'age'],
    'gls': ['float64', 'LSun', "", 'luminosity'],
    'teff': ['float64', 'K', "effective temperature", 'temperature'],
    'glsv': ['float64', 'LSun', "previous luminosity"],
    'teffv': ['float64', 'K', "previous effective temperature"],
    'dzeitj': ['float64', 'julianyr', "time step", 'time_step'],
    'dzeit': ['float64', 's', "time step"],
    'dzeitv': ['float64', 's', "previous time step"],
    'xmini': ['float64', 'MSun', "", 'initial_mass'],
    'summas': ['float64', 'MSun', "total mass"],
    'ab': ['float64', '', ""],
    'dm_lost': ['float64', 'MSun', "total mass lost"],
}

# Structural properties (m layers)
GENEC_STAR_STRUCTURE = {
    # 'GENEC name: [dtype, unit, description, AMUSE name (optional)]
    'q': ['float64', '', ""],
    'p': ['float64', '', ""],
    't': ['float64', '', ""],
    'r': ['float64', '', ""],
    's': ['float64', '', ""],
    'x': ['float64', '', "H fraction"],
    'y3': ['float64', '', "He3 fraction"],
    'y': ['float64', '', "He fraction"],
    'xc12': ['float64', '', "C12 fraction"],
    'xc13': ['float64', '', "C13 fraction"],
    'xn14': ['float64', '', "N14 fraction"],
    'xn15': ['float64', '', ""],
    'xo16': ['float64', '', ""],
    'xo17': ['float64', '', ""],
    'xo18': ['float64', '', ""],
    'xne20': ['float64', '', ""],
    'xne22': ['float64', '', ""],
    'xmg24': ['float64', '', ""],
    'xmg25': ['float64', '', ""],
    'xmg26': ['float64', '', ""],
    'xf19': ['float64', '', ""],
    'xne21': ['float64', '', ""],
    'xna23': ['float64', '', ""],
    'xal27': ['float64', '', ""],
    'xsi28': ['float64', '', ""],
    'xc14': ['float64', '', ""],
    'xf18': ['float64', '', ""],
    'xal26': ['float64', '', ""],
    'xneut': ['float64', '', "Neutron fraction"],
    'xprot': ['float64', '', "Proton fraction"],
    'omegi': ['float64', '', "Rotation"],
    'xbid': ['float64', '', ""],
    'xbid1': ['float64', '', ""],
    'vp': ['float64', '', ""],
    'vt': ['float64', '', ""],
    'vr': ['float64', '', ""],
    'vs': ['float64', '', ""],
    'vx': ['float64', '', ""],
    'vy': ['float64', '', ""],
    'vy3': ['float64', '', ""],
    'vxc12': ['float64', '', ""],
    'vxc13': ['float64', '', ""],
    'vxn14': ['float64', '', ""],
    'vxn15': ['float64', '', ""],
    'vxo16': ['float64', '', ""],
    'vxo17': ['float64', '', ""],
    'vxo18': ['float64', '', ""],
    'vxne20': ['float64', '', ""],
    'vxne22': ['float64', '', ""],
    'vxmg24': ['float64', '', ""],
    'vxmg25': ['float64', '', ""],
    'vxmg26': ['float64', '', ""],
    'vxf19': ['float64', '', ""],
    'vxne21': ['float64', '', ""],
    'vxna23': ['float64', '', ""],
    'vxal27': ['float64', '', ""],
    'vxsi28': ['float64', '', ""],
    'vxc14': ['float64', '', ""],
    'vxf18': ['float64', '', ""],
    'vxal26g': ['float64', '', ""],
    'vxneut': ['float64', '', ""],
    'vxprot': ['float64', '', ""],
    'vomegi': ['float64', '', ""],
    'vxbid': ['float64', '', ""],
    'vxbid1': ['float64', '', ""],
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

    @legacy_function
    def finalize_stellar_model():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @remote_function
    def set_genec_path(genec_path='s'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_initialised(index_of_the_particle='i'):
        returns (initialised='bool')

    @remote_function(can_handle_array=True)
    def set_initialised(index_of_the_particle='i', initialised='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_starname(index_of_the_particle='i'):
        returns (starname='string')

    @remote_function(can_handle_array=True)
    def set_starname(index_of_the_particle='i', starname='string'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nwseq(index_of_the_particle='i'):
        returns (nwseq='int32')

    @remote_function(can_handle_array=True)
    def set_nwseq(index_of_the_particle='i', nwseq='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_modanf(index_of_the_particle='i'):
        returns (modanf='int32')

    @remote_function(can_handle_array=True)
    def set_modanf(index_of_the_particle='i', modanf='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nzmod(index_of_the_particle='i'):
        returns (nzmod='int32')

    @remote_function(can_handle_array=True)
    def set_nzmod(index_of_the_particle='i', nzmod='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_end_at_phase(index_of_the_particle='i'):
        returns (end_at_phase='int32')

    @remote_function(can_handle_array=True)
    def set_end_at_phase(index_of_the_particle='i', end_at_phase='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_end_at_model(index_of_the_particle='i'):
        returns (end_at_model='int32')

    @remote_function(can_handle_array=True)
    def set_end_at_model(index_of_the_particle='i', end_at_model='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_irot(index_of_the_particle='i'):
        returns (irot='int32')

    @remote_function(can_handle_array=True)
    def set_irot(index_of_the_particle='i', irot='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_isol(index_of_the_particle='i'):
        returns (isol='int32')

    @remote_function(can_handle_array=True)
    def set_isol(index_of_the_particle='i', isol='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_imagn(index_of_the_particle='i'):
        returns (imagn='int32')

    @remote_function(can_handle_array=True)
    def set_imagn(index_of_the_particle='i', imagn='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ialflu(index_of_the_particle='i'):
        returns (ialflu='int32')

    @remote_function(can_handle_array=True)
    def set_ialflu(index_of_the_particle='i', ialflu='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ianiso(index_of_the_particle='i'):
        returns (ianiso='int32')

    @remote_function(can_handle_array=True)
    def set_ianiso(index_of_the_particle='i', ianiso='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ipop3(index_of_the_particle='i'):
        returns (ipop3='int32')

    @remote_function(can_handle_array=True)
    def set_ipop3(index_of_the_particle='i', ipop3='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ibasnet(index_of_the_particle='i'):
        returns (ibasnet='int32')

    @remote_function(can_handle_array=True)
    def set_ibasnet(index_of_the_particle='i', ibasnet='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_phase(index_of_the_particle='i'):
        returns (phase='int32')

    @remote_function(can_handle_array=True)
    def set_phase(index_of_the_particle='i', phase='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_var_rates(index_of_the_particle='i'):
        returns (var_rates='bool')

    @remote_function(can_handle_array=True)
    def set_var_rates(index_of_the_particle='i', var_rates='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_bintide(index_of_the_particle='i'):
        returns (bintide='bool')

    @remote_function(can_handle_array=True)
    def set_bintide(index_of_the_particle='i', bintide='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_binm2(index_of_the_particle='i'):
        returns (binm2='float64')

    @remote_function(can_handle_array=True)
    def set_binm2(index_of_the_particle='i', binm2='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_periodini(index_of_the_particle='i'):
        returns (periodini='float64')

    @remote_function(can_handle_array=True)
    def set_periodini(index_of_the_particle='i', periodini='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_const_per(index_of_the_particle='i'):
        returns (const_per='bool')

    @remote_function(can_handle_array=True)
    def set_const_per(index_of_the_particle='i', const_per='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iprezams(index_of_the_particle='i'):
        returns (iprezams='int32')

    @remote_function(can_handle_array=True)
    def set_iprezams(index_of_the_particle='i', iprezams='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_zinit(index_of_the_particle='i'):
        returns (zinit='float64')

    @remote_function(can_handle_array=True)
    def set_zinit(index_of_the_particle='i', zinit='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_zsol(index_of_the_particle='i'):
        returns (zsol='float64')

    @remote_function(can_handle_array=True)
    def set_zsol(index_of_the_particle='i', zsol='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_z(index_of_the_particle='i'):
        returns (z='float64')

    @remote_function(can_handle_array=True)
    def set_z(index_of_the_particle='i', z='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iopac(index_of_the_particle='i'):
        returns (iopac='int32')

    @remote_function(can_handle_array=True)
    def set_iopac(index_of_the_particle='i', iopac='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ikappa(index_of_the_particle='i'):
        returns (ikappa='int32')

    @remote_function(can_handle_array=True)
    def set_ikappa(index_of_the_particle='i', ikappa='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idiff(index_of_the_particle='i'):
        returns (idiff='int32')

    @remote_function(can_handle_array=True)
    def set_idiff(index_of_the_particle='i', idiff='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iadvec(index_of_the_particle='i'):
        returns (iadvec='int32')

    @remote_function(can_handle_array=True)
    def set_iadvec(index_of_the_particle='i', iadvec='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_istati(index_of_the_particle='i'):
        returns (istati='int32')

    @remote_function(can_handle_array=True)
    def set_istati(index_of_the_particle='i', istati='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_icoeff(index_of_the_particle='i'):
        returns (icoeff='int32')

    @remote_function(can_handle_array=True)
    def set_icoeff(index_of_the_particle='i', icoeff='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fenerg(index_of_the_particle='i'):
        returns (fenerg='float64')

    @remote_function(can_handle_array=True)
    def set_fenerg(index_of_the_particle='i', fenerg='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_richac(index_of_the_particle='i'):
        returns (richac='float64')

    @remote_function(can_handle_array=True)
    def set_richac(index_of_the_particle='i', richac='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_igamma(index_of_the_particle='i'):
        returns (igamma='int32')

    @remote_function(can_handle_array=True)
    def set_igamma(index_of_the_particle='i', igamma='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_frein(index_of_the_particle='i'):
        returns (frein='float64')

    @remote_function(can_handle_array=True)
    def set_frein(index_of_the_particle='i', frein='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_K_Kawaler(index_of_the_particle='i'):
        returns (K_Kawaler='float64')

    @remote_function(can_handle_array=True)
    def set_K_Kawaler(index_of_the_particle='i', K_Kawaler='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Omega_saturation(index_of_the_particle='i'):
        returns (Omega_saturation='float64')

    @remote_function(can_handle_array=True)
    def set_Omega_saturation(index_of_the_particle='i', Omega_saturation='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_rapcrilim(index_of_the_particle='i'):
        returns (rapcrilim='float64')

    @remote_function(can_handle_array=True)
    def set_rapcrilim(index_of_the_particle='i', rapcrilim='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_vwant(index_of_the_particle='i'):
        returns (vwant='float64')

    @remote_function(can_handle_array=True)
    def set_vwant(index_of_the_particle='i', vwant='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xfom(index_of_the_particle='i'):
        returns (xfom='float64')

    @remote_function(can_handle_array=True)
    def set_xfom(index_of_the_particle='i', xfom='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_omega(index_of_the_particle='i'):
        returns (omega='float64')

    @remote_function(can_handle_array=True)
    def set_omega(index_of_the_particle='i', omega='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xdial(index_of_the_particle='i'):
        returns (xdial='float64')

    @remote_function(can_handle_array=True)
    def set_xdial(index_of_the_particle='i', xdial='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idialo(index_of_the_particle='i'):
        returns (idialo='int32')

    @remote_function(can_handle_array=True)
    def set_idialo(index_of_the_particle='i', idialo='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idialu(index_of_the_particle='i'):
        returns (idialu='int32')

    @remote_function(can_handle_array=True)
    def set_idialu(index_of_the_particle='i', idialu='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Add_Flux(index_of_the_particle='i'):
        returns (Add_Flux='bool')

    @remote_function(can_handle_array=True)
    def set_Add_Flux(index_of_the_particle='i', Add_Flux='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_diff_only(index_of_the_particle='i'):
        returns (diff_only='bool')

    @remote_function(can_handle_array=True)
    def set_diff_only(index_of_the_particle='i', diff_only='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_B_initial(index_of_the_particle='i'):
        returns (B_initial='float64')

    @remote_function(can_handle_array=True)
    def set_B_initial(index_of_the_particle='i', B_initial='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_add_diff(index_of_the_particle='i'):
        returns (add_diff='float64')

    @remote_function(can_handle_array=True)
    def set_add_diff(index_of_the_particle='i', add_diff='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_n_mag(index_of_the_particle='i'):
        returns (n_mag='int32')

    @remote_function(can_handle_array=True)
    def set_n_mag(index_of_the_particle='i', n_mag='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alpha_F(index_of_the_particle='i'):
        returns (alpha_F='float64')

    @remote_function(can_handle_array=True)
    def set_alpha_F(index_of_the_particle='i', alpha_F='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nsmooth(index_of_the_particle='i'):
        returns (nsmooth='int32')

    @remote_function(can_handle_array=True)
    def set_nsmooth(index_of_the_particle='i', nsmooth='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_qminsmooth(index_of_the_particle='i'):
        returns (qminsmooth='bool')

    @remote_function(can_handle_array=True)
    def set_qminsmooth(index_of_the_particle='i', qminsmooth='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_imloss(index_of_the_particle='i'):
        returns (imloss='int32')

    @remote_function(can_handle_array=True)
    def set_imloss(index_of_the_particle='i', imloss='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fmlos(index_of_the_particle='i'):
        returns (fmlos='float64')

    @remote_function(can_handle_array=True)
    def set_fmlos(index_of_the_particle='i', fmlos='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ifitm(index_of_the_particle='i'):
        returns (ifitm='int32')

    @remote_function(can_handle_array=True)
    def set_ifitm(index_of_the_particle='i', ifitm='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fitm(index_of_the_particle='i'):
        returns (fitm='float64')

    @remote_function(can_handle_array=True)
    def set_fitm(index_of_the_particle='i', fitm='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_fitmi(index_of_the_particle='i'):
        returns (fitmi='float64')

    @remote_function(can_handle_array=True)
    def set_fitmi(index_of_the_particle='i', fitmi='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_deltal(index_of_the_particle='i'):
        returns (deltal='float64')

    @remote_function(can_handle_array=True)
    def set_deltal(index_of_the_particle='i', deltal='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_deltat(index_of_the_particle='i'):
        returns (deltat='float64')

    @remote_function(can_handle_array=True)
    def set_deltat(index_of_the_particle='i', deltat='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nndr(index_of_the_particle='i'):
        returns (nndr='int32')

    @remote_function(can_handle_array=True)
    def set_nndr(index_of_the_particle='i', nndr='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_RSG_Mdot(index_of_the_particle='i'):
        returns (RSG_Mdot='int32')

    @remote_function(can_handle_array=True)
    def set_RSG_Mdot(index_of_the_particle='i', RSG_Mdot='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_SupraEddMdot(index_of_the_particle='i'):
        returns (SupraEddMdot='bool')

    @remote_function(can_handle_array=True)
    def set_SupraEddMdot(index_of_the_particle='i', SupraEddMdot='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_Be_mdotfrac(index_of_the_particle='i'):
        returns (Be_mdotfrac='float64')

    @remote_function(can_handle_array=True)
    def set_Be_mdotfrac(index_of_the_particle='i', Be_mdotfrac='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_start_mdot(index_of_the_particle='i'):
        returns (start_mdot='float64')

    @remote_function(can_handle_array=True)
    def set_start_mdot(index_of_the_particle='i', start_mdot='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iledou(index_of_the_particle='i'):
        returns (iledou='int32')

    @remote_function(can_handle_array=True)
    def set_iledou(index_of_the_particle='i', iledou='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idifcon(index_of_the_particle='i'):
        returns (idifcon='int32')

    @remote_function(can_handle_array=True)
    def set_idifcon(index_of_the_particle='i', idifcon='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iover(index_of_the_particle='i'):
        returns (iover='int32')

    @remote_function(can_handle_array=True)
    def set_iover(index_of_the_particle='i', iover='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_elph(index_of_the_particle='i'):
        returns (elph='float64')

    @remote_function(can_handle_array=True)
    def set_elph(index_of_the_particle='i', elph='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_my(index_of_the_particle='i'):
        returns (my='int32')

    @remote_function(can_handle_array=True)
    def set_my(index_of_the_particle='i', my='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dovhp(index_of_the_particle='i'):
        returns (dovhp='float64')

    @remote_function(can_handle_array=True)
    def set_dovhp(index_of_the_particle='i', dovhp='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iunder(index_of_the_particle='i'):
        returns (iunder='int32')

    @remote_function(can_handle_array=True)
    def set_iunder(index_of_the_particle='i', iunder='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dunder(index_of_the_particle='i'):
        returns (dunder='float64')

    @remote_function(can_handle_array=True)
    def set_dunder(index_of_the_particle='i', dunder='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gkorm(index_of_the_particle='i'):
        returns (gkorm='float64')

    @remote_function(can_handle_array=True)
    def set_gkorm(index_of_the_particle='i', gkorm='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alph(index_of_the_particle='i'):
        returns (alph='float64')

    @remote_function(can_handle_array=True)
    def set_alph(index_of_the_particle='i', alph='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_agdr(index_of_the_particle='i'):
        returns (agdr='float64')

    @remote_function(can_handle_array=True)
    def set_agdr(index_of_the_particle='i', agdr='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_faktor(index_of_the_particle='i'):
        returns (faktor='float64')

    @remote_function(can_handle_array=True)
    def set_faktor(index_of_the_particle='i', faktor='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrp(index_of_the_particle='i'):
        returns (dgrp='float64')

    @remote_function(can_handle_array=True)
    def set_dgrp(index_of_the_particle='i', dgrp='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrl(index_of_the_particle='i'):
        returns (dgrl='float64')

    @remote_function(can_handle_array=True)
    def set_dgrl(index_of_the_particle='i', dgrl='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgry(index_of_the_particle='i'):
        returns (dgry='float64')

    @remote_function(can_handle_array=True)
    def set_dgry(index_of_the_particle='i', dgry='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgrc(index_of_the_particle='i'):
        returns (dgrc='float64')

    @remote_function(can_handle_array=True)
    def set_dgrc(index_of_the_particle='i', dgrc='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgro(index_of_the_particle='i'):
        returns (dgro='float64')

    @remote_function(can_handle_array=True)
    def set_dgro(index_of_the_particle='i', dgro='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dgr20(index_of_the_particle='i'):
        returns (dgr20='float64')

    @remote_function(can_handle_array=True)
    def set_dgr20(index_of_the_particle='i', dgr20='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nbchx(index_of_the_particle='i'):
        returns (nbchx='int32')

    @remote_function(can_handle_array=True)
    def set_nbchx(index_of_the_particle='i', nbchx='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_nrband(index_of_the_particle='i'):
        returns (nrband='int32')

    @remote_function(can_handle_array=True)
    def set_nrband(index_of_the_particle='i', nrband='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xcn(index_of_the_particle='i'):
        returns (xcn='float64')

    @remote_function(can_handle_array=True)
    def set_xcn(index_of_the_particle='i', xcn='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_islow(index_of_the_particle='i'):
        returns (islow='int32')

    @remote_function(can_handle_array=True)
    def set_islow(index_of_the_particle='i', islow='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_icncst(index_of_the_particle='i'):
        returns (icncst='int32')

    @remote_function(can_handle_array=True)
    def set_icncst(index_of_the_particle='i', icncst='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_tauH_fit(index_of_the_particle='i'):
        returns (tauH_fit='int32')

    @remote_function(can_handle_array=True)
    def set_tauH_fit(index_of_the_particle='i', tauH_fit='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_display_plot(index_of_the_particle='i'):
        returns (display_plot='bool')

    @remote_function(can_handle_array=True)
    def set_display_plot(index_of_the_particle='i', display_plot='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iauto(index_of_the_particle='i'):
        returns (iauto='int32')

    @remote_function(can_handle_array=True)
    def set_iauto(index_of_the_particle='i', iauto='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iprn(index_of_the_particle='i'):
        returns (iprn='int32')

    @remote_function(can_handle_array=True)
    def set_iprn(index_of_the_particle='i', iprn='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_iout(index_of_the_particle='i'):
        returns (iout='int32')

    @remote_function(can_handle_array=True)
    def set_iout(index_of_the_particle='i', iout='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_itmin(index_of_the_particle='i'):
        returns (itmin='int32')

    @remote_function(can_handle_array=True)
    def set_itmin(index_of_the_particle='i', itmin='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xyfiles(index_of_the_particle='i'):
        returns (xyfiles='bool')

    @remote_function(can_handle_array=True)
    def set_xyfiles(index_of_the_particle='i', xyfiles='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_idebug(index_of_the_particle='i'):
        returns (idebug='int32')

    @remote_function(can_handle_array=True)
    def set_idebug(index_of_the_particle='i', idebug='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_itests(index_of_the_particle='i'):
        returns (itests='int32')

    @remote_function(can_handle_array=True)
    def set_itests(index_of_the_particle='i', itests='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_verbose(index_of_the_particle='i'):
        returns (verbose='bool')

    @remote_function(can_handle_array=True)
    def set_verbose(index_of_the_particle='i', verbose='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_stop_deg(index_of_the_particle='i'):
        returns (stop_deg='bool')

    @remote_function(can_handle_array=True)
    def set_stop_deg(index_of_the_particle='i', stop_deg='bool'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_n_snap(index_of_the_particle='i'):
        returns (n_snap='int32')

    @remote_function(can_handle_array=True)
    def set_n_snap(index_of_the_particle='i', n_snap='int32'):
        returns ()

    # end Parameters

    # begin Properties
    @remote_function(can_handle_array=True)
    def get_m(index_of_the_particle='i'):
        returns (m='int32')

    @remote_function(can_handle_array=True)
    def set_m(index_of_the_particle='i', m='int32'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gms(index_of_the_particle='i'):
        returns (gms='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_gms(index_of_the_particle='i', gms='float64' | units.MSun):
        returns ()

    @remote_function(can_handle_array=True)
    def get_alter(index_of_the_particle='i'):
        returns (alter='float64' | units.julianyr)

    @remote_function(can_handle_array=True)
    def set_alter(index_of_the_particle='i', alter='float64' | units.julianyr):
        returns ()

    @remote_function(can_handle_array=True)
    def get_gls(index_of_the_particle='i'):
        returns (gls='float64' | units.LSun)

    @remote_function(can_handle_array=True)
    def set_gls(index_of_the_particle='i', gls='float64' | units.LSun):
        returns ()

    @remote_function(can_handle_array=True)
    def get_teff(index_of_the_particle='i'):
        returns (teff='float64' | units.K)

    @remote_function(can_handle_array=True)
    def set_teff(index_of_the_particle='i', teff='float64' | units.K):
        returns ()

    @remote_function(can_handle_array=True)
    def get_glsv(index_of_the_particle='i'):
        returns (glsv='float64' | units.LSun)

    @remote_function(can_handle_array=True)
    def set_glsv(index_of_the_particle='i', glsv='float64' | units.LSun):
        returns ()

    @remote_function(can_handle_array=True)
    def get_teffv(index_of_the_particle='i'):
        returns (teffv='float64' | units.K)

    @remote_function(can_handle_array=True)
    def set_teffv(index_of_the_particle='i', teffv='float64' | units.K):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeitj(index_of_the_particle='i'):
        returns (dzeitj='float64' | units.julianyr)

    @remote_function(can_handle_array=True)
    def set_dzeitj(index_of_the_particle='i', dzeitj='float64' | units.julianyr):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeit(index_of_the_particle='i'):
        returns (dzeit='float64' | units.s)

    @remote_function(can_handle_array=True)
    def set_dzeit(index_of_the_particle='i', dzeit='float64' | units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dzeitv(index_of_the_particle='i'):
        returns (dzeitv='float64' | units.s)

    @remote_function(can_handle_array=True)
    def set_dzeitv(index_of_the_particle='i', dzeitv='float64' | units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_xmini(index_of_the_particle='i'):
        returns (xmini='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_xmini(index_of_the_particle='i', xmini='float64' | units.MSun):
        returns ()

    @remote_function(can_handle_array=True)
    def get_summas(index_of_the_particle='i'):
        returns (summas='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_summas(index_of_the_particle='i', summas='float64' | units.MSun):
        returns ()

    @remote_function(can_handle_array=True)
    def get_ab(index_of_the_particle='i'):
        returns (ab='float64')

    @remote_function(can_handle_array=True)
    def set_ab(index_of_the_particle='i', ab='float64'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_dm_lost(index_of_the_particle='i'):
        returns (dm_lost='float64' | units.MSun)

    @remote_function(can_handle_array=True)
    def set_dm_lost(index_of_the_particle='i', dm_lost='float64' | units.MSun):
        returns ()

    # end Properties

    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.OUT,
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
            'starname', dtype='string', direction=function.IN,
            default='AmuseStar', description="The star's name")
        # function.addParameter(
        #     'age_tag', dtype='float64', direction=function.IN,
        #     description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
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
        function.can_handle_array = False
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.OUT,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        for parameter in {
            **GENEC_STAR_PARAMETERS,
            **GENEC_STAR_PROPERTIES,
            **GENEC_STAR_STRUCTURE
        }.items():
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

        function.result_type = 'int32'
        return function

    @legacy_function
    def get_stellar_model():
        """
        Get an existing model
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        for parameter in {
            **GENEC_STAR_PARAMETERS,
            **GENEC_STAR_PROPERTIES,
            **GENEC_STAR_STRUCTURE
        }.items():
            function.addParameter(
                parameter[0],
                dtype=parameter[1][0],
                unit=parameter[1][1],
                description=parameter[1][2],
                direction=function.OUT,
            )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_firstlast_zone():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first', dtype='int32', direction=function.OUT)
        function.addParameter('last', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_surface_velocity():
        "Retrieve the surface velocity of the star."
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of"
        )
        function.addParameter(
            'surface_velocity', dtype='float64', direction=function.OUT
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_luminosity_at_zone():
        """
        Retrieve the luminosity at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of"
        )
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of"
        )
        function.addParameter(
            'lum_i', dtype='float64', direction=function.OUT,
            unit=units.LSun,
            description=(
                "The luminosity at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_of_species():
        """
        Retrieve the mass number of the chemical abundance variable of the
        star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'species', dtype='int32', direction=function.IN,
            description="The species of the star to get the mass number of")
        function.addParameter(
            'species_mass', dtype='float64', direction=function.OUT,
            description=(
                "The mass number of the chemical abundance variable of "
                "the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_at_zone():
        """
        Retrieve eps at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps', dtype='float64', direction=function.OUT,
            description=(
                "eps at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_epsy_at_zone():
        """
        Retrieve epsy at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'epsy', dtype='float64', direction=function.OUT,
            description=(
                "epsy at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_c_adv_at_zone():
        """
        Retrieve eps_c_adv at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_c_adv', dtype='float64', direction=function.OUT,
            description=(
                "eps_c_adv at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_ne_adv_at_zone():
        """
        Retrieve eps_ne_adv at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_ne_adv', dtype='float64', direction=function.OUT,
            description=(
                "eps_ne_adv at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_o_adv_at_zone():
        """
        Retrieve eps_o_adv at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_o_adv', dtype='float64', direction=function.OUT,
            description=(
                "eps_o_adv at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_si_adv_at_zone():
        """
        Retrieve eps_si_adv at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_si_adv', dtype='float64', direction=function.OUT,
            description=(
                "eps_si_adv at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_grav_at_zone():
        """
        Retrieve eps_grav at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_grav', dtype='float64', direction=function.OUT,
            description=(
                "eps_grav at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_eps_nu_at_zone():
        """
        Retrieve eps_nu at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'eps_nu', dtype='float64', direction=function.OUT,
            description=(
                "eps_nu at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_nabla_rad_at_zone():
        """
        Retrieve nabla_rad at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'nabla_rad', dtype='float64', direction=function.OUT,
            description=(
                "nabla_rad at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_nabla_ad_at_zone():
        """
        Retrieve nabla_rad at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'nabla_ad', dtype='float64', direction=function.OUT,
            description=(
                "nabla_ad at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_nabla_mu_at_zone():
        """
        Retrieve nabla_rad at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'nabla_mu', dtype='float64', direction=function.OUT,
            description=(
                "nabla_mu at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_fraction_at_zone():
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'dq_i', dtype='float64', direction=function.OUT,
            description=(
                "The mass fraction at the specified zone/mesh-cell of "
                "the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_omegi_at_zone():
        """
        Retrieve the rotation rate at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'omegi_i', dtype='float64', direction=function.OUT,
            description=(
                "The rotation rate at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_fraction_of_h_at_zone():
        """
        Retrieve the fractional abundance of h at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_he3_at_zone():
        """
        Retrieve the fractional abundance of he3 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_he_at_zone():
        """
        Retrieve the fractional abundance of he at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c12_at_zone():
        """
        Retrieve the fractional abundance of c12 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c13_at_zone():
        """
        Retrieve the fractional abundance of c13 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_n14_at_zone():
        """
        Retrieve the fractional abundance of n14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_n15_at_zone():
        """
        Retrieve the fractional abundance of n15 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o16_at_zone():
        """
        Retrieve the fractional abundance of o16 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o17_at_zone():
        """
        Retrieve the fractional abundance of o17 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_o18_at_zone():
        """
        Retrieve the fractional abundance of o18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne20_at_zone():
        """
        Retrieve the fractional abundance of ne20 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne22_at_zone():
        """
        Retrieve the fractional abundance of ne22 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg24_at_zone():
        """
        Retrieve the fractional abundance of mg24 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg25_at_zone():
        """
        Retrieve the fractional abundance of mg25 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_mg26_at_zone():
        """
        Retrieve the fractional abundance of mg26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_c14_at_zone():
        """
        Retrieve the fractional abundance of c14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_f18_at_zone():
        """
        Retrieve the fractional abundance of f18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_f19_at_zone():
        """
        Retrieve the fractional abundance of f19 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_ne21_at_zone():
        """
        Retrieve the fractional abundance of ne21 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_na23_at_zone():
        """
        Retrieve the fractional abundance of na23 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_al26_at_zone():
        """
        Retrieve the fractional abundance of al26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_al27_at_zone():
        """
        Retrieve the fractional abundance of al27 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_si28_at_zone():
        """
        Retrieve the fractional abundance of si28 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_neut_at_zone():
        """
        Retrieve the fractional abundance of neut at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_prot_at_zone():
        """
        Retrieve the fractional abundance of prot at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_bid_at_zone():
        """
        Retrieve the fractional abundance of bid at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def get_mass_fraction_of_bid1_at_zone():
        """
        Retrieve the fractional abundance of bid1 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.OUT,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_h_at_zone():
        """
        Set the fractional abundance of h at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_he3_at_zone():
        """
        Set the fractional abundance of he3 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_he_at_zone():
        """
        Set the fractional abundance of he at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_c12_at_zone():
        """
        Set the fractional abundance of c12 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_c13_at_zone():
        """
        Set the fractional abundance of c13 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_n14_at_zone():
        """
        Set the fractional abundance of n14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_n15_at_zone():
        """
        Set the fractional abundance of n15 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_o16_at_zone():
        """
        Set the fractional abundance of o16 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_o17_at_zone():
        """
        Set the fractional abundance of o17 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_o18_at_zone():
        """
        Set the fractional abundance of o18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_ne20_at_zone():
        """
        Set the fractional abundance of ne20 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_ne22_at_zone():
        """
        Set the fractional abundance of ne22 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_mg24_at_zone():
        """
        Set the fractional abundance of mg24 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_mg25_at_zone():
        """
        Set the fractional abundance of mg25 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_mg26_at_zone():
        """
        Set the fractional abundance of mg26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_c14_at_zone():
        """
        Set the fractional abundance of c14 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_f18_at_zone():
        """
        Set the fractional abundance of f18 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_f19_at_zone():
        """
        Set the fractional abundance of f19 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_ne21_at_zone():
        """
        Set the fractional abundance of ne21 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_na23_at_zone():
        """
        Set the fractional abundance of na23 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_al26_at_zone():
        """
        Set the fractional abundance of al26 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_al27_at_zone():
        """
        Set the fractional abundance of al27 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_si28_at_zone():
        """
        Set the fractional abundance of si28 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_neut_at_zone():
        """
        Set the fractional abundance of neut at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_prot_at_zone():
        """
        Set the fractional abundance of prot at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_bid_at_zone():
        """
        Set the fractional abundance of bid at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


    @legacy_function
    def set_mass_fraction_of_bid1_at_zone():
        """
        Set the fractional abundance of bid1 at the specified
        zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to set the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to set the value of")
        function.addParameter(
            'Xj_i', dtype='float64', direction=function.IN,
            description=(
                "The fractional chemical abundance variable at the specified "
                "zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function



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

        # for set_name in ['fullparticles']:
        for set_name in ['particles', 'fullparticles']:
            handler.define_set(set_name, 'index_of_the_star')
            InternalStellarStructure.define_particle_sets(
                self, handler, set_name=set_name
            )
            handler.set_new(set_name, 'new_particle')

            for parameter in {
                **GENEC_STAR_PARAMETERS,
                **GENEC_STAR_PROPERTIES,
                # **GENEC_STAR_STRUCTURE
            }.items():
                if len(parameter[1]) == 4:
                    names = (
                        parameter[1][3],
                    )
                else:
                    names = (parameter[0],)
                handler.add_getter(
                    set_name,
                    f'get_{parameter[0]}',
                    names=names,
                )
            handler.add_getter(set_name, 'get_radius')
            handler.add_getter(set_name, 'get_mass')
            handler.add_getter(set_name, 'get_age')
            handler.add_getter(set_name, 'get_luminosity')
            handler.add_getter(set_name, 'get_temperature')
            handler.add_getter(set_name, 'get_time_step', names=('time_step',))
            handler.add_getter(
                set_name, 'get_number_of_species', names=('n_species',)
            )

            # handler.add_method(set_name, 'get_radius_profile')
            # handler.add_method(set_name, 'get_temperature_profile')
            # handler.add_method(set_name, 'get_luminosity_profile')
            handler.add_method(set_name, 'get_mass_profile')
            handler.add_method(set_name, 'get_eps_profile')
            handler.add_method(set_name, 'get_epsy_profile')
            handler.add_method(set_name, 'get_eps_c_adv_profile')
            handler.add_method(set_name, 'get_eps_ne_adv_profile')
            handler.add_method(set_name, 'get_eps_o_adv_profile')
            handler.add_method(set_name, 'get_eps_si_adv_profile')
            handler.add_method(set_name, 'get_eps_grav_profile')
            handler.add_method(set_name, 'get_eps_nu_profile')
            handler.add_method(set_name, 'get_cumulative_mass_profile')
            handler.add_getter(
                set_name, 'get_surface_velocity', names=('surface_velocity',)
            )

            handler.add_method(set_name, 'evolve_one_step')
            handler.add_method(set_name, 'evolve_for')
            handler.set_delete(set_name, 'delete_star')

            handler.add_gridded_getter(
                set_name,
                'get_nabla_rad_at_zone', 'get_firstlast_zone',
                names=('nabla_rad_profile',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_nabla_ad_at_zone', 'get_firstlast_zone',
                names=('nabla_ad_profile',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_nabla_mu_at_zone', 'get_firstlast_zone',
                names=('nabla_mu_profile',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_eps_at_zone', 'get_firstlast_zone',
                names=('eps_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_epsy_at_zone', 'get_firstlast_zone',
                names=('epsy_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_c_adv_at_zone', 'get_firstlast_zone',
                names=('eps_c_adv_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_ne_adv_at_zone', 'get_firstlast_zone',
                names=('eps_ne_adv_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_o_adv_at_zone', 'get_firstlast_zone',
                names=('eps_o_adv_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_si_adv_at_zone', 'get_firstlast_zone',
                names=('eps_si_adv_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_grav_at_zone', 'get_firstlast_zone',
                names=('eps_grav_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_eps_nu_at_zone', 'get_firstlast_zone',
                names=('eps_nu_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_radius_at_zone', 'get_firstlast_zone',
                names=('radius_profile',)
            )
            handler.add_gridded_setter(
                set_name,
                'set_radius_at_zone', 'get_firstlast_zone',
                names=('radius_profile',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_temperature_at_zone', 'get_firstlast_zone',
                names=('temperature_profile',)
            )
            handler.add_gridded_setter(
                set_name,
                'set_temperature_at_zone', 'get_firstlast_zone',
                names=('temperature_profile',)
            )

            handler.add_gridded_getter(
                set_name,
                'get_density_at_zone', 'get_firstlast_zone',
                names=('density_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'set_density_at_zone', 'get_firstlast_zone',
                names=('density_profile',)
            )
            handler.add_gridded_getter(
                set_name,
                'get_luminosity_at_zone', 'get_firstlast_zone',
                names=('luminosity_profile',),
            )

            handler.add_gridded_getter(
                set_name,
                'get_pressure_at_zone', 'get_firstlast_zone',
                names=('pressure_profile',)
            )
            # handler.add_gridded_setter(
            #     set_name,
            #     'set_pressure_at_zone', 'get_firstlast_zone',
            #     names=('pressure_profile',)
            # )

            for species in SPECIES_NAMES:
                handler.add_gridded_getter(
                    set_name,
                    f'get_mass_fraction_of_{species}_at_zone',
                    'get_firstlast_zone',
                    names=(f'abundance_{species}',)
                )

                handler.add_gridded_setter(
                    set_name,
                    f'set_mass_fraction_of_{species}_at_zone',
                    'get_firstlast_zone',
                    names=(f'abundance_{species}',)
                )

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
        StellarEvolution.define_state(self, handler)
        # InternalStellarStructure.define_state(self, handler)

        # Only allow setting of starname in EDIT or UPDATE states
        # I.e. must do initialize_code and commit_parameters FIRST!

        # Initialized (initialize_code)
        # handler.add_method

        # -> Edit (commit_parameters)
        # handler.add_method('EDIT', 'set_starname')
        # handler.add_method('EDIT', 'new_particle')

        # -> Run (commit_particles)
        handler.add_transition('EDIT', 'RUN', 'commit_particles')

        for state in ["UPDATE"]:
            for parameter in {
                **GENEC_STAR_PARAMETERS,
                **GENEC_STAR_PROPERTIES,
            }:
                handler.add_method(state, f'get_{parameter[0]}')
            handler.add_method(state, 'get_chemical_abundance_profiles')
            handler.add_method(state, 'get_mass_fraction_of_species_at_zone')
            handler.add_method(state, 'get_mu_at_zone')
            handler.add_method(state, 'get_number_of_species')
            handler.add_method(state, 'get_pressure_at_zone')
            handler.add_method(state, 'get_radius')
            handler.add_method(state, 'get_radius_at_zone')
            handler.add_method(state, 'get_eps_at_zone')
            handler.add_method(state, 'get_epsy_at_zone')
            handler.add_method(state, 'get_eps_c_adv_at_zone')
            handler.add_method(state, 'get_eps_ne_adv_at_zone')
            handler.add_method(state, 'get_eps_o_adv_at_zone')
            handler.add_method(state, 'get_eps_si_adv_at_zone')
            handler.add_method(state, 'get_eps_grav_at_zone')
            handler.add_method(state, 'get_eps_nu_at_zone')
            handler.add_method(state, 'get_surface_velocity')
            handler.add_method(state, 'get_temperature')
            handler.add_method(state, 'get_temperature_at_zone')
            handler.add_method(state, 'get_density_at_zone')
            handler.add_method(state, 'get_luminosity')
            # handler.add_method(state, 'get_luminosity_at_zone')
            handler.add_method(state, 'get_time_step')
            handler.add_method(state, 'get_mass')
            handler.add_method(state, 'get_age')
            handler.add_method(state, 'get_mass_fraction_at_zone')
            for species in SPECIES_NAMES:
                handler.add_method(
                    state, f'get_mass_fraction_of_{species}_at_zone'
                )

        # -> Update
        handler.add_transition('RUN', 'UPDATE', 'finalize_stellar_model')

        handler.add_method('UPDATE', 'set_n_snap')
        #handler.add_method('UPDATE', 'set_ipoly')
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

        #handler.add_method('UPDATE', 'set_starname')
        # handler.add_method('UPDATE', 'new_particle')

    def define_methods(self, handler):
        InternalStellarStructure.define_methods(self, handler)
        StellarEvolution.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (units.MSun, handler.NO_UNIT, handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
        )

        handler.add_method(
            "read_genec_model",
            (handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
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

    def new_particle_from_model(
        self, internal_structure, current_age=0 | units.julianyr, key=None
    ):
        self.new_stellar_model(
            internal_structure[0], internal_structure[1],
            internal_structure[2], internal_structure[3],
            internal_structure[4], internal_structure[5],
            internal_structure[6], internal_structure[7],
            internal_structure[8], internal_structure[9],
            internal_structure[10], internal_structure[11],
            internal_structure[12], internal_structure[13],
            internal_structure[14], internal_structure[15],
            internal_structure[16], internal_structure[17],
            internal_structure[18], internal_structure[19],
            internal_structure[20], internal_structure[21],
            internal_structure[22], internal_structure[23],
            internal_structure[24], internal_structure[25],
            internal_structure[26], internal_structure[27],
            internal_structure[28], internal_structure[29],
            internal_structure[30], internal_structure[31],
            internal_structure[32], internal_structure[33],
            internal_structure[34], internal_structure[35],
            internal_structure[36], internal_structure[37],
            internal_structure[38], internal_structure[39],
            internal_structure[40], internal_structure[41],
            internal_structure[42], internal_structure[43],
            internal_structure[44], internal_structure[45],
            internal_structure[46], internal_structure[47],
            internal_structure[48], internal_structure[49],
            internal_structure[50], internal_structure[51],
            internal_structure[52], internal_structure[53],
            internal_structure[54], internal_structure[55],
            internal_structure[56], internal_structure[57],
            internal_structure[58], internal_structure[59],
            internal_structure[60], internal_structure[61],
            internal_structure[62], internal_structure[63],
            internal_structure[64], internal_structure[65],
            internal_structure[66], internal_structure[67],
            internal_structure[68], internal_structure[69],
            internal_structure[70], internal_structure[71],
            internal_structure[72], internal_structure[73],
            internal_structure[74], internal_structure[75],
            internal_structure[76], internal_structure[77],
            internal_structure[78], internal_structure[79],
            internal_structure[80], internal_structure[81],
            internal_structure[82], internal_structure[83],
            internal_structure[84], internal_structure[85],
            internal_structure[86], internal_structure[87],
            internal_structure[88], internal_structure[89],
            internal_structure[90], internal_structure[91],
            internal_structure[92], internal_structure[93],
            internal_structure[94], internal_structure[95],
            internal_structure[96], internal_structure[97],
            internal_structure[98], internal_structure[99],
            internal_structure[100], internal_structure[101],
            internal_structure[102], internal_structure[103],
            internal_structure[104], internal_structure[105],
            internal_structure[106], internal_structure[107],
            internal_structure[108], internal_structure[109],
            internal_structure[110], internal_structure[111],
            internal_structure[112], internal_structure[113],
            internal_structure[114], internal_structure[115],
            internal_structure[116], internal_structure[117],
            internal_structure[118], internal_structure[119],
            internal_structure[120], internal_structure[121],
            internal_structure[122], internal_structure[123],
            internal_structure[124], internal_structure[125],
            internal_structure[126], internal_structure[127],
            internal_structure[128], internal_structure[129],
            internal_structure[130], internal_structure[131],
            internal_structure[132], internal_structure[133],
            internal_structure[134], internal_structure[135],
            internal_structure[136], internal_structure[137],
            internal_structure[138], internal_structure[139],
            internal_structure[140], internal_structure[141],
            internal_structure[142], internal_structure[143],
            internal_structure[144], internal_structure[145],
            internal_structure[146], internal_structure[147],
            internal_structure[148], internal_structure[149],
            internal_structure[150], internal_structure[151],
            internal_structure[152], internal_structure[153],
            internal_structure[154], internal_structure[155],
            internal_structure[156], internal_structure[157],
            internal_structure[158], internal_structure[159],
            internal_structure[160], internal_structure[161],
            internal_structure[162], internal_structure[163],
            internal_structure[164], internal_structure[165],
            internal_structure[166], internal_structure[167],
            internal_structure[168], internal_structure[169],
            internal_structure[170], internal_structure[171],
            internal_structure[172], internal_structure[173],
            internal_structure[174], internal_structure[175],
            internal_structure[176],
            )
        return
