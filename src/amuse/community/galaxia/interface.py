"""
Interface Definition for the galaxy models
"""
from amuse.units import units
from amuse.community import CodeInterface
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community.interface.gd import GravityFieldCode

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification

# note: angle units can be added (unit.rad) - however breaks current scripts


class BarAndSpiralsInterface(CodeInterface, CommonCodeInterface):
    """
    Galactic model of the Milky Way. The components of the Galaxy are:

    -- Axisymmetric part --
    The axisymmetric component of the Galaxy has a bulge, disk and a dark matter halo.
    the potentials and parameters associated with them are taken from
    Allen & Santillan (1991). In this code, the parameters of
    the axisymmetric part of the Galaxy are defined as:

    CENTRAL BULGE (spherical potential)
    mass_bulge: mass of the bulge
    b_bulge: scale height

    DISK(Miyamoto-Nagai potential)
    mass_disk: mass of the disk
    a_disk: scale lenght
    b_disk: scale height

    DARK MATTER HALO (logaritmic potential)
    mass_halo: mass of the halo
    a_halo: scale lenght

    -- Central bar---
    The central bar of the Milky Way is modelled as a Ferrers bar
    with n=1 ( Romero-Gomez et al. 2007, A&A, 472, 63R; Romero-Gomez et al. 2011, MNRAS, 418, 1176R and references therein)

    The parameters to set the bar potential are:
    bar_contribution= True (default: false)

    The parameters of the bar that can be set in the code are:
    aaxis_bar: semi-major axis of the bar
    axis_ratio_bar: b/a
    bar_phase: initial orientation of the bar
    mass_bar: mass of the bar
    nbt: number of revolutions to obtain a bar. This is used to set an adiabatic growing bar (default: 0, bar already present in the potential)
    omega_bar: Pattern speed of the bar
    tgrowth_bar: if nbt is not 0, growing time of the bar. This value can be only obtained

    -- Spiral arms--
    There are two models for the spiral arms, 2D TWA and 3D CG02

    The parameters to set the spiral arm potential are:
    spiral_contribution= True (default: false)
    spiral_phase: Initial orientation of the arms
    rsigma: scale lenght of the spiral arms.
    m: number of arms

       -- 2D TWA--
    the Tight Winding Approximation. The model and default parameters are set
    according to locus 2 in Antoja et al. 2011, MNRAS, 418, 1423A.

    The parameters of the TWA model that can be set in the code are:
    amplitude: amplitude of the perturbation
    rsp: starting radius of the spiral locus

      -- 3D spiral arms, Cox & Gomez 2002 (CG02)
    Model from Cox & Gomez 2002, ApJS, 142, 261

    The parameters to set this spiral arm model are:
    spiral_model=1 (default:0 = 2D TWA)

    The parameters of the CG02 potential that can be set are:
    fiducial_radius
    scale_height
    spiral_density_amplitude

      --Transient spirals--
    By now the prescription is to have several non overlapping
    spiral events. The transient spirals only can be set when the bar potential is also set.

    Parameters:
    transient_spiral= True (default: false)
    sigma_s: duration of the transient event
    t_sim: total simulation time

    """

    use_modules = ['BarAndSpiralsInterface']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="GalaxyModel_worker", **options)

    @legacy_function
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate the force on
        bodies at those points, multiply with the mass of the bodies
        """
        function = LegacyFunctionSpecification()
        for x in ['eps', 'x', 'y', 'z']:
            function.addParameter(x, dtype='float64', direction=function.IN, unit=units.kpc)
        for x in ['ax', 'ay', 'az']:
            function.addParameter(
               x,
               dtype='float64',
               direction=function.OUT,
               unit=100*units.km**2 * units.s**-2/units.kpc
               )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()
        for x in ['eps', 'x', 'y', 'z']:
            function.addParameter(
               x,
               dtype='float64',
               direction=function.IN,
               unit=units.kpc
               )
        for x in ['phi']:
            function.addParameter(
               x,
               dtype='float64',
               direction=function.OUT,
               unit=100*units.km**2 * units.s**-2
               )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_local_density():
        """
        Retrieve the local stellar density at a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('t', dtype='float64', direction=function.IN, description="time", unit=97781310.5721*units.yr)
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('density', dtype='float64', direction=function.OUT, description="local density", unit=2.32e7*units.MSun/units.kpc**3)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_velcirc():
        """
        Retrieve the circular velocity due to the axisymmetric potetial at a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('vel_circ', dtype='float64', direction=function.OUT, description="circular velocity", unit=10*units.km/units.s)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def get_epifreq():
        """
        Retrieve the epicyclic frequency due to the axisymmetric potetial at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('k', dtype='float64', direction=function.OUT, description="epicyclic freq", unit=10*units.kms)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spiral_density():
        """
        Retrieve the density of the 3D spiral arms at a given point
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('dens', dtype='float64', direction=function.OUT, description="epicyclic freq", unit=2.32e7*units.MSun/units.kpc**3)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tidal_tensor():
        """
        Retrieve the second derivatives of the total force at a given point
        and at a given time.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('t', dtype='float64', direction=function.IN, description="time", unit=97781310.5721*units.yr)
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('Fxx', dtype='float64', direction=function.OUT, description="fxx", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fyx', dtype='float64', direction=function.OUT, description="fyx", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fzx', dtype='float64', direction=function.OUT, description="fzx", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fxy', dtype='float64', direction=function.OUT, description="fxy", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fyy', dtype='float64', direction=function.OUT, description="fyy", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fzy', dtype='float64', direction=function.OUT, description="fzy", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fxz', dtype='float64', direction=function.OUT, description="fxz", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fyz', dtype='float64', direction=function.OUT, description="fyz", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('Fzz', dtype='float64', direction=function.OUT, description="fzz", unit=100*units.kms**2/units.kpc**2)

        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eigen_values():
        """
        Retrieve the eigen values of the tidal tensor.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('t', dtype='float64', direction=function.IN, description="time", unit=97781310.5721*units.yr)
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('lambda1', dtype='float64', direction=function.OUT, description="eigen values", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('lambda2', dtype='float64', direction=function.OUT, description="eigen values", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('lambda3', dtype='float64', direction=function.OUT, description="eigen values", unit=100*units.kms**2/units.kpc**2)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'

        return function

    @legacy_function
    def get_tidal_radius():
        """
        Retrieve the tidal radius of a star cluster
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('t', dtype='float64', direction=function.IN, description="time", unit=97781310.5721*units.yr)
        function.addParameter('x', dtype='float64', direction=function.IN, description="x position", unit=units.kpc)
        function.addParameter('y', dtype='float64', direction=function.IN, description="y position", unit=units.kpc)
        function.addParameter('z', dtype='float64', direction=function.IN, description="z position", unit=units.kpc)
        function.addParameter('mc', dtype='float64', direction=function.IN, description="initial cluster mass", unit=2.32e7*units.MSun)
        function.addParameter('rt', dtype='float64', direction=function.OUT, description="tidal radius", unit=units.kpc)

        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    # The following functions set the constants from Amuse

    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_omega_sys():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_system', dtype='float64', direction=function.OUT, unit=10*units.km/(units.s*units.kpc))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_initial_phase():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_phase', dtype='float64', direction=function.OUT)  # unit=units.rad
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('xflag', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_flag():
        function = LegacyFunctionSpecification()
        function.addParameter('xflag', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    # BAR
    @legacy_function
    def set_bar_phase():
        function = LegacyFunctionSpecification()
        function.addParameter('bar_phase', dtype='float64', direction=function.IN)  # unit=units.rad
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bar_phase():
        function = LegacyFunctionSpecification()
        function.addParameter('bar_phase', dtype='float64', direction=function.OUT)  # unit=units.rad
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_bar', dtype='float64', direction=function.IN, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_bar', dtype='float64', direction=function.OUT, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_aaxis_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('aaxis_bar', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_aaxis_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('aaxis_bar', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_caxis_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('caxis_bar', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_caxis_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('caxis_bar', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_axis_ratio_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('axis_ratio_bar', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_axis_ratio_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('axis_ratio_bar', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_omega_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_bar', dtype='float64', direction=function.IN, unit=10.*units.km/(units.s*units.kpc))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_omega_bar():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_bar', dtype='float64', direction=function.OUT, unit=10.*units.km/(units.s*units.kpc))
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nbt():
        function = LegacyFunctionSpecification()
        function.addParameter('nbt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nbt():
        function = LegacyFunctionSpecification()
        function.addParameter('nbt', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tgrowth():
        function = LegacyFunctionSpecification()
        function.addParameter('tgrowth_bar', dtype='float64', direction=function.OUT, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    # SPIRAL
    @legacy_function
    def set_spiral_phase():
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_phase', dtype='float64', direction=function.IN)  # unit=units.rad
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spiral_phase():
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_phase', dtype='float64', direction=function.OUT)  # unit=units.rad
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_m():
        function = LegacyFunctionSpecification()
        function.addParameter('m', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_m():
        function = LegacyFunctionSpecification()
        function.addParameter('m', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tan_pitch_angle():
        function = LegacyFunctionSpecification()
        function.addParameter('tan_pitch_angle', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tan_pitch_angle():
        function = LegacyFunctionSpecification()
        function.addParameter('tan_pitch_angle', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_r_sigma():
        function = LegacyFunctionSpecification()
        function.addParameter('r_sigma', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_r_sigma():
        function = LegacyFunctionSpecification()
        function.addParameter('r_sigma', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_omega_spiral():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_spiral', dtype='float64', direction=function.IN, unit=10*units.kms/units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_omega_spiral():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_spiral', dtype='float64', direction=function.OUT, unit=10*units.kms/units.kpc)
        function.result_type = 'int32'
        return function

    # TWA

    @legacy_function
    def set_N():
        function = LegacyFunctionSpecification()
        function.addParameter('N', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_N():
        function = LegacyFunctionSpecification()
        function.addParameter('N', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rsp():
        function = LegacyFunctionSpecification()
        function.addParameter('rsp', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_rsp():
        function = LegacyFunctionSpecification()
        function.addParameter('rsp', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amplitude():
        function = LegacyFunctionSpecification()
        function.addParameter('amplitude', dtype='float64', direction=function.IN, unit=100*units.kms**2/units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_amplitude():
        function = LegacyFunctionSpecification()
        function.addParameter('amplitude', dtype='float64', direction=function.OUT, unit=100*units.kms**2/units.kpc)
        function.result_type = 'int32'
        return function

    # CG02 spiral arms
    @legacy_function
    def set_spiral_density_amplitude():
        """
        CG02 model spiral arms density amplitude, measured at fiducial radius
          rho_0 in eqn(1--3) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_density_amplitude', dtype='float64',
                              direction=function.IN, unit=2.32e7*units.MSun/units.kpc**3,
                              description='CG02 model spiral arms density amplitude, measured at fiducial radius')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spiral_density_amplitude():
        """
        CG02 model spiral arms density amplitude, measured at fiducial radius
          rho_0 in eqn(1--3) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_density_amplitude', dtype='float64',
                              direction=function.OUT, unit=2.32e7*units.MSun/units.kpc**3,
                              description='CG02 model spiral arms density amplitude, measured at fiducial radius')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_fiducial_radius():
        """
        CG02 model fiducial radius where the density amplitude and initial phase are measured
          R_0 in eqn(1--3) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('fiducial_radius', dtype='float64', direction=function.IN, unit=units.kpc,
                              description='CG02 model fiducial radius where the density amplitude and initial phase are measured')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_fiducial_radius():
        """
        CG02 model fiducial radius where the density amplitude and initial phase are measured
          R_0 in eqn(1--3) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('fiducial_radius', dtype='float64', direction=function.OUT, unit=units.kpc,
                              description='CG02 model fiducial radius where the density amplitude and initial phase are measured')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_scale_height():
        """
        CG02 model spiral arms scale height
          H in eqn(1) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('scale_height', dtype='float64',
                              direction=function.IN, unit=units.kpc,
                              description='CG02 model spiral arms scale height')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_scale_height():
        """
        CG02 model spiral arms scale height
          H in eqn(1) of CG02
        """
        function = LegacyFunctionSpecification()
        function.addParameter('scale_height', dtype='float64',
                              direction=function.OUT, unit=units.kpc,
                              description='CG02 model spiral arms scale height')
        function.result_type = 'int32'
        return function

    # transcient spirals

    @legacy_function
    def set_t_sim():
        function = LegacyFunctionSpecification()
        function.addParameter('t_sim', dtype='float64', direction=function.IN, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_t_sim():
        function = LegacyFunctionSpecification()
        function.addParameter('t_sim', dtype='float64', direction=function.OUT, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_sigma_s():
        function = LegacyFunctionSpecification()
        function.addParameter('sigma_s', dtype='float64', direction=function.IN, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_sigma_s():
        function = LegacyFunctionSpecification()
        function.addParameter('sigma_s', dtype='float64', direction=function.OUT, unit=97781310.5721*units.yr)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_omega_spiral2():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_spiral2', dtype='float64', direction=function.IN, unit=10.*units.km/(units.s*units.kpc))
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_omega_spiral2():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_spiral2', dtype='float64', direction=function.OUT, unit=10.*units.km/(units.s*units.kpc))
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_amplitude2():
        function = LegacyFunctionSpecification()
        function.addParameter('amplitude2', dtype='float64', direction=function.IN, unit=100*units.kms**2/units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_amplitude2():
        function = LegacyFunctionSpecification()
        function.addParameter('amplitude2', dtype='float64', direction=function.OUT, unit=100*units.kms**2/units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_tan_pitch_angle2():
        function = LegacyFunctionSpecification()
        function.addParameter('tan_pitch_angle2', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_tan_pitch_angle2():
        function = LegacyFunctionSpecification()
        function.addParameter('tan_pitch_angle2', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_m2():
        function = LegacyFunctionSpecification()
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_m2():
        function = LegacyFunctionSpecification()
        function.addParameter('m2', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_phi21():
        function = LegacyFunctionSpecification()
        function.addParameter('phi21_spiral', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_phi21():
        function = LegacyFunctionSpecification()
        function.addParameter('phi21_spiral', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    # _____________________ AXI_____________________________________________
    @legacy_function
    def set_mass_bulge():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_bulge', dtype='float64', direction=function.IN, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass_bulge():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_bulge', dtype='float64', direction=function.OUT, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_b_bulge():
        function = LegacyFunctionSpecification()
        function.addParameter('b_bulge', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_b_bulge():
        function = LegacyFunctionSpecification()
        function.addParameter('b_bulge', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_disk', dtype='float64', direction=function.IN, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_disk', dtype='float64', direction=function.OUT, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_a_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('a_disk', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_a_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('a_disk', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_b_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('b_disk', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_b_disk():
        function = LegacyFunctionSpecification()
        function.addParameter('b_disk', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass_halo():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_halo', dtype='float64', direction=function.IN, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass_halo():
        function = LegacyFunctionSpecification()
        function.addParameter('mass_halo', dtype='float64', direction=function.OUT, unit=2.32e7*units.MSun)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_a_halo():
        function = LegacyFunctionSpecification()
        function.addParameter('a_halo', dtype='float64', direction=function.IN, unit=units.kpc)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_a_halo():
        function = LegacyFunctionSpecification()
        function.addParameter('a_halo', dtype='float64', direction=function.OUT, unit=units.kpc)
        function.result_type = 'int32'
        return function

    # The following function sets the force of the spiral TWA or bar

    @legacy_function
    def set_spiral_contribution():
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_contribution', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spiral_contribution():
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_contribution', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_spiral_model():
        """
        set model of spiral arms
        0 -- TWA (2D)
        1 -- Cox and Gomez (3D)
        2 -- Lepine (2D)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_model', dtype='int32', direction=function.IN,
                              description='model of spiral arms (default: 0; TWA 2D arms)')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_spiral_model():
        """
        set model of spiral arms
        0 -- TWA (2D)
        1 -- Cox and Gomez (3D)
        2 -- Lepine (2D)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('spiral_model', dtype='int32', direction=function.OUT,
                              description='model of spiral arms (default: 0; TWA 2D arms)')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bar_contribution():
        function = LegacyFunctionSpecification()
        function.addParameter('bar_contribution', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_bar_contribution():
        function = LegacyFunctionSpecification()
        function.addParameter('bar_contribution', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_transient_spiral():
        function = LegacyFunctionSpecification()
        function.addParameter('transient_spiral', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_transient_spiral():
        function = LegacyFunctionSpecification()
        function.addParameter('transient_spiral', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    def before_set_parameter(self):
        pass

    def before_get_parameter(self):
        pass


class BarAndSpiralsDoc(object):
    def __get__(self, instance, owner):
        return instance.interface_doc+"\n\n"+instance.parameters.__doc__


class BarAndSpirals3D(CommonCode, GravityFieldCode):
    NBODY = object()

    __doc__ = BarAndSpiralsDoc()

    def __init__(self, unit_converter=None,  **options):
        self.unit_converter = unit_converter
        legacy_interface = BarAndSpiralsInterface(**options)
        self.interface_doc = legacy_interface.__doc__
        CommonCode.__init__(self, legacy_interface, **options)

    def define_parameters(self, handler):
        handler.add_method_parameter(
           "get_time",
           "set_time",
           "time",
           "Evolution time in the model",
           default_value=0 | 97781310.5721*units.yr
           )

        handler.add_method_parameter(
           "get_omega_sys",
           None,
           "omega_system",
           "pattern speed of the system",
           default_value=0 | 10*units.km/(units.s*units.kpc)
           )

        handler.add_method_parameter(
           "get_initial_phase",
           None,
           "initial_phase",
           "phase of the system. To convert between inertial and rotating frames",
           default_value=0  # | units.rad
           )
        handler.add_method_parameter(
           "get_flag",
           "set_flag",
           "xflag",
           "flag for the ref system to compute tidal tensor. 1-> corrotating sys, 2-> inertial sys",
           default_value=0
           )

        # BAR
        handler.add_method_parameter(
           "get_bar_phase",
           "set_bar_phase",
           "bar_phase",
           "Initial phase of the bar",
           default_value=0  # | units.rad
           )

        handler.add_method_parameter(
           "get_mass_bar",
           "set_mass_bar",
           "mass_bar",
           "The mass of the bar in the model",
           default_value=431 | 2.32e7*units.MSun
           )

        handler.add_method_parameter(
           "get_aaxis_bar",
           "set_aaxis_bar",
           "aaxis_bar",
           "The semimajor axis of the bar in the model",
           default_value=3.13 | units.kpc
           )

        handler.add_method_parameter(
           "get_axis_ratio_bar",
           "set_axis_ratio_bar",
           "axis_ratio_bar",
           "b/a of the bar",
           default_value=0.32
           )

        handler.add_method_parameter(
           "get_caxis_bar",
           "set_caxis_bar",
           "caxis_bar",
           "The vertical axis of the bar in the model",
           default_value=0 | units.kpc
           )

        handler.add_method_parameter(
           "get_omega_bar",
           "set_omega_bar",
           "omega_bar",
           "The pattern speed of the bar in the model",
           default_value=5 | 10.*units.km/(units.s*units.kpc)
           )

        handler.add_method_parameter(
           "get_nbt",
           "set_nbt",
           "nbt",
           "The number of rotations of the bar in the model. This is to set tgrow",
           default_value=0
           )

        handler.add_method_parameter(
           "get_tgrowth",
           None,
           "tgrowth_bar",
           "Growing time of the bar",
           default_value=0 | 97781310.5721*units.yr
           )

        # SPIRAL
        handler.add_method_parameter(
           "get_spiral_phase",
           "set_spiral_phase",
           "spiral_phase",
           "Initial phase of the spiral arms",
           default_value=0  # | units.rad
           )

        handler.add_method_parameter(
           "get_N",
           "set_N",
           "N",
           "How sharply spiral change to bar in the central region",
           default_value=100
           )

        handler.add_method_parameter(
           "get_tan_pitch_angle",
           "set_tan_pitch_angle",
           "tan_pitch_angle",
           " tangent of the pitch angle",
           default_value=0.277
           )

        handler.add_method_parameter(
           "get_rsp",
           "set_rsp",
           "rsp",
           "radius when the spiral starts",
           default_value=1.5 | units.kpc
           )

        handler.add_method_parameter(
           "get_amplitude",
           "set_amplitude",
           "amplitude",
           "Amplitude of the perturbation",
           default_value=8.5 | 100*units.kms**2/units.kpc
           )

        handler.add_method_parameter(
           "get_r_sigma",
           "set_r_sigma",
           "r_sigma",
           "scalelength rsigma of the spiral arms",
           default_value=2.5 | units.kpc
           )

        handler.add_method_parameter(
           "get_omega_spiral",
           "set_omega_spiral",
           "omega_spiral",
           "Pattern speed of the spiral",
           default_value=2. | 10.*units.kms/units.kpc
           )

        handler.add_method_parameter(
           "get_m",
           "set_m",
           "m",
           "Number of spirals",
           default_value=2
           )

        # CG02 3D spiral model
        #  default values set as in onriginal CG02 paper

        handler.add_method_parameter(
           "get_spiral_density_amplitude",
           "set_spiral_density_amplitude",
           "spiral_density_amplitude",
           "CG02 model spiral arms density amplitude",
           default_value=1.35633 | 2.32e7*units.MSun/units.kpc**3
           )

        handler.add_method_parameter(
           "get_fiducial_radius",
           "set_fiducial_radius",
           "fiducial_radius",
           "CG02 spiral arms fiducial radius",
           default_value=8.0 | units.kpc
           )

        handler.add_method_parameter(
           "get_scale_height",
           "set_scale_height",
           "scale_height",
           "CG02 spiral arms scale height",
           default_value=0.18 | units.kpc
           )

        # transcient structure

        handler.add_method_parameter(
            "get_sigma_s",
            "set_sigma_s",
            "sigma_s",
            "Duration of the sp force",
            default_value=1.02269032206 | 97781310.5721*units.yr
            )

        handler.add_method_parameter(
            "get_t_sim",
            "set_t_sim",
            "t_sim",
            "simulation time",
            default_value=47.0437548146 | 97781310.5721*units.yr
            )

        # LEPINE MODEL
        handler.add_method_parameter(
           "get_spiral_model",
           "set_spiral_model",
           "spiral_model",
           "0 is TWA, 1 is C&G, 2 is Lepine ",
           default_value=0
           )

        handler.add_method_parameter(
           "get_omega_spiral2",
           "set_omega_spiral2",
           "omega_spiral2",
           "Pattern speed of the second SA ",
           default_value=2. | 10.*units.kms/units.kpc
           )

        handler.add_method_parameter(
           "get_amplitude2",
           "set_amplitude2",
           "amplitude2",
           "Amplitude of the second SA pattern",
           default_value=6.8 | 100*units.kms**2/units.kpc
           )

        handler.add_method_parameter(
           "get_tan_pitch_angle2",
           "set_tan_pitch_angle2",
           "tan_pitch_angle2",
           " tangent of the pitch angle second SA pattern",
           default_value=-0.1227845
           )

        handler.add_method_parameter(
           "get_m2",
           "set_m2",
           "m2",
           "Number of spirals of the second SA pattern",
           default_value=2
           )

        handler.add_method_parameter(
          "get_phi21",
          "set_phi21",
          "phi21_spiral",
          "Initial phase of the second pattern with respect to the primarly one",
          default_value=-3.49065850399
            )

        # AXI
        handler.add_method_parameter(
           "get_mass_bulge",
           "set_mass_bulge",
           "mass_bulge",
           "mass of the total central component in the model",
           default_value=606 | 2.32e7*units.MSun
           )

        handler.add_method_parameter(
           "get_b_bulge",
           "set_b_bulge",
           "b_bulge",
           "b constant of the bulge's potential",
           default_value=0.3873 | units.kpc
           )

        handler.add_method_parameter(
           "get_mass_disk",
           "set_mass_disk",
           "mass_disk",
           "The mass of the disk in the model",
           default_value=3690 | 2.32e7*units.MSun
           )

        handler.add_method_parameter(
           "get_a_disk",
           "set_a_disk",
           "a_disk",
           "The constant a in the potential of the disk",
           default_value=5.3178 | units.kpc
           )

        handler.add_method_parameter(
           "get_b_disk",
           "set_b_disk",
           "b_disk",
           "The constant b in the potential of the disk",
           default_value=0.25 | units.kpc
           )

        handler.add_method_parameter(
           "get_mass_halo",
           "set_mass_halo",
           "mass_halo",
           "The mass of the dark matter halo in the model",
           default_value=4615 | 2.32e7*units.MSun
           )

        handler.add_method_parameter(
           "get_a_halo",
           "set_a_halo",
           "a_halo",
           "The constant a in the potential of the dark matter halo",
           default_value=12.0 | units.kpc
           )

        handler.add_boolean_parameter(
           "get_spiral_contribution",
           "set_spiral_contribution",
           "spiral_contribution",
           "Flag whether to include a spiral in the model",
           False
           )
        handler.add_boolean_parameter(
           "get_bar_contribution",
           "set_bar_contribution",
           "bar_contribution",
           "Flag whether to include a bar in the model",
           False
           )

        handler.add_boolean_parameter(
           "get_transient_spiral",
           "set_transient_spiral",
           "transient_spiral",
           "Flag whether to include transient spirals in the model",
           False
           )

    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED', 'RUN', 'commit_parameters')
        handler.add_transition('RUN', 'CHANGE_PARAMETERS_RUN', 'before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN', 'RUN', 'recommit_parameters')
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        GravityFieldCode.define_state(self, handler)
        handler.add_method('RUN', 'get_local_density')
        handler.add_method('RUN', 'get_velcirc')
        handler.add_method('RUN', 'get_epifreq')
        handler.add_method('RUN', 'get_spiral_density')
        handler.add_method('RUN', 'get_tidal_tensor')
        handler.add_method('RUN', 'get_eigen_values')
        handler.add_method('RUN', 'get_tidal_radius')

    def before_set_parameter(self):
        pass

    def before_get_parameter(self):
        pass

    def evolve_model(self, t_end):
        self.parameters.time = t_end

    @property
    def model_time(self):
        return self.parameters.time

    def get_density_at_point(self, x, y, z):
        return self.get_local_density(self.parameters.time, x, y, z)


Barandspirals3d = BarAndSpirals3D
Galaxia = BarAndSpirals3D
