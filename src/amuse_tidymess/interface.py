import warnings

from amuse.community.interface.gd import (
    GravitationalDynamics,
    GravitationalDynamicsInterface,
    GravityFieldInterface,
    GravityFieldCode
)
from amuse.community.interface.stopping_conditions import (
    StoppingConditions, StoppingConditionInterface
)
from amuse.rfi.core import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
    remote_function
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.units import constants, nbody_system, units as u


particle_attributes = [
    'mass',
    'radius',
    'x', 'y', 'z',
    'vx', 'vy', 'vz',
    'xi', 'kf', 'tau',
    'wx', 'wy', 'wz',
    'a_mb'
]

class TidymessInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    GravitationalDynamicsInterface,
    GravityFieldInterface,
    StoppingConditionInterface,
):
    """
    Tidymess is a N-Body code with tides

    .. [#] Boekholt & Correia (MNRAS 2023, vol. 522, pp. 2885–2900)
    """

    include_headers = ['tidymess_worker.h', 'stopcond.h']

    def __init__(self, **kwargs):
        CodeInterface.__init__(
            self, name_of_the_worker='tidymess_worker', **kwargs
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def new_particle():
        """
        Add a particle to Tidymess.

        Parameters
        ----------
        mass : float
            Particle mass.
        radius : float
            Particle radius.
        x, y, z : float
            Particle position.
        vx, vy, vz : float
            Particle velocity.
        xi : float
            Particle moment of inertia.
        kf : float
            Particle fluid love number.
        tau : float
            Particle fluid relaxation time.
        wx, wy, wz : float
            Particle spin vector.
        a_mb : float
            Particle magnetic breaking coefficient.

        Returns
        -------
        index_of_the_particle : int
            Particle index.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.OUT)
        for name in particle_attributes:
            function.addParameter(name, dtype='d', direction=function.IN)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            particle was created and added to Tidymess
        -1 - ERROR
            particle could not be created
        """
        return function

    @legacy_function
    def get_state():
        """
        Retieve the state of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        mass : float
            Particle mass.
        radius : float
            Particle radius.
        x, y, z : float
            Particle position.
        vx, vy, vz : float
            Particle velocity.
        xi : float
            Particle moment of inertia.
        kf : float
            Particle fluid love number.
        tau : float
            Particle fluid relaxation time.
        wx, wy, wz : float
            Particle spin vector.
        a_mb : float
            Particle magnetic breaking coefficient.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for name in particle_attributes:
            function.addParameter(name, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            state was retrieved from particle
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_state():
        """
        Set the state of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        mass : float
            Particle mass.
        radius : float
            Particle radius.
        x, y, z : float
            Particle position.
        vx, vy, vz : float
            Particle velocity.
        xi : float
            Particle moment of inertia.
        kf : float
            Particle fluid love number.
        tau : float
            Particle fluid relaxation time.
        wx, wy, wz : float
            Particle spin vector.
        a_mb : float
            Particle magnetic breaking coefficient.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for name in particle_attributes:
            function.addParameter(name, dtype='d', direction=function.IN)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @remote_function
    def get_xi(index_of_the_particle='i'):
        """
        Retrieve the moment of inertia of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        xi : float
            Moment of inertia of the particle.
        """
        returns (xi='d')

    @remote_function
    def set_xi(index_of_the_particle='i', xi='d'):
        """
        Set the moment of inertia of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        xi : float
            Moment of inertia of the particle.
        """
        returns ()

    @remote_function
    def get_kf(index_of_the_particle='i'):
        """
        Retrieve the fluid love number of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        kf : float
            Fluid love number of the particle.
        """
        returns (kf='d')

    @remote_function
    def set_kf(index_of_the_particle='i', kf='d'):
        """
        Set the fluid love number of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        kf : float
            Fluid love number of the particle.
        """
        returns ()

    @remote_function
    def get_tau(index_of_the_particle='i'):
        """
        Retrieve the fluid relaxation time of a particle.

        Paramters
        ---------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        tau : float
            Fluid relaxation time of the particle.
        """
        returns (tau='d')

    @remote_function
    def set_tau(index_of_the_particle='i', tau='d'):
        """
        Retrieve the fluid relaxation time of a particle.

        Paramters
        ---------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        tau : float
            Fluid relaxation time of the particle.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_spin(index_of_the_particle='i'):
        """
        Retrieve the spin vector of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        wx, wy, wz : float
            Spin vector of the particle.
        """
        returns (wx='d', wy='d', wz='d')

    @remote_function(can_handle_array=True)
    def set_spin(index_of_the_particle='i', wx='d', wy='d', wz='d'):
        """
        Set the spin vector of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        wx, wy, wz : float
            Spin vector of the particle.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_amb(index_of_the_particle='i'):
        """
        Retrieve the magnetic braking coefficient of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        a_mb : float
            Magnetic braking coefficient of the particle.
        """
        returns (a_mb='d')

    @remote_function(can_handle_array=True)
    def set_amb(index_of_the_particle='i', a_mb='d'):
        """
        Set the magnetic braking coefficient of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        a_mb : float
            Magnetic braking coefficient of the particle.
        """
        returns ()

    @remote_function
    def get_tidal_model():
        """
        Get Tidymess tidal model.

        Returns
        -------
        tidal_model : int
            Tidal model. One of: 0=none, 1=conservative, 2=linear,
            3=creep direct, 4=creep tidymess.
        """
        returns (tidal_model='i')

    @remote_function
    def set_tidal_model(tidal_model='i'):
        """
        Set Tidymess tidal model.

        Parameters
        ----------
        tidal_model : int
            Tidal model. One of: 0=none, 1=conservative, 2=linear,
            3=creep direct, 4=creep tidymess.
        """
        returns ()

    @remote_function
    def get_pn_order():
        """
        Get Tidymess pn order parameter.

        Returns
        -------
        pn_order : int
            Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn
        """
        returns (pn_order='i')

    @remote_function
    def set_pn_order(pn_order='i'):
        """
        Set Tidymess pn order parameter.

        Parameters
        ----------
        pn_order : int
            Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn
        """
        returns ()

    @remote_function
    def get_magnetic_braking():
        """
        Get magnetic braking parameter.

        Returns
        -------
        magnetic_braking : int
            Magnetic braking: 0=off, 1=on
        """
        returns (magnetic_braking='i')

    @remote_function
    def set_magnetic_braking(magnetic_braking='i'):
        """
        Set magnetic braking parameter.

        Parameters
        ----------
        magnetic_braking : int
            Magnetic braking: 0=off, 1=on
        """
        returns ()

    @remote_function
    def get_speed_of_light():
        """
        Retrieve Tidymess speed of light parameter.

        Returns
        -------
        speed_of_light : float
            Speed of light. Only used in conjunction with
            N-body units and `pn_order>0`, otherwise equal to c.
        """
        returns (speed_of_light='d' | nbody_system.speed)

    @remote_function
    def set_speed_of_light(speed_of_light='d' | nbody_system.speed):
        """
        Set Tidymess speed of light parameter.

        Parameters
        ----------
        speed_of_light : float
            Speed of light. Only used in conjunction with
            N-body units and `pn_order>0`, otherwise equal to c.
        """
        returns ()

    @remote_function
    def get_dt_mode():
        """
        Retrieve Tidymess dt mode parameter.

        Returns
        -------
        dt_mode : int
            timestep mode: 0=constant dt, 1=adaptive dt, 2=adaptive weighted dt.
        """
        returns (dt_mode='i')

    @remote_function
    def set_dt_mode(dt_mode='i'):
        """
        Set Tidymess dt mode parameter.

        Parameters
        ----------
        dt_mode : int
            timestep mode: 0=constant dt, 1=adaptive dt, 2=adaptive weighted dt.
        """
        returns ()

    @remote_function
    def get_dt_const():
        """
        Retrieve Tidymess constant dt parameter.

        Only applicable if `dt_mode = 0`.

        Returns
        -------
        dt_const : float
            Constant dt value.
        """
        returns (dt_const='d' | nbody_system.time)

    @remote_function
    def set_dt_const(dt_const='d' | nbody_system.time):
        """
        Set Tidymess constant dt parameter.

        Only applicable if `dt_mode = 0`.

        Parameters
        ----------
        dt_const : float
            Constant dt value.
        """
        returns ()

    @remote_function
    def get_eta():
        """
        Retrieve Tidymess eta (accuracy) parameter.

        Returns
        -------
        eta : float
            Integration accuracy parameter.
        """
        returns (eta='d')

    @remote_function
    def set_eta(eta='d'):
        """
        Set Tidymess eta (accuracy) parameter.

        Parameters
        ----------
        eta : float
            Integration accuracy parameter.
        """
        returns ()

    @remote_function
    def get_n_iter():
        """
        Retrieve Tidymess n_iter parameter. This is the
        number of iterations to improve reversibility.

        Returns
        -------
        n_iter : int
            Number of integration iterations.
        """
        returns (n_iter='i')

    @remote_function
    def set_n_iter(n_iter='i'):
        """
        Set Tidymess n_iter parameter. This is the
        number of iterations to improve reversibility.

        Parameters
        ----------
        n_iter : int
            Number of integration iterations.
        """
        returns ()

    @remote_function
    def get_collision_mode():
        """
        Get Tidymess collision mode parameter.

        Returns
        -------
        collision_mode : int
            Collision mode: 0=off, 1=flag, 2=exception, 3=replace.
        """
        returns (collision_mode='i')

    @remote_function
    def set_collision_mode(collision_mode='i'):
        """
        Set Tidymess collision mode parameter.

        Parameters
        ----------
        collision_mode : int
            Collision mode: 0=off, 1=flag, 2=exception, 3=replace.
        """
        returns ()

    @remote_function
    def get_roche_mode():
        """
        Retrieve Tidymess roche mode parameter.

        Returns
        -------
        roche_mode : int
            Roche mode parameter: 0=off, 1=flag, 2=exception.
        """
        returns (roche_mode='i')

    @remote_function
    def set_roche_mode(roche_mode='i'):
        """
        Set Tidymess roche mode parameter.

        Parameters
        ----------
        roche_mode : int
            Roche mode parameter: 0=off, 1=flag, 2=exception.
        """
        returns ()

    @remote_function
    def get_breakup_mode():
        """
        Retrieve Tidymess breakup mode parameter.

        Returns
        -------
        breakup_mode : int
            Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception
        """
        returns (breakup_mode='i')

    @remote_function
    def set_breakup_mode(breakup_mode='i'):
        """
        Set Tidymess breakup mode parameter.

        Parameters
        ----------
        breakup_mode : int
            Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception
        """
        returns ()

    @remote_function
    def get_initial_shape():
        """
        Get Tidymess initial shape parameter.

        Returns
        -------
        initial_shape:
            Initial of bodies: 0=sphere, 1=equilibrium.
        """
        returns (initial_shape='i')

    @remote_function
    def set_initial_shape(initial_shape='i'):
        """
        Set Tidymess initial shape parameter.

        Parameters
        ----------
        initial_shape:
            Initial of bodies: 0=sphere, 1=equilibrium.
        """
        returns ()

    @remote_function
    def get_num_integration_step():
        """
        Retrieve Tidymess num_integration_step.
        """
        returns (num_integration_step='i')

    @remote_function
    def get_total_energy():
        """
        Retrieve total energy of the system.

        Returns
        -------
        total_energy : float
            System total energy.
        """
        returns (total_energy='d')

    @remote_function(can_handle_array=True)
    def convert_spin_vectors_to_inertial(
        lod='d', obl='d', psi='d'
    ):
        """
        Convert spin vector {length of day, obliquity, spin precession angle}
        to spin vector {wx, wy, wz}.

        Useful for converting initial conditions into the format readable
        by Tidymess.

        Parameters
        ----------
        lod, obl, psi : float
            Spin vector defined as {length of day, obliquity, spin precession angle}.

        Returns
        -------
        wx, wy, wz : float
            Converted spin vector in {wx, wy, wz} form.
        """
        returns (wx='d', wy='d', wz='d')


class Tidymess(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody=None, **options):

        legacy_interface = TidymessInterface(**options)
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

        if convert_nbody is not None:
            self.parameters.speed_of_light = constants.c

    def commit_parameters(self):
        """Check speed of light before committing parameters."""
        self._check_pn_speed_of_light()
        super().commit_parameters()

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)

        handler.add_method('RUN', 'get_spin')
        handler.add_transition('RUN', 'UPDATE', 'set_spin', False)
        handler.add_method('UPDATE', 'set_spin')

        self.stopping_conditions.define_state(handler)


    def define_methods(self, handler):
        """Map TidymessInterface legacy_functions into Tidymess user methods."""
        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            'new_particle',
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.NO_UNIT,
                handler.NO_UNIT,
                nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                nbody_system.time,
            ),
            (handler.INDEX, handler.ERROR_CODE)
        )

        handler.add_method(
            'get_state',
            (handler.INDEX),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.NO_UNIT,
                handler.NO_UNIT,
                nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                nbody_system.time,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            'set_state',
            (
                handler.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.NO_UNIT,
                handler.NO_UNIT,
                nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                1 / nbody_system.time,
                nbody_system.time,
            ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_xi',
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_xi',
            (handler.INDEX, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_kf',
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_kf',
            (handler.INDEX, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_tau',
            (handler.INDEX,),
            (nbody_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_tau',
            (handler.INDEX, nbody_system.time),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_spin',
            (handler.INDEX,),
            (
                1/nbody_system.time,
                1/nbody_system.time,
                1/nbody_system.time,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            'set_spin',
            (
                handler.INDEX,
                1/nbody_system.time,
                1/nbody_system.time,
                1/nbody_system.time,
            ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_amb',
            (handler.INDEX,),
            (nbody_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_amb',
            (handler.INDEX, nbody_system.time),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_num_integration_step',
            (),
            (handler.INDEX, handler.ERROR_CODE)
        )

        handler.add_method(
            'convert_spin_vectors_to_inertial',
            (
                nbody_system.time,
                u.rad,
                u.rad,
            ),
            (
                1/nbody_system.time,
                1/nbody_system.time,
                1/nbody_system.time,
                handler.ERROR_CODE
            )
        )

        self.stopping_conditions.define_methods(handler)


    def define_parameters(self, handler):
        """
        Define setters and getters as Tidymess.parameters.
        This is the standard API for users to access and modify
        Tidymess parameters. Users should use parameters rather
        than acessing the setters and getters defined in
        TidymessInterface directly. Accessors used here must be
        defined in TidymessInterface as a legacy_function, as
        well as the interface.cc.
        User access then follows the syntax:
            Tidymess.parameters.param_name = value
        Users should not use:
            Tidymess.set_param_name(value)
        """
        GravitationalDynamics.define_parameters(self, handler)
        handler.add_method_parameter(
            'get_tidal_model',
            'set_tidal_model',
            'tidal_model',
            '0=none (default), 1=conservative, 2=linear, 3=creep direct, 4=creep tidymess',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_pn_order',
            'set_pn_order',
            'pn_order',
            'Post-Newtonian order: 0=none (default), 1=1pn, 2=1+2pn, 25=1+2+2.5pn',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_magnetic_braking',
            'set_magnetic_braking',
            'magnetic_braking',
            'Magnetic braking. 0=off (default), 1=on',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_speed_of_light',
            'set_speed_of_light',
            'speed_of_light',
            'Speed of light. Should be set depending on unit system being used, or if using post newtonian.',
            default_value=1e100 | nbody_system.speed,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_dt_mode',
            'set_dt_mode',
            'dt_mode',
            '0=constant dt; 1=adaptive dt; 2=adaptive, weighted dt',
            default_value=2,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_dt_const',
            'set_dt_const',
            'dt_const',
            'constant time step',
            default_value=0.015625 | nbody_system.time,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_eta',
            'set_eta',
            'eta',
            'accuracy parameter; timestep multiplication factor, default=0.0625 (only used if dt_mode>0)',
            default_value=0.0625,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_n_iter',
            'set_n_iter',
            'n_iter',
            'Number of iterations to improve reversibility (default=1)',
            default_value=1,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_collision_mode',
            'set_collision_mode',
            'collision_mode',
            '0=off (default), 1=flag, 2=exception, 3=replace',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_roche_mode',
            'set_roche_mode',
            'roche_mode',
            '0=off, 1=flag, 2=exception',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_breakup_mode',
            'set_breakup_mode',
            'breakup_mode',
            'Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        handler.add_method_parameter(
            'get_initial_shape',
            'set_initial_shape',
            'initial_shape',
            'Initial shape, 0=sphere (default), 1=equilibrium',
            default_value=0,
            is_vector=False,
            must_set_before_get=False,
        )

        self.stopping_conditions.define_parameters(handler)


    def define_properties(self, handler):
        """Define read only properties of Tidymess"""
        GravitationalDynamics.define_properties(self, handler)
        handler.add_property('get_total_energy', public_name='total_energy')


    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)

        self.stopping_conditions.define_particle_set(handler)


    def _check_pn_speed_of_light(self):
        """
        Warn if PN is active and the speed of light has not been set.

        Notes
        -----
        - Only relevant when no unit converter is used (`convert_nbody=None`).
        - With a converter, `speed_of_light` is set in `__init__`.
        """
        if self.unit_converter is not None:
            return
        if self.parameters.pn_order == 0:
            return
        # default tidymess sentinel speed of light value
        if self.parameters.speed_of_light == 1e100 | nbody_system.speed:
            warnings.warn(
                "pn_order != 0 but speed_of_light is not set; "
                "set code.parameters.speed_of_light before running.",
                stacklevel=3,
            )
