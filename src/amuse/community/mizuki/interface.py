from amuse.community import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
    LiteratureReferencesMixIn,
)

from amuse.community.interface.gd import (
    GravitationalDynamicsInterface,
    GravitationalDynamics,
    # GravityFieldInterface,
    # GravityFieldCode,
)
from amuse.community.interface.stopping_conditions import (
    StoppingConditionInterface,
    StoppingConditions,
)
from amuse.units import nbody_system


class MizukiInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    GravitationalDynamicsInterface,
    StoppingConditionInterface,
    # GravityFieldInterface,
):
    """
    Mizuki: based on c++ nbody+sph example of FDPS

    References:
        .. [#] Iwasawa et al. (FDPS)
    """

    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self,
            name_of_the_worker="mizuki_worker",
            **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def new_sph_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
        )
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN, default=0.01,
        )
        function.result_type = 'int32'
        return function


    def get_state(self, index_of_the_particle):
        return self.get_state_sph(index_of_the_particle)

    @legacy_function
    def get_state_sph():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle to get the state from. This
            index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.OUT,
            description="The current mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'u', dtype='float64', direction=function.OUT,
            description="The current specific energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The current smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found
        -1 - ERROR
            particle could not be found"""
        return function

    @legacy_function
    def set_state_sph():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="Identifier of the particle"
        )
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The new mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'u', dtype='float64', direction=function.IN,
            description="The new specific energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description="The new smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_mass():
        """
        Retrieve the mass of a particle. Mass is a scalar property of a
        particle, this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=(
                "Index of the particle to get the state from."
                " This index must have been returned by an earlier call to"
                " :meth:`new_particle`"
            )
        )
        function.addParameter(
            'mass', dtype='float64', direction=function.OUT,
            description="The current mass of the particle"
        )
        function.result_type = 'int32'
        # function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'x', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.addParameter(
            'y', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.addParameter(
            'z', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'vx', dtype='float64', direction=function.OUT,
            unit=nbody_system.speed,
        )
        function.addParameter(
            'vy', dtype='float64', direction=function.OUT,
            unit=nbody_system.speed,
        )
        function.addParameter(
            'vz', dtype='float64', direction=function.OUT,
            unit=nbody_system.speed,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'density', dtype='float64', direction=function.OUT,
            description="The density of the particle",
            unit=nbody_system.density,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'pressure', dtype='float64', direction=function.OUT,
            description="The pressure of the particle",
            unit=nbody_system.pressure,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="Smoothing length",
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description="Smoothing length",
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'u', dtype='float64', direction=function.OUT,
            description='Specific energy',
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'radius', dtype='float64', direction=function.OUT,
            description="Smoothing length",
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_radius():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter(
            'radius', dtype='float64', direction=function.IN,
            description="Smoothing length",
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        return function


class Mizuki(GravitationalDynamics,):
    __interface__ = MizukiInterface

    def __init__(
        self,
        convert_nbody=None,
        **options
    ):
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            MizukiInterface(**options),
            convert_nbody,
            **options
        )

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        # GravityFieldCode.define_state(self, handler)
        # self.stopping_conditions.define_state(handler)

        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        handler.add_method('END', 'initialize_code')

        handler.add_method('EDIT', 'new_sph_particle')
        handler.add_method('UPDATE', 'new_sph_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        handler.add_method('RUN', 'get_velocity')
        handler.add_method('RUN', 'get_acceleration')
        handler.add_method('RUN', 'get_internal_energy')
        # handler.add_method('RUN', 'get_dinternal_energy_dt')
        handler.add_method('RUN', 'get_smoothing_length')
        handler.add_method('RUN', 'get_density')
        handler.add_method('RUN', 'get_pressure')
        handler.add_method('RUN', 'get_state_sph')

        self.stopping_conditions.define_state(handler)

    def define_parameters(self, handler):
        # handler.add_method_parameter(
        #     "get_eps2",
        #     "set_eps2",
        #     "epsilon_squared",
        #     "smoothing parameter for gravity calculations",
        #     default_value=0.0 | nbody_system.length * nbody_system.length
        # )

        # handler.add_method_parameter(
        #     "get_dtime",
        #     "set_dtime",
        #     "timestep",
        #     "timestep for system",
        #     default_value=1.0/8 | nbody_system.time
        # )

        # handler.add_boolean_parameter(
        #     "get_use_hydro",
        #     "set_use_hydro",
        #     "use_hydro_flag",
        #     "Hydrodynamics flag. True means: SPH hydro included,"
        #     " False means: gravity only.",
        #     True
        # )
        pass

    def define_particle_sets(self, handler):
        handler.define_super_set(
            'particles',
            ['gas_particles', ],
            index_to_default_set=0
        )

        # handler.define_set('dm_particles', 'index_of_the_particle')
        # handler.set_new('dm_particles', 'new_dm_particle')
        # handler.set_delete('dm_particles', 'delete_particle')
        # handler.add_setter('dm_particles', 'set_state')
        # handler.add_getter('dm_particles', 'get_state')
        # handler.add_setter('dm_particles', 'set_mass')
        # handler.add_getter('dm_particles', 'get_mass', names = ('mass',))
        # handler.add_setter('dm_particles', 'set_position')
        # handler.add_getter('dm_particles', 'get_position')
        # handler.add_setter('dm_particles', 'set_radius')
        # handler.add_getter('dm_particles', 'get_radius')
        # handler.add_setter('dm_particles', 'set_velocity')
        # handler.add_getter('dm_particles', 'get_velocity')

        handler.define_set('gas_particles', 'index_of_the_particle')
        handler.set_new('gas_particles', 'new_sph_particle')
        handler.set_delete('gas_particles', 'delete_particle')
        handler.add_getter('gas_particles', 'get_state_sph')
        handler.add_setter('gas_particles', 'set_state_sph')
        handler.add_getter('gas_particles', 'get_mass')
        handler.add_setter('gas_particles', 'set_mass')
        # handler.add_getter('gas_particles', 'get_radius')
        # handler.add_setter('gas_particles', 'set_radius')
        handler.add_getter('gas_particles', 'get_position')
        handler.add_setter('gas_particles', 'set_position')
        handler.add_getter('gas_particles', 'get_velocity')
        handler.add_setter('gas_particles', 'set_velocity')
        handler.add_getter('gas_particles', 'get_acceleration')
        handler.add_getter('gas_particles', 'get_internal_energy')
        handler.add_setter('gas_particles', 'set_internal_energy')
        # handler.add_getter('gas_particles', 'get_dinternal_energy_dt')
        handler.add_getter('gas_particles', 'get_smoothing_length')
        handler.add_setter('gas_particles', 'set_smoothing_length')
        handler.add_getter('gas_particles', 'get_density', names=('rho',))
        handler.add_getter('gas_particles', 'get_density', names=('density',))
        handler.add_getter('gas_particles', 'get_pressure')

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)

        # handler.add_method(
        #     "new_dm_particle",
        #     (
        #         nbody_system.mass,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.length,
        #     ),
        #     (
        #         handler.INDEX,
        #         handler.ERROR_CODE,
        #     )
        # )

        handler.add_method(
            "new_sph_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_state_sph",
            (
                handler.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_state_sph",
            (
                handler.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_density",
            (
                handler.INDEX,
            ),
            (
                nbody_system.density,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_smoothing_length",
            (
                handler.INDEX,
            ),
            (
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_smoothing_length",
            (
                handler.INDEX,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_pressure",
            (
                handler.INDEX,
            ),
            (
                nbody_system.pressure,
                handler.ERROR_CODE,
            )
        )

        # handler.add_method(
        #     "set_velocity",
        #     (
        #         handler.INDEX,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #     ),
        #     (
        #         handler.ERROR_CODE
        #     )
        # )

        # handler.add_method(
        #     "get_velocity",
        #     (
        #         handler.INDEX,
        #     ),
        #     (
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         handler.ERROR_CODE
        #     )
        # )

        # handler.add_method(
        #     "get_acceleration",
        #     (
        #         handler.INDEX,
        #     ),
        #     (
        #         nbody_system.speed / nbody_system.time,
        #         nbody_system.speed / nbody_system.time,
        #         nbody_system.speed / nbody_system.time,
        #         handler.ERROR_CODE
        #     )
        # )

        handler.add_method(
            "get_internal_energy",
            (
                handler.INDEX,
            ),
            (
                nbody_system.specific_energy,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_internal_energy",
            (
                handler.INDEX,
                nbody_system.specific_energy,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        # handler.add_method(
        #     "get_dinternal_energy_dt",
        #     (
        #         handler.INDEX,
        #     ),
        #     (
        #         nbody_system.specific_energy/nbody_system.time,
        #         handler.ERROR_CODE
        #     )
        # )

        # handler.add_method(
        #     "new_star_particle",
        #     (
        #         nbody_system.mass,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.time,
        #         nbody_system.length,
        #     ),
        #     (
        #         handler.INDEX,
        #         handler.ERROR_CODE,
        #     )
        # )
        # handler.add_method(
        #     "get_state_star",
        #     (
        #         handler.INDEX,
        #     ),
        #     (
        #         nbody_system.mass,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.time,
        #         nbody_system.length,
        #         handler.ERROR_CODE
        #     )
        # )
        # handler.add_method(
        #     "set_state_star",
        #     (
        #         handler.INDEX,
        #         nbody_system.mass,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.length,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.speed,
        #         nbody_system.time,
        #         nbody_system.length,
        #     ),
        #     (
        #         handler.ERROR_CODE,
        #     )
        # )
        # handler.add_method(
        #     "set_star_tform",
        #     (
        #         handler.INDEX,
        #         nbody_system.time,
        #     ),
        #     (
        #         handler.ERROR_CODE,
        #     )
        # )
        # handler.add_method(
        #     "get_star_tform",
        #     (
        #         handler.INDEX,
        #     ),
        #     (
        #         nbody_system.time,
        #         handler.ERROR_CODE
        #     )
        # )

        # handler.add_method(
        #     'get_hydro_state_at_point',
        #     (
        #         nbody_system.length, nbody_system.length,
        #         nbody_system.length, nbody_system.speed, nbody_system.speed,
        #         nbody_system.speed
        #     ),
        #     (
        #         nbody_system.density, nbody_system.momentum_density,
        #         nbody_system.momentum_density, nbody_system.momentum_density,
        #         nbody_system.energy_density, handler.ERROR_CODE
        #     )
        # )

        handler.add_method(
            "get_eps2",
            (
            ),
            (
                nbody_system.length * nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_eps2",
            (
                nbody_system.length * nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        # handler.add_method(
        #     "get_dtime",
        #     (),
        #     (nbody_system.time, handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "set_dtime",
        #     (nbody_system.time, ),
        #     (handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "get_verbosity",
        #     (),
        #     (handler.NO_UNIT, handler.ERROR_CODE,)
        # )
        #
        # handler.add_method(
        #     "set_verbosity",
        #     (handler.NO_UNIT, ),
        #     (handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "get_pboxsize",
        #     (),
        #     (nbody_system.length, handler.ERROR_CODE,)
        # )
        #
        # handler.add_method(
        #     "set_pboxsize",
        #     (nbody_system.length, ),
        #     (handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "get_gamma",
        #     (),
        #     (handler.NO_UNIT, handler.ERROR_CODE,)
        # )
        #
        # handler.add_method(
        #     "set_gamma",
        #     (handler.NO_UNIT, ),
        #     (handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "get_alpha",
        #     (),
        #     (handler.NO_UNIT, handler.ERROR_CODE,)
        # )
        #
        # handler.add_method(
        #     "set_alpha",
        #     (handler.NO_UNIT, ),
        #     (handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "get_beta",
        #     (),
        #     (handler.NO_UNIT, handler.ERROR_CODE,)
        # )

        # handler.add_method(
        #     "set_beta",
        #     (handler.NO_UNIT, ),
        #     (handler.ERROR_CODE,)
        # )
        # handler.add_method(
        #     "get_thermal_energy",
        #     (),
        #     (nbody_system.energy, handler.ERROR_CODE,)
        # )
        #
        # handler.add_method(
        #     "get_total_energy",
        #     (),
        #     (nbody_system.energy, handler.ERROR_CODE,)
        # )

        # self.stopping_conditions.define_methods(handler)
