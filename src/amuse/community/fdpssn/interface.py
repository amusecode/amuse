"""
Interface for the SPH+Nbody FDPS code
"""
from amuse.rfi.core import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.community.interface.gd import (
    GravitationalDynamicsInterface,
    GravitationalDynamics,
)
from amuse.units import nbody_system
# from amuse.support.interface import InCodeComponentImplementation


class FDPSSNInterface(
        CodeInterface,
        LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        # StoppingConditionsInterface,
):
    """
    Interface for the FDPS Nbody+SPH

    .. [#] Iwasawa et al.

    """

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self,
            name_of_the_worker="fdpssn_worker",
            **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    # Particle (attribute) setters/getters

    @legacy_function
    def new_dm_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
        )
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    def new_particle(self, mass, x, y, z, vx, vy, vz):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz)

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
            'h_smooth', dtype='float64', direction=function.IN, default=0.)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state_dm():
        """
        Retrieve the current state of a DM particle. The mass, position and
        velocity are returned.
        """
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
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    def get_state(self, index_of_the_particle):
        return self.get_state_dm(index_of_the_particle)

    @legacy_function
    def get_state_sph():
        """
        Retrieve the current state of an SPH particle. The mass, position,
        velocity, internal energy and smoothing length are returned.
        """
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
            description="The current internal energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The current smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_state_sph():
        """
        Update the current state of an SPH particle. The mass, position,
        velocity, internal energy and smoothing length are updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle for which the state is to be
            updated. This index must have been returned by an earlier call to
            :meth:`new_particle`""")
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
            description="The new internal energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description="The new smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    def set_state(self, index_of_the_particle, mass, x, y, z, vx, vy, vz):
        return self.set_state_dm(
            index_of_the_particle, mass, x, y, z, vx, vy, vz)

    @legacy_function
    def set_state_dm():
        """
        Update the current state of a DM particle. The mass, position and
        velocity are updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle for which the state is to be
            updated. This index must have been returned by an earlier call to
            :meth:`new_particle`""")
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
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'u', dtype='float64', direction=function.IN,
            description=''
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description=''
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'pressure', dtype='float64', direction=function.OUT,
            description=''
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'density', dtype='float64', direction=function.OUT,
            description="The current density of the particle"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'u', dtype='float64', direction=function.OUT,
            description="The current internal energy of the particle"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The current smoothing length of the particle"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    # Parameter setters/getters

    @legacy_function
    def get_epsilon_squared():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'epsilon_squared', dtype='float64', direction=function.OUT,
            description='',
            unit=nbody_system.length**2,
        )
        function.result_type = 'int32'
        function.result_doc = ""
        return function

    @legacy_function
    def set_epsilon_squared():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'epsilon_squared', dtype='float64', direction=function.IN,
            description='',
            unit=nbody_system.length**2,
        )
        function.result_type = 'int32'
        function.result_doc = ""
        return function


class FDPSSN(
        # InCodeComponentImplementation,
        GravitationalDynamics,
):
    """
    Low-level interface to FDPS-SPHNbody
    """

    def __init__(self, convert_nbody=None, **options):
        GravitationalDynamics.__init__(
            self,
            FDPSSNInterface(**options),
            convert_nbody,
            **options
        )

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        # GravityFieldCode.define_state(self, object)

        object.add_method('EDIT', 'new_dm_particle')
        object.add_method('UPDATE', 'new_dm_particle')
        object.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        object.add_method('EDIT', 'new_sph_particle')
        object.add_method('UPDATE', 'new_sph_particle')
        object.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        object.add_method('RUN', 'get_internal_energy')
        object.add_method('RUN', 'get_smoothing_length')
        object.add_method('RUN', 'get_density')
        object.add_method('RUN', 'get_pressure')
        object.add_method('RUN', 'get_state_sph')

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_epsilon_squared",
            "set_epsilon_squared",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value=0.01 | nbody_system.length**2,
        )

    def define_particle_sets(self, object):
        object.define_super_set(
            'particles',
            ['dm_particles', 'gas_particles'],
            index_to_default_set=0,
        )

        object.define_set('dm_particles', 'index_of_the_particle')
        object.set_new('dm_particles', 'new_dm_particle')
        object.set_delete('dm_particles', 'delete_particle')
        object.add_getter('dm_particles', 'get_state_dm')
        object.add_setter('dm_particles', 'set_state_dm')
        object.add_getter('dm_particles', 'get_mass')
        object.add_setter('dm_particles', 'set_mass')
        object.add_getter('dm_particles', 'get_position')
        object.add_setter('dm_particles', 'set_position')
        object.add_getter('dm_particles', 'get_velocity')
        object.add_setter('dm_particles', 'set_velocity')

        object.define_set('gas_particles', 'index_of_the_particle')
        object.set_new('gas_particles', 'new_sph_particle')
        object.set_delete('gas_particles', 'delete_particle')
        object.add_getter('gas_particles', 'get_state_sph')
        object.add_setter('gas_particles', 'set_state_sph')
        object.add_getter('gas_particles', 'get_mass')
        object.add_setter('gas_particles', 'set_mass')
        object.add_getter('gas_particles', 'get_position')
        object.add_setter('gas_particles', 'set_position')
        object.add_getter('gas_particles', 'get_velocity')
        object.add_setter('gas_particles', 'set_velocity')
        object.add_getter('gas_particles', 'get_internal_energy')
        object.add_setter('gas_particles', 'set_internal_energy')
        object.add_getter('gas_particles', 'get_smoothing_length')
        object.add_setter('gas_particles', 'set_smoothing_length')
        object.add_getter('gas_particles', 'get_density', names=('rho',))
        object.add_getter('gas_particles', 'get_density', names=('density',))
        object.add_getter('gas_particles', 'get_pressure')

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            "new_dm_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )

        object.add_method(
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
                object.INDEX,
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "get_state_sph",
            (
                object.INDEX,
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
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "set_state_sph",
            (
                object.INDEX,
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
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "set_internal_energy",
            (
                object.INDEX,
                nbody_system.specific_energy,
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_internal_energy",
            (
                object.INDEX,
            ),
            (
                nbody_system.specific_energy,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_density",
            (
                object.INDEX,
            ),
            (
                nbody_system.density,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_pressure",
            (
                object.INDEX,
            ),
            (
                nbody_system.pressure,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "set_smoothing_length",
            (
                object.INDEX,
                nbody_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_smoothing_length",
            (
                object.INDEX,
            ),
            (
                nbody_system.length,
                object.ERROR_CODE
            )
        )
