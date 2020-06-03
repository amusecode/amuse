from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

class EticsInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn):
    """
        .. [#] Meiron, Y., Li, B., Holley-Bockelmann, K., & Spurzem, R. 2014, ApJ, 792, 98:
        .. [#] ... "Expansion techniques for collisionless stellar dynamical simulations"
        .. [#] [2014ApJ...792...98M]
    """

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker='etics_worker', **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN,
            description = 'The mass of the particle')
        function.addParameter('x', dtype='float64', direction=function.IN,
            description = 'The initial position vector of the particle')
        function.addParameter('y', dtype='float64', direction=function.IN,
            description = 'The initial position vector of the particle')
        function.addParameter('z', dtype='float64', direction=function.IN,
            description = 'The initial position vector of the particle')
        function.addParameter('vx', dtype='float64', direction=function.IN,
            description = 'The initial velocity vector of the particle')
        function.addParameter('vy', dtype='float64', direction=function.IN,
            description = 'The initial velocity vector of the particle')
        function.addParameter('vz', dtype='float64', direction=function.IN,
            description = 'The initial velocity vector of the particle')
        function.addParameter('radius', dtype='float64', direction=function.IN,
            description = 'The radius of the particle', default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN,
                 description = 'The mass of the particle')
        function.addParameter('radius', dtype='float64', direction=function.IN,
                 description = 'The radius of the particle')
        function.addParameter('x', dtype='float64', direction=function.IN,
                 description = 'The initial position vector of the particle')
        function.addParameter('y', dtype='float64', direction=function.IN,
                 description = 'The initial position vector of the particle')
        function.addParameter('z', dtype='float64', direction=function.IN,
                 description = 'The initial position vector of the particle')
        function.addParameter('vx', dtype='float64', direction=function.IN,
                 description = 'The initial velocity vector of the particle')
        function.addParameter('vy', dtype='float64', direction=function.IN,
                 description = 'The initial velocity vector of the particle')
        function.addParameter('vz', dtype='float64', direction=function.IN,
                 description = 'The initial velocity vector of the particle')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT,
                 description = 'The mass of the particle')
        function.addParameter('radius', dtype='float64', direction=function.OUT,
                 description = 'The radius of the particle')
        function.addParameter('x', dtype='float64', direction=function.OUT,
                 description = 'The initial position vector of the particle')
        function.addParameter('y', dtype='float64', direction=function.OUT,
                 description = 'The initial position vector of the particle')
        function.addParameter('z', dtype='float64', direction=function.OUT,
                 description = 'The initial position vector of the particle')
        function.addParameter('vx', dtype='float64', direction=function.OUT,
                 description = 'The initial velocity vector of the particle')
        function.addParameter('vy', dtype='float64', direction=function.OUT,
                 description = 'The initial velocity vector of the particle')
        function.addParameter('vz', dtype='float64', direction=function.OUT,
                 description = 'The initial velocity vector of the particle')
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('dt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT, description = 'number of particles')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN, description = 'time step')
        function.result_type = 'int32'
        return function

    @legacy_function
    def update_force_potential_arrays():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('dt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

class Etics(GravitationalDynamics):
    def __init__(self, convert_nbody = None, **keyword_arguments):
        legacy_interface = EticsInterface(**keyword_arguments)
        GravitationalDynamics.__init__(self, legacy_interface, convert_nbody, **keyword_arguments)

    def define_parameters(self, handler):
        handler.add_method_parameter(
            'get_begin_time',
            'set_begin_time',
            'begin_time',
            'model time to start the simulation at',
            default_value = 0.0 | nbody_system.time
        )
        handler.add_method_parameter(
            'get_time_step',
            'set_time_step',
            'time_step',
            'constant timestep for iteration',
            default_value = 0.001953125 | nbody_system.time
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        # Define some shortcuts for better readability.
        M = nbody_system.mass
        L = nbody_system.length
        V = nbody_system.speed
        T = nbody_system.time
        handler.add_method('new_particle', (M,L,L,L,V,V,V,L), (handler.INDEX, handler.ERROR_CODE))
        handler.add_method('set_state', (handler.INDEX, M,L,L,L,L,V,V,V), (handler.ERROR_CODE))
        handler.add_method('get_state', (handler.INDEX), (M,L,L,L,L,V,V,V, handler.ERROR_CODE))
        handler.add_method('set_time_begin', (T), (handler.ERROR_CODE))
        handler.add_method('get_time_begin', (), (T, handler.ERROR_CODE))
        handler.add_method('get_number_of_particles', (), (units.none, handler.ERROR_CODE))
        handler.add_method('get_time_step', (), (T, handler.ERROR_CODE))
        handler.add_method('set_time_step', (T), (handler.ERROR_CODE))
        handler.add_method('update_force_potential_arrays', (T), (handler.ERROR_CODE))

    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)


