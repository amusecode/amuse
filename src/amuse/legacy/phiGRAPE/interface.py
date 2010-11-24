import numpy

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics

class PhiGRAPEInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface, StoppingConditionInterface):
    """
        .. [#] Harfst, S., Gualandris, A., Merritt, D., Spurzem, R., Portegies Zwart, S., & Berczik, P. 2007, New Astronomy, 12, 357
    """

    MODE_G6LIB = 'g6lib'
    MODE_GPU   = 'gpu'
    MODE_GRAPE = 'grape'

    def __init__(self, mode = MODE_G6LIB, **options):
        LegacyInterface.__init__(self, name_of_the_worker = self.name_of_the_muse_worker(mode), **options)
        LiteratureRefs.__init__(self)

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'phigrape_worker'
        elif mode == self.MODE_GPU:
            return 'phigrape_worker_gpu'
        elif mode == self.MODE_GRAPE:
            return 'phigrape_worker_grape'
        

    

    


    def initialize_particles(self, time):
        return self.commit_particles()

    def reinitialize_particles():
        return self.recommit_particles()

    @legacy_function
    def get_time_step():
        function = LegacyFunctionSpecification()
        function.result_type = 'd'
        return function

    @legacy_function
    def set_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('etas', dtype='d', direction=function.IN)
        function.addParameter('eta', dtype='d', direction=function.IN)
        return function

    @legacy_function
    def set_eta_s():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta1():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta_s():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def set_initialize_once():
        """
        Sets the initialize GPU/GRAPE once parameter. When
        set the GPU will be initialized during the
        :func:`commit_parameters` call and released
        during the :func:`cleanup_code` call.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='i', direction=function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_initialize_once():
        """
        Returns the current value of the initialize once
        parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='i', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_energy_error():
        function = LegacyFunctionSpecification()
        function.result_type = 'd'
        return function

    @legacy_function
    def find_colliding_secondary():
        function = LegacyFunctionSpecification()
        function.addParameter('id1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_gravity_at_point():
        """
        Determine the gravitational force on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
            description = "The smoothing parameter")
        function.addParameter('x', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('y', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('z', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('forcex', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcey', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcez', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('number_of_points', dtype='int32', direction=function.LENGTH,
            description = "number of points to determine the force for")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
         0 - OK
            Force could be calculated
        -1 - ERROR
            No force calculation supported
        """
        return function
        
    
    @legacy_function
    def get_potential_at_point():
        """
        Determine the potential on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
         description = "The smoothing factor, may be ignored by the code")
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('phi', dtype='float64', direction=function.OUT)
        function.addParameter('number_of_points', dtype='int32', direction=function.LENGTH,
            description = "number of points to determine the force for")
        function.must_handle_array = True
        function.result_type = 'int32'
        return function


    

class PhiGRAPEInterfaceGL(PhiGRAPEInterface):

    def __init__(self, mode = PhiGRAPEInterface.MODE_G6LIB):
        PhiGRAPEInterface.__init__(self, mode = mode)

    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()
        return function

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'glworker_code'
        if mode == self.MODE_GPU:
            return 'glworker_code_gpu'
        if mode == self.MODE_GRAPE:
            return 'glworker_code_grape'



class PhiGRAPE(GravitationalDynamics):

    def __init__(self, convert_nbody = None, mode = PhiGRAPEInterface.MODE_G6LIB, use_gl = False, **options):
        nbody_interface = None
        if use_gl:
            nbody_interface = PhiGRAPEInterfaceGL(mode, **options)
        else:
            nbody_interface = PhiGRAPEInterface(mode, **options)

        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            nbody_system.length * nbody_system.length,
            0.0 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_eta",
            "set_eta1",
            "timestep_parameter",
            "timestep parameter",
            units.none,
            0.02 |  units.none
        )

        object.add_method_parameter(
            "get_eta_s",
            "set_eta_s",
            "initial_timestep_parameter",
            "parameter to determine the initial timestep",
            units.none,
            0.01 |  units.none
        )
        
        object.add_method_parameter(
            "get_initialize_once",
            "set_initialize_once",
            "initialize_gpu_once",
            "set to 1 if the gpu must only be initialized once, 0 if it can be initialized for every call\nIf you want to run multiple instances of the code on the same gpu this parameter needs to be 0 (default)",
            units.none,
            0 |  units.none
        )

        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )

        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )

        self.stopping_conditions.define_methods(object)

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object, 'particles')
        
    def stop(self):
        self.cleanup_code()
        self.overridden().stop()
        
