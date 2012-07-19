import numpy

from amuse.units import nbody_system
from amuse.units import units
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class PhiGRAPEInterface(
    CodeInterface, 
    LiteratureReferencesMixIn, 
    GravitationalDynamicsInterface, 
    StoppingConditionInterface,
    GravityFieldInterface):
    """
        .. [#] Harfst, S., Gualandris, A., Merritt, D., Spurzem, R., Portegies Zwart, S., & Berczik, P. 2007, New Astronomy, 12, 357
    """

    MODE_G6LIB = 'g6lib'
    MODE_GPU   = 'gpu'
    MODE_GRAPE = 'grape'

    def __init__(self, mode = MODE_G6LIB, number_of_workers = 1, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(mode, number_of_workers), number_of_workers = number_of_workers, **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self, mode, number_of_workers):
        if number_of_workers > 1:
            if mode == self.MODE_G6LIB:
                return 'phigrape_worker_mpi'
            elif mode == self.MODE_GPU:
                return 'phigrape_worker_gpu'
            elif mode == self.MODE_GRAPE:
                return 'phigrape_worker_grape'
            else:
                return 'phigrape_worker_mpi'
        else:
            if mode == self.MODE_G6LIB:
                return 'phigrape_worker'
            elif mode == self.MODE_GPU:
                return 'phigrape_worker_gpu'
            elif mode == self.MODE_GRAPE:
                return 'phigrape_worker_grape'
            else:
                return 'phigrape_worker'

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

class PhiGRAPEInterfaceGL(PhiGRAPEInterface):

    def __init__(self, mode = PhiGRAPEInterface.MODE_G6LIB):
        PhiGRAPEInterface.__init__(self, mode = mode)

    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()
        return function

    def name_of_the_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'phigrape_worker_gl'
        if mode == self.MODE_GPU:
            return 'phigrape_worker_gl_gpu'
        if mode == self.MODE_GRAPE:
            return 'phigrape_worker_gl_grape'
        else:
            return 'phigrape_worker_gl'



class PhiGRAPE(GravitationalDynamics, GravityFieldCode):

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
                                       
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_eta",
            "set_eta1",
            "timestep_parameter",
            "timestep parameter",
            default_value = 0.02
        )

        object.add_method_parameter(
            "get_eta_s",
            "set_eta_s",
            "initial_timestep_parameter",
            "parameter to determine the initial timestep",
            default_value = 0.01
        )
        
        object.add_method_parameter(
            "get_initialize_once",
            "set_initialize_once",
            "initialize_gpu_once",
            "set to 1 if the gpu must only be initialized once, 0 if it can be initialized for every call\nIf you want to run multiple instances of the code on the same gpu this parameter needs to be 0 (default)",
            default_value = 0
        )
        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_eta",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eta1",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_eta_s",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            "set_eta_s",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_initialize_once",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            "set_initialize_once",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        self.stopping_conditions.define_methods(object)

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)
        
