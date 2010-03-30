import numpy

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs

class PhiGRAPEInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface):
    """
        .. [#] Harfst, S., Gualandris, A., Merritt, D., Spurzem, R., Portegies Zwart, S., & Berczik, P. 2007, New Astronomy, 12, 357
    """

    MODE_G6LIB = 'g6lib'
    MODE_GPU   = 'gpu'
    MODE_GRAPE = 'grape'
    MODE_PG    = 'pg'

    def __init__(self, convert_nbody = None, mode = MODE_G6LIB):
        LegacyInterface.__init__(self, name_of_the_worker = self.name_of_the_muse_worker(mode))
        LiteratureRefs.__init__(self)

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'worker_code'
        elif mode == self.MODE_GPU:
            return 'worker_code_gpu'
        elif mode == self.MODE_GRAPE:
            return 'worker_code_grape'
        elif mode == self.MODE_PG:
            return 'worker_code_phantom_grape'

    def setup_module(self):
        return self.initialize_code()

    def cleanup_module(self):
        return self.cleanup_code()


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

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'glworker_code'
        if mode == self.MODE_GPU:
            return 'glworker_code_gpu'
        if mode == self.MODE_GRAPE:
            return 'glworker_code_grape'



class PhiGRAPE(GravitationalDynamics):


    def __init__(self, convert_nbody = None, mode = PhiGRAPEInterface.MODE_G6LIB, use_gl = False):

        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()

        nbody_interface = None
        if use_gl:
            nbody_interface = PhiGRAPEInterfaceGL(mode)
        else:
            nbody_interface = PhiGRAPEInterface(mode)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
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
            units.none ,
            0.02 |  units.none
        )

        object.add_method_parameter(
            "get_eta_s",
            "set_eta_s",
            "initial_timestep_parameter",
            "parameter to determine the initial timestep",
            units.none ,
            0.01 |  units.none
        )



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
