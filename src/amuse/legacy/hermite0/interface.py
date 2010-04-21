from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.support.lit import LiteratureRefs


class HermiteInterface(LegacyInterface, LiteratureRefs, GravitationalDynamicsInterface):
    """
    N-body integration module with shared but variable time step
    (the same for all particles but its size changing in time),
    using the Hermite integration scheme.


    .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """
    include_headers = ['worker_code.h', 'parameters.h', 'local.h']

    t = legacy_global(name='t',id=20,dtype='d')
    dt_param = legacy_global(name='dt_param',id=21,dtype='d')
    dt_dia = legacy_global(name='dt_dia',id=22,dtype='d')
    eps2 = legacy_global(name='eps2',id=23,dtype='d')
    flag_collision = legacy_global(name='flag_collision',id=24,dtype='i')

    def __init__(self, **options):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code", **options)
        LiteratureRefs.__init__(self)

    def setup_module(self):
        self.commit_parameters()
        self.commit_particles()


    def cleanup_module(self):
        self.cleanup_code

    def reinitialize_particles(self):
        self.recommit_particles()

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_particle', dtype='int32', direction=function.IN,
            description = "Index of particle to delete")
        function.result_type = 'int32'
        function.resutl_doc = """
        0 - OK
           The particle was deleted
        """
        return function

class HermiteDoc(object):

    def __get__(self, instance, owner):
        return instance.parameters.__doc__

class Hermite(GravitationalDynamics):

    __doc__ = HermiteDoc()

    def __init__(self, convert_nbody = None, **options):

        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()


        legacy_interface = HermiteInterface(**options)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_attribute_parameter(
            "eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            nbody_system.length * nbody_system.length,
            0.3 | nbody_system.length * nbody_system.length
        )


    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )

        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )
