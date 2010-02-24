from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.support.lit import LiteratureRefs


class HermiteInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """  
    N-body integration module with shared but variable time step
    (the same for all particles but its size changing in time),
    using the Hermite integration scheme.

    
    .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """
    include_headers = ['hermite_code.h', 'parameters.h', 'local.h']

    t = legacy_global(name='t',id=20,dtype='d')
    dt_param = legacy_global(name='dt_param',id=21,dtype='d')
    dt_dia = legacy_global(name='dt_dia',id=22,dtype='d')
    eps2 = legacy_global(name='eps2',id=23,dtype='d')
    flag_collision = legacy_global(name='flag_collision',id=24,dtype='i')

    def __init__(self, convert_nbody = None):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")       
        LiteratureRefs.__init__(self)

    @legacy_function
    def setup_module():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

    @legacy_function
    def cleanup_module():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function

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
        
class Hermite(GravitationalDynamicsInterface):
    
    
    def __init__(self, convert_nbody = None):
        
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
       
        
        legacy_interface = HermiteInterface()
        
        GravitationalDynamicsInterface.__init__(
            self,
            legacy_interface,
            convert_nbody,
        )     
            
    def setup_parameters(self, object):
        object.add_attribute_parameter(
            "eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.3 | nbody_system.length * nbody_system.length
        )     
           
    def setup_methods(self, object):
        GravitationalDynamicsInterface.setup_methods(self, object)
        
        object.add_method(
            "new_particle", 
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ), 
            ( 
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "delete_particle", 
            (
                object.NO_UNIT,
            ), 
            ( 
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state", 
            (
                object.NO_UNIT,
            ), 
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state", 
            (
                object.NO_UNIT,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ), 
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_mass", 
            (
                object.NO_UNIT,
                nbody_system.mass,
            ), 
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass", 
            (
                object.NO_UNIT,
            ), 
            (
                nbody_system.mass,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_radius", 
            (
                object.NO_UNIT,
                nbody_system.length,
            ), 
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_radius", 
            (
                object.NO_UNIT,
            ), 
            (
                nbody_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_position", 
            (
                object.NO_UNIT,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ), 
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position", 
            (
                object.NO_UNIT,
            ), 
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                object.ERROR_CODE
            )
        )
    
    def setup_all_particles(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_state')
        object.add_getter('particles', 'get_state')
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_setter('particles', 'set_position')
        object.add_getter('particles', 'get_position')

