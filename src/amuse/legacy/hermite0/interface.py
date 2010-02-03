from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface.gd import NBodyGravitationalDynamicsBinding
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.data.core import Particles,ParticlesWithUnitsConverted
from amuse.support.data.binding import InCodeAttributeStorage
from amuse.support.data import binding

class HermiteInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """ 
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
        
        

class HermiteInCodeAttributeStorage(InCodeAttributeStorage):
    new_particle_method = binding.NewParticleMethod(
        "new_particle", 
        (
            ("mass", "mass", nbody_system.mass),
            ("radius", "radius", nbody_system.length),
            ("x", "x", nbody_system.length),
            ("y", "y", nbody_system.length),
            ("z", "z", nbody_system.length),
            ("vx", "vx", nbody_system.speed),
            ("vy", "vy", nbody_system.speed),
            ("vz", "vz", nbody_system.speed),
        )
    )
    
    getters = (
        binding.ParticleGetAttributesMethod(
            "get_state",
            (
                ("mass", "mass", nbody_system.mass),
                ("radius", "radius", nbody_system.length),
                ("x", "x", nbody_system.length),
                ("y", "y", nbody_system.length),
                ("z", "z", nbody_system.length),
                ("vx", "vx", nbody_system.speed),
                ("vy", "vy", nbody_system.speed),
                ("vz", "vz", nbody_system.speed),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_mass",
            (
                ("mass", "mass", nbody_system.mass),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_radius",
            (
                ("radius", "radius", nbody_system.length),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_position",
            (
                ("x", "x", nbody_system.length),
                ("y", "y", nbody_system.length),
                ("z", "z", nbody_system.length),
            )
        ),
    )
    
    
    setters = (
        binding.ParticleSetAttributesMethod(
            "set_mass",
            (
                ("mass", "mass", nbody_system.mass),
            )
        ),
        binding.ParticleSetAttributesMethod(
            "set_radius",
            (
                ("radius", "radius", nbody_system.length),
            )
        ),
        binding.ParticleSetAttributesMethod(
            "set_position",
            (
                ("x", "x", nbody_system.length),
                ("y", "y", nbody_system.length),
                ("z", "z", nbody_system.length),
            )
        ),
        
    )
            
class HermiteBinding(NBodyGravitationalDynamicsBinding):
    parameter_definitions = [
        parameters.ModuleAttributeParameterDefinition(
            "eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        )
    ]
    
    attribute_definitions = [ ]

    def __init__(self, convert_nbody = None):
        NBodyGravitationalDynamicsBinding.__init__(self, convert_nbody)
        
        self.nbody_particles = Particles(storage = HermiteInCodeAttributeStorage(self))
        self.particles = ParticlesWithUnitsConverted(self.nbody_particles, self.convert_nbody.as_converter_from_si_to_nbody())
    
    def current_model_time(self):
        return self.convert_nbody.to_si( self.t | nbody_system.time)
            
    
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time))
        return result
        
class Hermite(HermiteInterface, HermiteBinding):
    """  
    N-body integration module with shared but variable time step
    (the same for all particles but its size changing in time),
    using the Hermite integration scheme.

    
    .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """	
    
    def __init__(self, convert_nbody = None):
        HermiteInterface.__init__(self)
        HermiteBinding.__init__(self, convert_nbody)
        

