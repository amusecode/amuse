from amuse.legacy import *
from amuse.legacy.interface.se import StellarEvolution
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.data.binding import InterfaceWithParametersBinding, InterfaceWithObjectsBinding

from amuse.support.data.core import Particles
from amuse.support.data.binding import InCodeAttributeStorage2
from amuse.support.data import binding



import os

class EVtwinInterface(LegacyInterface, LiteratureRefs, StellarEvolution):
    """
    Need to have docs
    """
    use_modules = ['twin_library_v2']
    
    def __init__(self):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")
        LiteratureRefs.__init__(self)
         
    @property
    def default_path_to_ev_database(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'src')
        
    

    @legacy_function
    def get_maximum_number_of_stars():
        """
        Retrieve the maximum number of stars that can be
        handled by this instance.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('metallicity', dtype='int32', direction=function.OUT,
            description = "The current value of the maximum number of stars")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of was retrieved
        """
        return function
        
    
    @legacy_function
    def set_maximum_number_of_stars():
        """
        Update the maximum number of stars that can be
        handled by this instance. Need to set this number
        before calling :method:`initialize_code`. Cannot be
        changed once initialize_code has been called.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('metallicity', dtype='int32', direction=function.IN,
            description = "The new value of the maximum number of stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            The code cannot update the maximum number of stars
        """
        return function
        
    
    @legacy_function
    def set_ev_path():
        """
        Update the path to the EVtwin database.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "Name of the the directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function
        
    @legacy_function
    def set_init_dat_name():
        """
        Update name of the init.dat file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "File in the evpath directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
        
    @legacy_function
    def set_init_run_name():
        """
        Update name of the init.run file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "File in the evpath directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
        
class EVtwinInCodeAttributeStorage(InCodeAttributeStorage2):
    name_of_delete_particle = "delete_star"
    
    new_particle_method = binding.NewParticleMethod(
        "new_particle", 
        (
            ("mass", "mass", units.MSun),
            ("radius", "radius", units.RSun),
        )
    )
    
    getters = (
        binding.ParticleGetAttributesMethod(
            "get_mass",
            (
                ("mass", "mass", units.MSun),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_radius",
            (
                ("radius", "radius", units.RSun),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_stellar_type",
            (
                ("type", "stellar_type", units.stellar_type),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_age",
            (
                ("age", "age", units.Myr),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_luminosity",
            (
                ("luminosity", "luminosity", units.LSun),
            )
        ),
        
    )
    
class EVtwinBinding(InterfaceWithParametersBinding, InterfaceWithObjectsBinding):
    
    def __init__(self):
        InterfaceWithParametersBinding.__init__(self)
        InterfaceWithObjectsBinding.__init__(self)
        
        self.particles = Particles()
        self.particles._private.attribute_storage = EVtwinInCodeAttributeStorage(self)
        
    parameter_definitions = [
        parameters.ModuleMethodParameterDefinition_Next(
            "get_maximum_number_of_stars",
            "set_maximum_number_of_stars",
            "maximum_number_of_stars", 
            "Maximum number of stars that can be allocated", 
            units.none, 
            10 | units.none
        ),
        
        parameters.ModuleMethodParameterDefinition_Next(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stats", 
            units.percentage, 
            0.02 | units.percentage
        ),
        
        
        parameters.ModuleMethodParameterDefinition_Next(
            None,
            "set_ev_path",
            "path_to_data", 
            "Path to the data directory", 
            units.string, 
            "src" | units.string
        ),
        
    ]
    
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.set_ev_path(self.default_path_to_ev_database)
        self.initialize_code()
        
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def evolve_model(self, end_time = None):
        if end_time is None:
            self._evolve_particles(self.particles)
            return
            
        end_times = end_time.as_vector_with_length(len(self.particles))
        
        particles_set = particles.to_set()
        while len(particles_set) > 0:
            self._evolve_particles(particles_set)
            particles_set = particles_set.select(lambda x : x < end_time, ["age"])
            
                
    def _evolve_particles(self, particles):
        print "EVO1"
        for particle in particles:
            index = self.particles._private.attribute_storage.mapping_from_particle_key_to_index_in_the_code[particle.key]
            
            errorcode = self.evolve(index)
            print index, errorcode
            if errorcode != 0:
                raise Exception("Error during evolve of particle")
            print particle
        
    def current_model_time(self):
        return self._current_model_time
        
class EVtwin(EVtwinInterface, EVtwinBinding):
    """  
    """
    
    def __init__(self):
        EVtwinInterface.__init__(self)
        EVtwinBinding.__init__(self)
        
