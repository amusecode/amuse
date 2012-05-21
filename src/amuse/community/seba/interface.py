from amuse.community import *
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset

class SeBaInterface(CodeInterface, LiteratureReferencesMixIn):
    
    """
    Stellar evolution is performed by the rapid single-star evolution
    and binary evolution using SeBa.This is a package of
    semi-analytical formulae which covers all phases of evolution from
    the zero-age main-sequence up to and including remnant phases. It
    is valid for masses in the range 0.01-1000 Msun with variable
    metallicity.  SeBa includes prescriptions for mass loss by stellar
    winds, supernova and supports binary evolution.
    
        .. [#] Portegies Zwart S.F. & Verbunt F., 1996, A&A, 309, 179:
        .. [#] ... "Population synthesis of high-mass binaries"
        .. [#] Toonen, S., Nelemans, G., Portegies Zwart S.F., 2012 submitted to A&A (arXiv 1101.2787)
        .. [#] ... "New population synthesis model: Preliminary results for close double white dwarf populations"
    """

    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="seba_worker", **options)
        LiteratureReferencesMixIn.__init__(self)

    
    @legacy_function
    def evolve_star():
        function = LegacyFunctionSpecification()  
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.addParameter('endtime', dtype='float64', direction=function.IN)
        function.addParameter('metal', dtype='float64', direction=function.IN)
        function.addParameter('resulttime', dtype='float64', direction=function.OUT)
        function.addParameter('end_mass', dtype='float64', direction=function.OUT)
        function.addParameter('end_radius', dtype='float64', direction=function.OUT)
        function.addParameter('end_luminosity', dtype='float64', direction=function.OUT)
        function.addParameter('end_temperature', dtype='float64', direction=function.OUT)
        function.addParameter('time_step', dtype='float64', direction=function.OUT)
        function.addParameter('stellar_type', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

class SeBaParticles(Particles):
    
    def __init__(self, code_interface, storage = None):
        Particles.__init__(self, storage = storage)
        self._private.code_interface = code_interface 
    
    def add_particles_to_store(self, keys, attributes = [], values = []):
        if len(keys) == 0:
            return
            
        all_attributes = []
        all_attributes.extend(attributes)
        all_values = []
        all_values.extend(values)
        
        given_attributes = set(attributes)
        
        if not "initial_mass" in given_attributes:
            index_of_mass_attibute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attibute] * 1.0)
        
        super(SeBaParticles, self).add_particles_to_store(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface._evolve_particles(added_particles, 1e-08 | units.yr)

class SeBa(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SeBaInterface(**options), **options)
    
        self.model_time = 0.0 | units.yr

    def define_parameters(self, object):
    
        object.add_caching_parameter(
        "initialize",
        "z_in",
        "metallicity",
        "Metallicity of all stars",
        0.02
    )

    def define_methods(self, object):
        object.add_method(
            "evolve_star",
            (units.MSun, units.Myr, units.none),
            (units.Myr, units.MSun, units.RSun, units.LSun, units.K, units.Myr,units.stellar_type, object.ERROR_CODE)
        )
    
    
    def define_particle_sets(self, object):
        object.define_inmemory_set('particles', SeBaParticles)
        
    def _evolve_particles(self, particles, end_time):
        attributes = (
            "age",
            "mass",
            "radius",
            "luminosity",
            "temperature",
            "time_step",
            "stellar_type",
        )
        
        result = self.evolve_star(
            particles.initial_mass,
            end_time.as_vector_with_length(len(particles)),
            self.parameters.metallicity,
        )
        
        particles.set_values_in_store(particles.get_all_keys_in_store(), attributes, result)

    def evolve_model(self, end_time = None):
        if end_time is None:
            raise exceptions.CodeException("Cannot determine end_time from the code yet, so end_time must be provided!")
            
        self._evolve_particles(self.particles, end_time)
        self.model_time = end_time

    def commit_particles(self):
        self.evolve_model(1|units.yr)
