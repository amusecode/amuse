from amuse.community import *
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset
from amuse.community.interface import se

class SeBaInterface(CodeInterface, se.StellarEvolutionInterface, LiteratureReferencesMixIn):
    
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
        
    @legacy_function
    def evolve_system():
        """
        Evolve the model until the given time, or until a stopping condition is set.
        Need to call this evolve_system as evolve_model is overriden in se.StellarEvolution
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "Model time to evolve the code to. The model will be "
                "evolved until this time is reached exactly or just after.")
        function.result_type = 'int32'
        return function

    def evolve_model(self, time):
        return self.evolve_system(time)

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

class SeBa(se.StellarEvolution):

    def __init__(self, **options):
        se.StellarEvolution.__init__(self,  SeBaInterface(**options), **options)
    
        self.model_time = 0.0 | units.yr

    
    def evolve_model(self, end_time):
        return self.evolve_system(end_time);

    def define_methods(self, object):
        se.StellarEvolution.define_methods(self, object)
        
        object.add_method(
            "evolve_star",
            (units.MSun, units.Myr, units.none),
            (units.Myr, units.MSun, units.RSun, units.LSun, units.K, units.Myr,units.stellar_type, object.ERROR_CODE)
        )
        object.add_method(
            "evolve_system",
            (units.Myr,),
            (object.ERROR_CODE,)
        )
    
    def define_particle_sets(self, object):
       
        object.define_set('particles', 'index_of_the_star')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_star')
        
        object.add_getter('particles', 'get_radius', names = ('radius',))
        object.add_getter('particles', 'get_stellar_type', names = ('stellar_type',))
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_getter('particles', 'get_age', names = ('age',))
        object.add_getter('particles', 'get_time_step', names = ('time_step',))
        #object.add_getter('particles', 'get_spin', names = ('spin',))
        object.add_getter('particles', 'get_luminosity', names = ('luminosity',))
        object.add_getter('particles', 'get_temperature', names = ('temperature',))
        object.add_method('particles', 'evolve_one_step')
        object.add_method('particles', 'evolve_for')
    
