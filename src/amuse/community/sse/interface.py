from operator import itemgetter
from amuse.community import *

from amuse.support.units import units
from amuse.support.units import constants

from amuse.support.data.core import Particles, ParticlesSubset
from amuse.support.interface import InCodeComponentImplementation
from amuse.community.interface import common

class SSEInterface(CodeInterface, common.CommonCodeInterface , LiteratureReferencesMixIn): 
    """
    Stellar evolution is performed by the rapid single-star evolution (SSE)
    algorithm. This is a package of analytical formulae fitted to the detailed 
    models of Pols et al. (1998)  that covers all phases of evolution from the 
    zero-age main-sequence up to and including remnant phases. It is valid for 
    masses in the range 0.1-100 Msun and metallicity can be varied. The SSE 
    package contains a prescription for mass loss by stellar winds. It also 
    follows the evolution of rotational angular momentum for the star. Full 
    details can be found in the SSE paper:
    
        .. [#] Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543:
        .. [#] ... "Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity"
    """
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="sse_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
    
    @legacy_function   
    def initialize():
        function = LegacyFunctionSpecification()  
        function.addParameter('z_in', dtype='d', direction=function.IN)
        function.addParameter('neta_in', dtype='d', direction=function.IN)
        function.addParameter('bwind_in', dtype='d', direction=function.IN)
        function.addParameter('hewind_in', dtype='d', direction=function.IN)
        function.addParameter('sigma_in', dtype='d', direction=function.IN)
        function.addParameter('ifflag_in', dtype='i', direction=function.IN)
        function.addParameter('wdflag_in', dtype='i', direction=function.IN)
        function.addParameter('bhflag_in', dtype='i', direction=function.IN)
        function.addParameter('nsflag_in', dtype='i', direction=function.IN)
        function.addParameter('mxns_in', dtype='d', direction=function.IN)
        function.addParameter('pts1_in', dtype='d', direction=function.IN)
        function.addParameter('pts2_in', dtype='d', direction=function.IN)
        function.addParameter('pts3_in', dtype='d', direction=function.IN)
        function.addParameter('status', dtype='i', direction=function.OUT)
        return function
        
    @legacy_function     
    def evolve():
        function = LegacyFunctionSpecification()  
        function.name = 'evolve0'
        function.can_handle_array = True 
        function.addParameter('kw', dtype='i', direction=function.INOUT)
        function.addParameter('mass', dtype='d', direction=function.INOUT)
        function.addParameter('mt', dtype='d', direction=function.INOUT)
        function.addParameter('r', dtype='d', direction=function.INOUT)
        function.addParameter('lum', dtype='d', direction=function.INOUT)
        function.addParameter('mc', dtype='d', direction=function.INOUT)
        function.addParameter('rc', dtype='d', direction=function.INOUT)
        function.addParameter('menv', dtype='d', direction=function.INOUT)
        function.addParameter('renv', dtype='d', direction=function.INOUT)
        function.addParameter('ospin', dtype='d', direction=function.INOUT)
        function.addParameter('epoch', dtype='d', direction=function.INOUT)
        function.addParameter('tm', dtype='d', direction=function.INOUT)
        function.addParameter('tphys', dtype='d', direction=function.INOUT)
        function.addParameter('tphysf', dtype='d', direction=function.INOUT)
        return function
        
    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True  
        function.addParameter('kw', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age', dtype='d', direction=function.IN)
        function.addParameter('mt', dtype='d', direction=function.IN)
        function.addParameter('tm', dtype='d', direction=function.IN)
        function.addParameter('epoch', dtype='d', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.OUT)
        return function
        
    def initialize_code(self):
        pass
        
    def commit_parameters(self):
        pass
        
    def recommit_parameters(self):
        pass
        
    def cleanup_code(self):
        pass
        
    def commit_particles(self):
        pass
        
    
        
class SSEParticles(Particles):
    
    def __init__(self, code_interface, storage = None):
        Particles.__init__(self, storage = storage)
        self._private.code_interface = code_interface 
        self.add_calculated_attribute("temperature", self.calculate_effective_temperature, ["luminosity", "radius"])
    
    def calculate_effective_temperature(self, luminosity, radius):
        return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        if len(keys) == 0:
            return
            
        all_attributes = []
        all_attributes.extend(attributes)
        all_values = []
        all_values.extend(values)
        
        mapping_from_attribute_to_default_value = {
            "stellar_type" : 1 | units.stellar_type,
            "radius": 0 | units.RSun,
            "luminosity":  0 | units.LSun,
            "core_mass": 0 | units.MSun,
            "core_radius": 0 | units.RSun,
            "envelope_mass": 0 | units.MSun,
            "envelope_radius": 0 | units.RSun,
            "epoch": 0 | units.Myr,
            "spin": 0 | units.none,
            "main_sequence_lifetime": 0 | units.Myr,
            "age": 0 | units.Myr
        }
        
        given_attributes = set(attributes)
        
        if not "initial_mass" in given_attributes:
            index_of_mass_attibute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attibute] * 1.0)
        
        for attribute, default_value in mapping_from_attribute_to_default_value.iteritems():
            if not attribute in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))
        
        super(SSEParticles, self).add_particles_to_store(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface._evolve_particles(added_particles, 1e-08 | units.yr)
    
    def _state_attributes(self):
        return ["mass", "radius"]
        
        
    

    

    
        
        

class SSE(common.CommonCode):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, SSEInterface(**options), **options)
        
        self.parameters.set_defaults()
        self.model_time = 0.0 | units.yr
        
    
    def define_parameters(self, object):
    
        object.add_caching_parameter(
            "initialize",
            "z_in",
            "metallicity",
            "Metallicity of all stars",
            0.02 | units.none
        )
                
        object.add_caching_parameter(
            "initialize",
            "neta_in",
            "reimers_mass_loss_coefficient",
            "Reimers mass-loss coefficient (neta*4x10^-13; 0.5 normally)",
            0.5 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "bwind_in",
            "binary_enhanced_mass_loss_parameter",
            "The binary enhanced mass loss parameter (inactive for single).",
            0.0 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "hewind_in",
            "helium_star_mass_loss_factor",
            "Helium star mass loss factor",
            1.0 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "sigma_in",
            "SN_kick_speed_dispersion",
            "The dispersion in the Maxwellian for the SN kick speed (190 km/s).",
            190.0 | units.km / units.s
        )
        
        object.add_caching_parameter(
            "initialize",
            "ifflag_in",
            "white_dwarf_IFMR_flag", 
            "ifflag > 0 uses white dwarf IFMR (initial-final mass relation) of HPE, 1995, MNRAS, 272, 800 (0).",
            0 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "wdflag_in",
            "white_dwarf_cooling_flag", 
            "wdflag > 0 uses modified-Mestel cooling for WDs (0).",
            1 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "bhflag_in",
            "black_hole_kick_flag",
            "bhflag > 0 allows velocity kick at BH formation (0).",
            0 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "nsflag_in",
            "neutron_star_mass_flag",
            "nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1).",
            1 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "mxns_in",
            "maximum_neutron_star_mass",
            "The maximum neutron star mass (1.8, nsflag=0; 3.0, nsflag=1).",
            3.0 | units.MSun
        )
        
        object.add_caching_parameter(
            "initialize",
            "pts1_in",
            "fractional_time_step_1", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: MS (0.05)",
            0.05 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "pts2_in",
            "fractional_time_step_2", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: GB, CHeB, AGB, HeGB (0.01)",
            0.01 | units.none
        )
        
        object.add_caching_parameter(
            "initialize",
            "pts3_in",
            "fractional_time_step_3", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: HG, HeMS (0.02)",
            0.02 | units.none
        )
    
    def define_state(self, object):
        common.CommonCode.define_state(self, object)
        object.add_transition('INITIALIZED','RUN','commit_parameters')
        object.add_method('RUN', 'evolve')
    
    def define_methods(self, object):
        
        object.add_method( 
            "get_time_step", 
            (
                units.stellar_type, 
                units.MSun, 
                units.Myr, 
                units.MSun,  
                units.Myr, 
                units.Myr
            ),
            (units.Myr,)
        )
        
        object.add_method( 
            "evolve", 
            (
                units.stellar_type, 
                units.MSun, 
                units.MSun, 
                units.RSun, 
                units.LSun,
                units.MSun, 
                units.RSun, 
                units.MSun, 
                units.RSun, 
                units.none,
                units.Myr,
                units.Myr,
                units.Myr,
                units.Myr,
            ),
            (
                units.stellar_type, 
                units.MSun, 
                units.MSun, 
                units.RSun, 
                units.LSun,
                units.MSun, 
                units.RSun, 
                units.MSun, 
                units.RSun,
                units.none, 
                units.Myr,
                units.Myr,
                units.Myr,
                units.Myr,
            )
        )
        
        object.add_method(
            "initialize",
            (
                units.none,
                units.none, 
                units.none, 
                units.none, 
                units.km / units.s,
                units.none, 
                units.none, 
                units.none, 
                units.none, 
                units.MSun,
                units.none, 
                units.none, 
                units.none
            )
        )
         
    def define_particle_sets(self, object):
        object.define_inmemory_set('particles', SSEParticles)
        
        object.add_attribute(
            'particles',
            'time_step', 
            'get_time_step', 
            ('stellar_type', 'initial_mass', 'age', 
             'mass', 'main_sequence_lifetime', 'epoch')
        )
        
    def _evolve_particles(self, particles, end_time):
        attributes = (
            "stellar_type",
            "initial_mass", 
            "mass", 
            "radius",  
            "luminosity", 
            "core_mass", 
            "core_radius",
            "envelope_mass", 
            "envelope_radius", 
            "spin", 
            "epoch",
            "main_sequence_lifetime",
            "age", 
            "end_time"
        )
        
        result = self.evolve(
            particles.stellar_type,
            particles.initial_mass,
            particles.mass,
            particles.radius,
            particles.luminosity,
            particles.core_mass,
            particles.core_radius,
            particles.envelope_mass,
            particles.envelope_radius,
            particles.spin,
            particles.epoch,
            particles.main_sequence_lifetime,
            particles.age,
            end_time.as_vector_with_length(len(particles)))
        
        particles.set_values_in_store(particles.get_all_keys_in_store(), attributes, result)
        
    def evolve_model(self, end_time = None, keep_synchronous = True):
        if not keep_synchronous:
            self._evolve_particles(self.particles, self.particles.time_step + self.particles.age)
            return
        
        if end_time is None:
            end_time = self.model_time + min(self.particles.time_step)
        self._evolve_particles(self.particles, end_time - self.model_time + self.particles.age)
        self.model_time = end_time
    
    def _evolve_model_old(self, end_time = None, keep_synchronous = True):
        """
        This is the old implementation of evolve_model. Even with (keep_synchronous = True) 
        it is unable to evolve all stars to a common age, since it relies on the 
        individual timesteps as determined by the community code. Furthermore, it 
        is not suited to simulations with ongoing star formation, since it evolves 
        newly created stars to the same age as the old stars. 
        """
        if end_time is None:
            if keep_synchronous:
                ages = self.particles.age
                index, min_age = min(enumerate(ages), key=itemgetter(1))
                new_age = min_age + self.particles[index].time_step
                selection = self.particles.select(lambda x : x < new_age, ["age"])
                self._evolve_particles(selection, selection.time_step + selection.age)
                return
            end_time = self.particles.time_step + self.particles.age
            
        self._evolve_particles(self.particles, end_time)
    
        
    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
    
    def initialize_module_with_current_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.commit_parameters()

    def update_time_steps(self):
        pass
