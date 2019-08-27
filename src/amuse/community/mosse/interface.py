import numpy
from operator import itemgetter
from amuse.community import *
from amuse.units import units
from amuse.units import constants
from amuse.support.interface import InCodeComponentImplementation
from amuse.community.interface import common

from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset
class MOSSEInterface(CodeInterface, common.CommonCodeInterface , LiteratureReferencesMixIn): 
    """
    MOSSE (Massive Object in SSE) is an updated version of the **rapid** binary-star 
    evolution (SSE) algorithm. With respect to SSE, the major upgrades are that MOSSE 
    includes up-to-date equations for metal-dependent stellar winds and new prescriptions 
    for core-collapse supernova explosion (SNe). Moreover, MOSSE includes the dependence 
    of stellar winds on the Eddington factor: if a star approaches the Eddington limit 
    stellar winds become almost insensitive to metallicity. MOSSE includes also the effects 
    of Pulsation Pair Instability SNe and Pair Instability SNe, it has more precise 
    formulas to compute the core radii and it can use a different velocity distribution 
    to calculate the velocity kick due to electron-capture SNe. More details about MOSSE can
    be found in paper:

        .. [#] Nicola Giacobbo, Michela Mapelli & Mario Spera, 2018, MNRAS, 474, 2959:
        .. [#] ... Merging black hole binaries: the effects of progenitor s metallicity, mass-loss rate and Eddington factor

    Instead, the details about BSE can be found in the BSE paper:
        ..  Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897:
        ..  ... Evolution of binary stars and the effect of tides on binary populations
        ..  Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543:
        ..  ... Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity
    """
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="mosse_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
    
    @legacy_function   
    def initialize():
        function = LegacyFunctionSpecification()  
        function.addParameter('z_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('neta_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('bwind_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('hewind_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('sigma1_in', dtype='d', direction=function.IN, unit = units.km / units.s)
        function.addParameter('sigma2_in', dtype='d', direction=function.IN, unit = units.km / units.s)
        function.addParameter('ifflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('wdflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('bhflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('nsflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('piflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('mxns_in', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('idum_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('pts1_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('pts2_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('pts3_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('status', dtype='i', direction=function.OUT, unit = NO_UNIT)
        return function
        
    @legacy_function     
    def evolve_star():
        function = LegacyFunctionSpecification()  
        function.name = 'evolve0'
        function.can_handle_array = True 
        function.addParameter('kw', dtype='i', direction=function.INOUT, unit = units.stellar_type)
        function.addParameter('mass', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('mt', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('r', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('lum', dtype='d', direction=function.INOUT, unit = units.LSun)
        function.addParameter('mc', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('rc', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('menv', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('renv', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('ospin', dtype='d', direction=function.INOUT, unit = units.yr**-1)
        function.addParameter('epoch', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('tm', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('tphys', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('tphysf', dtype='d', direction=function.INOUT, unit = units.Myr)
        return function
        
    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True  
        function.addParameter('kw', dtype='i', direction=function.IN, unit = units.stellar_type)
        function.addParameter('mass', dtype='d', direction=function.IN, unit =  units.MSun)
        function.addParameter('age', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('mt', dtype='d', direction=function.IN, unit =  units.MSun)
        function.addParameter('tm', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('epoch', dtype='d', direction=function.IN, unit = units.Myr)
        function.addParameter('dt', dtype='d', direction=function.OUT, unit = units.Myr)
       
        return function
        
    @legacy_function      
    def get_mass_loss_wind():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True  
        function.addParameter('kw', dtype='i', direction=function.IN, unit = units.stellar_type)
        function.addParameter('lum', dtype='d', direction=function.IN, unit = units.LSun)
        function.addParameter('r', dtype='d', direction=function.IN, unit = units.RSun)
        function.addParameter('mt', dtype='d', direction=function.IN, unit =  units.MSun)
        function.addParameter('mc', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('mlout', dtype='d', direction=function.OUT, unit = units.MSun/units.yr)
       
        return function
        

    @legacy_function      
    def get_gyration_radius():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True  
        function.addParameter('kw', dtype='i', direction=function.IN, unit = units.stellar_type)
        function.addParameter('mass', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('mt', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('r', dtype='d', direction=function.IN, unit = units.RSun)
        function.addParameter('lum', dtype='d', direction=function.IN, unit = units.LSun)
        function.addParameter('epoch', dtype='d', direction=function.IN, unit = units.Myr)
        function.addParameter('tm', dtype='d', direction=function.IN, unit = units.Myr)
        function.addParameter('tphys', dtype='d', direction=function.IN, unit = units.Myr)
        function.addParameter('rg', dtype='d', direction=function.OUT, unit = NO_UNIT)
       
        return function
        
    def initialize_code(self):
        return 0
        
    def commit_parameters(self):
        return 0
        
    def recommit_parameters(self):
        return 0
        
    def cleanup_code(self):
        return 0
        
    def commit_particles(self):
        return 0
        
    
        
class MOSSEParticles(Particles):
    
    def __init__(self, code_interface, storage = None):
        Particles.__init__(self, storage = storage)
        self._private.code_interface = code_interface 
        self.add_calculated_attribute("temperature", self.calculate_effective_temperature, ["luminosity", "radius"])
        self.add_function_attribute("evolve_one_step", self.particleset_evolve_one_step, self.evolve_one_step)
        self.add_function_attribute("evolve_for", self.particleset_evolve_for, self.evolve_for)
    
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
            "CO_core_mass": 0 | units.MSun,
            "core_radius": 0 | units.RSun,
            "convective_envelope_mass": 0 | units.MSun,
            "convective_envelope_radius": 0 | units.RSun,
            "epoch": 0 | units.Myr,
            "spin": 0 | units.yr**-1,
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
        
        super(MOSSEParticles, self).add_particles_to_store(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface._evolve_particles(added_particles, 0 | units.yr)
    
    def evolve_one_step(self, particles, subset):
        self._private.code_interface._evolve_particles(subset.as_set(), subset.age + subset.time_step)
    
    def particleset_evolve_one_step(self, particles):
        self._private.code_interface._evolve_particles(particles, particles.age + particles.time_step)
    
    def evolve_for(self, particles, subset, delta_time):
        self._private.code_interface._evolve_particles(subset.as_set(), subset.age + delta_time)
    
    def particleset_evolve_for(self, particles, delta_time):
        self._private.code_interface._evolve_particles(particles, particles.age + delta_time)
    
    
    def get_defined_attribute_names(self):
        return ["mass", "radius"]
        
        
    

    

    
        
        

class MOSSE(common.CommonCode):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MOSSEInterface(**options), **options)
        
        self.model_time = 0.0 | units.yr
        
    
    def define_parameters(self, handler):
    
        handler.add_caching_parameter(
            "initialize",
            "z_in",
            "metallicity",
            "Metallicity of all stars",
            0.02
        )
                
        handler.add_caching_parameter(
            "initialize",
            "neta_in",
            "reimers_mass_loss_coefficient",
            "Reimers mass-loss coefficient (neta*4x10^-13; 0.5 normally)",
            0.5
        )
        
        handler.add_caching_parameter(
            "initialize",
            "bwind_in",
            "binary_enhanced_mass_loss_parameter",
            "The binary enhanced mass loss parameter (inactive for single).",
            0.0
        )
        
        handler.add_caching_parameter(
            "initialize",
            "hewind_in",
            "helium_star_mass_loss_factor",
            "Helium star mass loss factor",
            1.0
        )
        
        handler.add_caching_parameter(
            "initialize",
            "sigma1_in",
            "SN_kick_speed_dispersion_ICS",
            "The dispersion in the Maxwellian for the SN kick speed (265 km/s from Hobbs et al. 2005).",
            265.0 | units.km / units.s
        )
        
        handler.add_caching_parameter(
            "initialize",
            "sigma2_in",
            "SN_kick_speed_dispersion_ECS",
            "The dispersion in the Maxwellian for the SN kick speed (7 km/s).",
            7.0 | units.km / units.s
        )
        
        handler.add_caching_parameter(
            "initialize",
            "ifflag_in",
            "white_dwarf_IFMR_flag",
            "ifflag > 0 uses white dwarf IFMR (initial-final mass relation) of HPE, 1995, MNRAS, 272, 800 (0).",
            0
        )
        
        handler.add_caching_parameter(
            "initialize",
            "wdflag_in",
            "white_dwarf_cooling_flag",
            "wdflag > 0 uses modified-Mestel cooling for WDs (0).",
            1
        )
        
        handler.add_caching_parameter(
            "initialize",
            "bhflag_in",
            "black_hole_kick_flag",
            "bhflag > 0 allows velocity kick at BH formation (0).",
            1
        )
        
        handler.add_caching_parameter(
            "initialize",
            "nsflag_in",
            "neutron_star_mass_flag",
            "nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1).",
            3
        )
        
        handler.add_caching_parameter(
            "initialize",
            "piflag_in",
            "pair_instability_flag",
            "piflag > 0  activates the PPISN and PISN from Spera & Mapelli 2017, MNRAS, 470, 4739 (1).",
            1
        )

        handler.add_caching_parameter(
            "initialize",
            "mxns_in",
            "maximum_neutron_star_mass",
            "The maximum neutron star mass (1.8, nsflag=0; 3.0, nsflag=1).",
            3.0 | units.MSun
        )
        
        handler.add_caching_parameter(
            "initialize",
            "idum_in",
            "SN_kick_random_seed",
            "The random number seed used in the kick routine.",
            29769
        )
            
        handler.add_caching_parameter(
            "initialize",
            "pts1_in",
            "fractional_time_step_1",
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: MS (0.05)",
            0.05
        )
        
        handler.add_caching_parameter(
            "initialize",
            "pts2_in",
            "fractional_time_step_2",
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: GB, CHeB, AGB, HeGB (0.01)",
            0.01
        )
        
        handler.add_caching_parameter(
            "initialize",
            "pts3_in",
            "fractional_time_step_3",
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: HG, HeMS (0.02)",
            0.02
        )
    
    def define_state(self, handler):
        common.CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','RUN','commit_parameters')
        handler.add_method('RUN', 'evolve_star')
        
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('RUN', 'before_set_parameter')
    
         
    def define_particle_sets(self, handler):
        handler.define_inmemory_set('particles', MOSSEParticles)
        
        handler.add_attribute(
            'particles',
            'time_step', 
            'get_time_step', 
            ('stellar_type', 'initial_mass', 'age', 
             'mass', 'main_sequence_lifetime', 'epoch')
        )
        
        handler.add_attribute(
            'particles',
            'mass_loss_wind', 
            'get_mass_loss_wind', 
            ('stellar_type', 'luminosity', 
             'radius', 'mass', 
             'CO_core_mass')
        )
        
        handler.add_attribute(
            'particles',
            'gyration_radius', 
            'get_gyration_radius', 
            ('stellar_type', 'initial_mass','mass','radius',
             'luminosity','epoch','main_sequence_lifetime',
             'age')
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
            "convective_envelope_mass", 
            "convective_envelope_radius", 
            "spin", 
            "epoch",
            "main_sequence_lifetime",
            "age", 
            "end_time"
        )
        
        result = self.evolve_star(
            particles.stellar_type,
            particles.initial_mass,
            particles.mass,
            particles.radius,
            particles.luminosity,
            particles.core_mass,
            particles.core_radius,
            particles.convective_envelope_mass,
            particles.convective_envelope_radius,
            particles.spin,
            particles.epoch,
            particles.main_sequence_lifetime,
            particles.age,
            end_time.as_vector_with_length(len(particles)))
        
        # For helium (and helium exhausted) stars, the core mass returned is actually the CO core mass
        type = result[0].value_in(units.stellar_type)
        helium_star_selection = (type > 6) & (type < 16) & (type != 10)
        helium_stars = particles[helium_star_selection]
        other_stars = particles - helium_stars
        if len(helium_stars):
            helium_stars_results = [sub[helium_star_selection] for sub in result]
            helium_stars_results.append(helium_stars_results[2])
            helium_stars.set_values_in_store(helium_stars.get_all_indices_in_store(), (
                "stellar_type", "initial_mass", "mass", "radius",  "luminosity", 
                "CO_core_mass", 
                "core_radius", "convective_envelope_mass", "convective_envelope_radius", "spin", "epoch",
                "main_sequence_lifetime", "age", "end_time", 
                "core_mass"), helium_stars_results)
        if len(other_stars):
            other_star_selection = numpy.logical_not(helium_star_selection)
            other_stars.set_values_in_store(other_stars.get_all_indices_in_store(), attributes, 
                [sub[other_star_selection] for sub in result])
    
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
