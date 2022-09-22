from amuse.community import *
from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import Quantity
from amuse.community.interface import common

from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset

import numpy

class MOBSEInterface(CodeInterface, common.CommonCodeInterface , LiteratureReferencesMixIn):
    """
    MOBSE (Massive Object in BSE) is an updated version of the **rapid** binary-star 
    evolution (BSE) algorithm. With respect to BSE, the major upgrades are that MOBSE 
    includes up-to-date equations for metal-dependent stellar winds and new prescriptions 
    for core-collapse supernova explosion (SNe). Moreover, MOBSE includes the dependence 
    of stellar winds on the Eddington factor: if a star approaches the Eddington limit 
    stellar winds become almost insensitive to metallicity. MOBSE includes also the effects 
    of Pulsation Pair Instability SNe and Pair Instability SNe, it has more precise 
    formulas to compute the core radii and it can use a different velocity distribution 
    to calculate the velocity kick due to electron-capture SNe. More details about MOBSE can
    be found in paper:

        .. [#] Nicola Giacobbo, Michela Mapelli & Mario Spera, 2018, MNRAS, 474, 2959:
        .. [#] ... Merging black hole binaries: the effects of progenitor s metallicity, mass-loss rate and Eddington factor
        .. [#] [2018MNRAS.474.2959G]

    The details about BSE can be found in the BSE paper:
        ..  Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897:
        ..  ... Evolution of binary stars and the effect of tides on binary populations
        ..  [2002MNRAS.329..897H]
        ..  Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543:
        ..  ... Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity
        ..  [2000MNRAS.315..543H]
    """
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="mobse_worker", **options)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function   
    def initialize():
        function = LegacyFunctionSpecification()  
        function.addParameter('z_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('neta_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('bwind_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('hewind_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('alpha1_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('CElambda_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('ceflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
        function.addParameter('tflag_in', dtype='i', direction=function.IN, unit = NO_UNIT)
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
        function.addParameter('sigma1_in', dtype='d', direction=function.IN, unit = units.km / units.s)
        function.addParameter('sigma2_in', dtype='d', direction=function.IN, unit = units.km / units.s)
        function.addParameter('beta_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('xi_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('acc2_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('epsnov_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('eddfac_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('gamma_in', dtype='d', direction=function.IN, unit = NO_UNIT)
        function.addParameter('status', dtype='i', direction=function.OUT, unit = NO_UNIT)
        return function     

    @legacy_function     
    def evolve_binary():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('type1', dtype='i', direction=function.INOUT, unit = units.stellar_type)
        function.addParameter('type2', dtype='i', direction=function.INOUT, unit = units.stellar_type)
        function.addParameter('initial_mass1', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('initial_mass2', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('mass1', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('mass2', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('radius1', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('radius2', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('luminosity1', dtype='d', direction=function.INOUT, unit = units.LSun)
        function.addParameter('luminosity2', dtype='d', direction=function.INOUT, unit = units.LSun)
        function.addParameter('core_mass1', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('core_mass2', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('core_radius1', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('core_radius2', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('convective_envelope_mass1', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('convective_envelope_mass2', dtype='d', direction=function.INOUT, unit = units.MSun)
        function.addParameter('convective_envelope_radius1', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('convective_envelope_radius2', dtype='d', direction=function.INOUT, unit = units.RSun)
        function.addParameter('spin1', dtype='d', direction=function.INOUT, unit = NO_UNIT)
        function.addParameter('spin2', dtype='d', direction=function.INOUT, unit = NO_UNIT)
        function.addParameter('epoch1', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('epoch2', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('MS_lifetime1', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('MS_lifetime2', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('age', dtype='d', direction=function.INOUT, unit = units.Myr)
        function.addParameter('orbital_period', dtype='d', direction=function.INOUT, unit = units.day)
        function.addParameter('eccentricity', dtype='d', direction=function.INOUT, unit = NO_UNIT)
        function.addParameter('end_time', dtype='d', direction=function.INOUT, unit = units.Myr)
        return function
        
    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('type1', dtype='i', direction=function.IN, unit = units.stellar_type)
        function.addParameter('type2', dtype='i', direction=function.IN, unit = units.stellar_type)
        function.addParameter('initial_mass1', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('initial_mass2', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('mass1', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('mass2', dtype='d', direction=function.IN, unit = units.MSun)
        function.addParameter('MS_lifetime1', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('MS_lifetime2', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('epoch1', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('epoch2', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('age', dtype='d', direction=function.IN, unit =  units.Myr)
        function.addParameter('time_step', dtype='d', direction=function.OUT, unit = units.Myr)
        
        return function
    
    def get_time_step_for_binary(self, binary):
        current_values = {}
        current_values['type1'] = binary.type1.value_in(units.stellar_type)
        current_values['type2'] = binary.type2.value_in(units.stellar_type)
        current_values['initial_mass1'] = binary.initial_mass1.value_in(units.MSun)
        current_values['initial_mass2'] = binary.initial_mass2.value_in(units.MSun)
        current_values['mass1'] = binary.mass1.value_in(units.MSun)
        current_values['mass2'] = binary.mass2.value_in(units.MSun)
        current_values['MS_lifetime1'] = binary.MS_lifetime1.value_in(units.Myr)
        current_values['MS_lifetime2'] = binary.MS_lifetime2.value_in(units.Myr)
        current_values['epoch1'] = binary.epoch1.value_in(units.Myr)
        current_values['epoch2'] = binary.epoch2.value_in(units.Myr)
        current_values['age'] = binary.age.value_in(units.Myr)
        
        result = self.get_time_step(**current_values)
        
        return result | units.Myr
        
        
    def evolve_particle(self, particle, time_end):
        t = particle.current_time
        if particle.stellar_type == 15:
            return
        while t < time_end:
            t0 = t
            t  = t0 + self.get_time_step_for_binary(particle)
            if t > time_end:
                t = time_end
            self.evolve_star(particle, t)
            t1 = particle.current_time
            dt = t1 - t0
            t0 = t1
            if dt.value_in(units.Myr) == 0.0:
                #print t, t0, t1, dt, "BREAK BREAK BREAK!"
                return
            if particle.stellar_type == 15:
                return
    
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

class MOBSEStars(Particles):
    
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
            "convective_envelope_mass": 0 | units.MSun,
            "convective_envelope_radius": 0 | units.RSun,
            "epoch": 0 | units.Myr,
            "spin": 0 | units.none,
            "main_sequence_lifetime": 0 | units.Myr,
            "age": 0 | units.Myr,
            "stellar_type":  0 | units.stellar_type #units.stellar_type("Main Sequence star"),
        }
        
        given_attributes = set(attributes)
        
        if not "initial_mass" in given_attributes:
            index_of_mass_attibute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attibute] * 1.0)
        
        for attribute, default_value in mapping_from_attribute_to_default_value.items():
            if not attribute in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))
        
        super(MOBSEStars, self).add_particles_to_store(keys, all_attributes, all_values)
        
    
    def get_defined_attribute_names(self):
        return ["mass", "radius"]       

class MOBSEBinaries(Particles):
    
    def __init__(self, code_interface, storage = None):
        Particles.__init__(self, storage = storage)
        self._private.code_interface = code_interface 
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        if len(keys) == 0:
            return
            
        given_attributes = set(attributes)
        
        if not "child1" in given_attributes:
            raise Exception("a binary must always have a child1 attribute")
            
        if not "child2" in given_attributes:
            raise Exception("a binary must always have a child2 attribute")
            
        all_attributes = []
        all_values = []
        for attribute, value in zip(attributes, values):
            all_attributes.append(attribute)
            if attribute == 'child1' or attribute == 'child2':
                value = value.copy_with_link_transfer(None, self._private.code_interface.particles)
                all_values.append(value)
            else:
                all_values.append(value)
        
        mapping_from_attribute_to_default_value = {
            "eccentricity": 0.0 | units.none,
            "age": 0 | units.Myr
        }
        
        
        
        for attribute, default_value in mapping_from_attribute_to_default_value.items():
            if not attribute in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))
        
        super(MOBSEBinaries, self).add_particles_to_store(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface._evolve_binaries(added_particles, 1e-08 | units.yr)
        
    def get_defined_attribute_names(self):
        return ["eccentricity", "orbital_period", "age", "child1", "child2"]

class MOBSE(common.CommonCode):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MOBSEInterface(**options), **options)

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
            "alpha1_in",
            "common_envelope_efficiency",
            "The common-envelope efficiency parameter",
            1.0
        )
            
        handler.add_caching_parameter(
            "initialize",
            "CElambda_in",
            "common_envelope_binding_energy_factor",
            "The binding energy factor for common envelope evolution",
            0.1
        )
        
        handler.add_caching_parameter(
            "initialize",
            "ceflag_in",
            "common_envelope_model_flag",
            "ceflag > 0 activates spin-energy correction in common-envelope. ceflag = 3 activates de Kool common-envelope model (0).",
            0
        )
            
        handler.add_caching_parameter(
            "initialize",
            "tflag_in",
            "tidal_circularisation_flag",
            "tflag > 0 activates tidal circularisation (1).",
            1
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
            "nsflag = 3 Delayed SNe model from Fryer et al. 2012, ApJ, 749, 91.",
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
            "beta_in",
            "wind_velocity_factor",
            "The wind velocity factor: proportional to vwind**2 (1/8).",
            0.125
        )
        
        handler.add_caching_parameter(
            "initialize",
            "xi_in",
            "wind_accretion_efficiency",
            "The wind accretion efficiency factor (1.0).",
            1.0
        )
        
        handler.add_caching_parameter(
            "initialize",
            "acc2_in",
            "wind_accretion_factor",
            "The Bondi-Hoyle wind accretion factor (3/2).",
            1.5
        )
        
        handler.add_caching_parameter(
            "initialize",
            "epsnov_in",
            "nova_retained_accreted_matter_fraction",
            "The fraction of accreted matter retained in nova eruption (0.001).",
            0.001
        )
        
        handler.add_caching_parameter(
            "initialize",
            "eddfac_in",
            "Eddington_mass_transfer_limit_factor",
            "The Eddington limit factor for mass transfer (1.0).",
            1.0
        )
        
        handler.add_caching_parameter(
            "initialize",
            "gamma_in",
            "Roche_angular_momentum_factor",
            "The angular momentum factor for mass lost during Roche (-1.0). ",
            -1.0
        )
    
    
    def define_state(self, handler):
        common.CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','RUN','commit_parameters')
        handler.add_method('RUN', 'evolve_binary')
        
        handler.add_method('RUN','before_get_parameter')
        handler.add_method('RUN','before_set_parameter')
    
        
        
         
    def define_particle_sets(self, handler):
        handler.define_inmemory_set('particles', MOBSEStars)
        handler.define_inmemory_set('binaries', MOBSEBinaries)
        
        handler.add_attribute(
            'binaries',
            'time_step', 
            '_get_time_step', 
            ('child1', 'child2', 'age')
            #('child1', 'type2', 
            # 'initial_mass1', 'initial_mass2', 
            # 'mass1', 'mass2',
            # 'MS_lifetime1', 'MS_lifetime2',
            # 'epoch1', 'epoch2',
            #'age')
        )
    
    def _get_time_step(self, child1, child2, age):
        child1 = child1.as_set()
        child2 = child2.as_set()
        return self.get_time_step(
            child1.stellar_type, child2.stellar_type,
            child1.initial_mass, child2.initial_mass,
            child1.mass, child2.mass,
            child1.age, child2.age,
            child1.epoch, child2.epoch,
            age
        )
        
    def orbital_period_to_semi_major_axis(self, orbital_period, mass1, mass2):
        mu = (mass1 + mass2) * constants.G
        return (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)
    
    def semi_major_axis_to_orbital_period(self, semi_major_axis, mass1, mass2):
        mu = (mass1 + mass2) * constants.G
        return 2.0 * numpy.pi * ((semi_major_axis**3/mu)**0.5)
        
    def _evolve_binaries(self, particles, end_time):
        binary_attributes = (
            "age", 
            "semi_major_axis", 
            "eccentricity"
        )
        
        single_attributes = (
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
            "age", 
        )
        
        children1 = particles.child1.as_set()
        children2 = particles.child2.as_set()
        children1_arguments = children1.get_values_in_store(children1.get_all_indices_in_store(), single_attributes)
        children2_arguments = children2.get_values_in_store(children2.get_all_indices_in_store(), single_attributes)
        
        binaries_arguments = particles.get_values_in_store(particles.get_all_indices_in_store(), binary_attributes)
        
        binaries_arguments[1] = self.semi_major_axis_to_orbital_period(binaries_arguments[1] , children1_arguments[2], children2_arguments[2])
        
        arguments = []
        for argument1, argument2 in zip(children1_arguments, children2_arguments):
            arguments.append(argument1)
            arguments.append(argument2)
    
        arguments.extend(binaries_arguments)
        arguments.append(end_time.as_vector_with_length(len(particles)))
        
        result = self.evolve_binary(*arguments)
        
        result[-3] = self.orbital_period_to_semi_major_axis(result[-3] , result[4], result[5])
        
        
        children1_results = []
        children2_results = []
        index = 0
        for dummy in range(len(children1_arguments)):
            children1_results.append(result[index])
            index += 1
            children2_results.append(result[index])
            index += 1
        
        children1.set_values_in_store(children1.get_all_indices_in_store(), single_attributes, children1_results)
        children2.set_values_in_store(children2.get_all_indices_in_store(), single_attributes, children2_results)
        particles.set_values_in_store(particles.get_all_indices_in_store(), binary_attributes, result[index:])
        
        
    def evolve_model(self, end_time = None, keep_synchronous = True):
        if not keep_synchronous:
            self._evolve_binaries(self.binaries, self.binaries.time_step + self.binaries.age)
            return
        
        if end_time is None:
            end_time = self.model_time + min(self.binaries.time_step)
        self._evolve_binaries(self.binaries, end_time - self.model_time + self.binaries.age)
        self.model_time = end_time
    
    def commit_particles(self):
        pass

    def update_time_steps(self):
        pass
        
    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def initialize_module_with_current_parameters(self):
        self.commit_parameters()
        
    def initialize_module_with_default_parameters(self):
        """
        * neta is the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
        * bwind is the binary enhanced mass loss parameter (inactive for single).
        * hewind is a helium star mass loss factor (1.0 normally, inactive). 
        * sigma1 is the dispersion in the Maxwellian for the ICSN kick speed (265 km/s). 
        * sigma2 is the dispersion in the Maxwellian for the ECSN kick speed (7 km/s). 
        *
        * ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
        * wdflag > 0 uses modified-Mestel cooling for WDs (0). 
        * bhflag > 0 allows velocity kick at BH formation (1). 
        * nsflag = 1 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407.",
        *        = 2 Rapid SNe model from Fryer et al. 2012, ApJ, 749, 91.",
        *        = 3 Delayed SNe model from Fryer et al. 2012, ApJ, 749, 91.",
        * mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
        * piflag > 0 activates the PPISN and PISN from Spera & Mapelli 2017, MNRAS, 470, 4739 (1).",
        * idum is the random number seed used in the kick routine. 
        *
        * Next come the parameters that determine the timesteps chosen in each
        * evolution phase:
        *                 pts1 - MS                  (0.05)
        *                 pts2 - GB, CHeB, AGB, HeGB (0.01)
        *                 pts3 - HG, HeMS            (0.02)
        * as decimal fractions of the time taken in that phase.
        """
        self.parameters.set_defaults()
        self.commit_parameters()        



Mobse = MOBSE
