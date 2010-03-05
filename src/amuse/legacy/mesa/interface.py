import os

from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *
from amuse.legacy.support.lit import LiteratureRefs
from amuse.legacy.interface.se import StellarEvolution

from amuse.support.data.values import Quantity
from amuse.support.data.binding import ParticleAttributesModifier
from amuse.support.data.core import Particles, ParticlesSubset
from amuse.support.interface import CodeInterfaceOld

class MESAInterface(LegacyInterface, LiteratureRefs, StellarEvolution): 
    """
        .. [#]  MESA is free software; you can redistribute it and/or modify
                it under the terms of the GNU General Library Public License.
    """
    def __init__(self):
        try:
            LegacyInterface.__init__(self, name_of_the_worker="worker_code")
            LiteratureRefs.__init__(self)
            self.MESA_exists = True
        except Exception:
            print "MESA was not built. Skipping initialization."
            self.MESA_exists = False

    @property
    def default_path_to_inlist(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'mesa_reqs', 'AMUSE_inlist')

    @property
    def default_path_to_MESA_data(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'src', 'data')

    @legacy_function   
    def initialize():
        function = LegacyFunctionSpecification()  
        #function.addParameter('gamma_in', dtype='d', direction=function.IN)
        function.addParameter('inlist_path', dtype='string', direction=function.IN,
            description = "Path to the inlist file.")
        function.addParameter('MESA_data_path', dtype='string', direction=function.IN,
            description = "Path to the MESA data directory.")
        function.addParameter('status', dtype='i', direction=function.OUT)
        return function
        
    @legacy_function     
    def evolve_to():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('end_time', dtype='d', direction=function.IN)
        return function
        
    @legacy_function     
    def new_zams_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('status', dtype='i', direction=function.OUT)
        return function
        
#~     @legacy_function      
#~     def get_time_step():
#~         function = LegacyFunctionSpecification() 
#~         function.can_handle_array = True
#~         function.addParameter('type1', dtype='i', direction=function.IN)
#~         function.addParameter('type2', dtype='i', direction=function.IN)
#~         function.addParameter('initial_mass1', dtype='d', direction=function.IN)
#~         function.addParameter('initial_mass2', dtype='d', direction=function.IN)
#~         function.addParameter('mass1', dtype='d', direction=function.IN)
#~         function.addParameter('mass2', dtype='d', direction=function.IN)
#~         function.addParameter('MS_lifetime1', dtype='d', direction=function.IN)
#~         function.addParameter('MS_lifetime2', dtype='d', direction=function.IN)
#~         function.addParameter('epoch1', dtype='d', direction=function.IN)
#~         function.addParameter('epoch2', dtype='d', direction=function.IN)
#~         function.addParameter('age', dtype='d', direction=function.IN)
#~         function.addParameter('time_step', dtype='d', direction=function.OUT)
#~         return function
        

class MESAParticles(Particles):
    
    def __init__(self, code_interface):
        Particles.__init__(self)
        self._private.code_interface = code_interface 
        
    def _set_particles(self, keys, attributes = [], values = []):
        if len(keys) == 0:
            return
            
        all_attributes = []
        all_attributes.extend(attributes)
        all_values = []
        all_values.extend(values)
        
        mapping_from_attribute_to_default_value = {
            "type1" : 1 | units.stellar_type,
            "type2" : 1 | units.stellar_type,
            "radius1":  0 | units.RSun,
            "radius2":  0 | units.RSun,
            "luminosity1":  0 | units.LSun,
            "luminosity2":  0 | units.LSun,
            "core_mass1": 0 | units.MSun,
            "core_mass2": 0 | units.MSun,
            "core_radius1": 0 | units.RSun,
            "core_radius2": 0 | units.RSun,
            "envelope_mass1": 0 | units.MSun,
            "envelope_mass2": 0 | units.MSun,
            "envelope_radius1": 0 | units.RSun,
            "envelope_radius2": 0 | units.RSun,
            "epoch1": 0 | units.Myr,
            "epoch2": 0 | units.Myr,
            "spin1": 0 | units.none,
            "spin2": 0 | units.none,
            "MS_lifetime1": 0 | units.Myr,
            "MS_lifetime2": 0 | units.Myr,
#            "orbital_period": 200.0 | units.day,
            "eccentricity": 0.0 | units.none,
            "age": 0 | units.Myr
        }
        
        given_attributes = set(attributes)
        
        if not "initial_mass1" in given_attributes:
            index_of_mass_attibute = attributes.index("mass1")
            all_attributes.append("initial_mass1")
            all_values.append(values[index_of_mass_attibute] * 1.0)
        if not "initial_mass2" in given_attributes:
            index_of_mass_attibute = attributes.index("mass2")
            all_attributes.append("initial_mass2")
            all_values.append(values[index_of_mass_attibute] * 1.0)
        
        for attribute, default_value in mapping_from_attribute_to_default_value.iteritems():
            if not attribute in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))
        
        super(MESAParticles, self)._set_particles(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface.evolve_particles_method._run(self._private.code_interface, added_particles, 1e-06 | units.Myr)
    
    def _state_attributes(self):
        return ["mass1", "mass2", "radius1", "radius2"]
        
class MESABinding(CodeInterfaceOld):
    
    def __init__(self):
        CodeInterfaceOld.__init__(self)
        self.particles = MESAParticles(self)
        self.parameters.set_defaults()
   
    parameter_definitions = [
        parameters.ModuleCachingParameterDefinition(
            "metallicity",
            "metallicity",
            "Metallicity of all stars",
            units.none,
            0.02 | units.none
        ),
                
        parameters.ModuleCachingParameterDefinition(
            "neta",
            "reimers_mass_loss_coefficient",
            "Reimers mass-loss coefficient (neta*4x10^-13; 0.5 normally)",
            units.none,
            0.5 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "bwind",
            "binary_enhanced_mass_loss_parameter",
            "The binary enhanced mass loss parameter (inactive for single).",
            units.none,
            0.0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "hewind",
            "helium_star_mass_loss_factor",
            "Helium star mass loss factor",
            units.none,
            1.0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "alpha1",
            "common_envelope_efficiency",
            "The common-envelope efficiency parameter",
            units.none,
            1.0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "CElambda",
            "common_envelope_binding_energy_factor",
            "The binding energy factor for common envelope evolution",
            units.none,
            0.5 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "ceflag",
            "common_envelope_model_flag", 
            "ceflag > 0 activates spin-energy correction in common-envelope. ceflag = 3 activates de Kool common-envelope model (0).",
            units.none, 
            0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "tflag",
            "tidal_circularisation_flag", 
            "tflag > 0 activates tidal circularisation (1).",
            units.none, 
            1 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "ifflag",
            "white_dwarf_IFMR_flag", 
            "ifflag > 0 uses white dwarf IFMR (initial-final mass relation) of HPE, 1995, MNRAS, 272, 800 (0).",
            units.none, 
            0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "wdflag",
            "white_dwarf_cooling_flag", 
            "wdflag > 0 uses modified-Mestel cooling for WDs (0).",
            units.none, 
            1 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "bhflag",
            "black_hole_kick_flag",
            "bhflag > 0 allows velocity kick at BH formation (0).",
            units.none,
            0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "nsflag",
            "neutron_star_mass_flag",
            "nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1).",
            units.none,
            1 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "mxns",
            "maximum_neutron_star_mass",
            "The maximum neutron star mass (1.8, nsflag=0; 3.0, nsflag=1).",
            units.MSun,
            3.0 | units.MSun
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "idum",
            "SN_kick_random_seed", 
            "The random number seed used in the kick routine.",
            units.none, 
            29769 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "pts1",
            "fractional_time_step_1", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: MS (0.05)",
            units.none, 
            0.05 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "pts2",
            "fractional_time_step_2", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: GB, CHeB, AGB, HeGB (0.01)",
            units.none, 
            0.01 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "pts3",
            "fractional_time_step_3", 
            "The timesteps chosen in each evolution phase as decimal fractions of the time taken in that phase: HG, HeMS (0.02)",
            units.none, 
            0.02 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "sigma",
            "SN_kick_speed_dispersion",
            "The dispersion in the Maxwellian for the SN kick speed (190 km/s).",
            units.km / units.s,
            190.0 | units.km / units.s
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "beta",
            "wind_velocity_factor",
            "The wind velocity factor: proportional to vwind**2 (1/8).",
            units.none,
            0.125 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "xi",
            "wind_accretion_efficiency",
            "The wind accretion efficiency factor (1.0).",
            units.none,
            1.0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "acc2",
            "wind_accretion_factor",
            "The Bondi-Hoyle wind accretion factor (3/2).",
            units.none,
            1.5 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "epsnov",
            "nova_retained_accreted_matter_fraction",
            "The fraction of accreted matter retained in nova eruption (0.001).",
            units.none,
            0.001 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "eddfac",
            "Eddington_mass_transfer_limit_factor",
            "The Eddington limit factor for mass transfer (1.0).",
            units.none,
            1.0 | units.none
        ),
        
        parameters.ModuleCachingParameterDefinition(
            "gamma",
            "Roche_angular_momentum_factor",
            "The angular momentum factor for mass lost during Roche (-1.0). ",
            units.none,
            -1.0 | units.none
        ),
        
    ]
    
    update_time_step_method = ParticleAttributesModifier(
        "get_time_step",
        (
            ("type1", "type1", units.stellar_type),
            ("type2", "type2", units.stellar_type),
            ("initial_mass1", "initial_mass1", units.MSun),
            ("initial_mass2", "initial_mass2", units.MSun),
            ("mass1", "mass1", units.MSun),
            ("mass2", "mass2", units.MSun),
            ("MS_lifetime1", "MS_lifetime1", units.Myr),
            ("MS_lifetime2", "MS_lifetime2", units.Myr),
            ("epoch1", "epoch1", units.Myr),
            ("epoch2", "epoch2", units.Myr),
            ("age", "age", units.Myr),
            ("time_step", "time_step", units.Myr),
        )
    )
    
    evolve_particles_method = ParticleAttributesModifier(
        "evolve",
        (
            ("type1", "type1", units.stellar_type),
            ("type2", "type2", units.stellar_type),
            ("initial_mass1", "initial_mass1", units.MSun),
            ("initial_mass2", "initial_mass2", units.MSun),
            ("mass1", "mass1", units.MSun),
            ("mass2", "mass2", units.MSun),
            ("radius1", "radius1", units.RSun),
            ("radius2", "radius2", units.RSun),
            ("luminosity1", "luminosity1", units.LSun),
            ("luminosity2", "luminosity2", units.LSun),
            ("core_mass1", "core_mass1", units.MSun),
            ("core_mass2", "core_mass2", units.MSun),
            ("core_radius1", "core_radius1", units.RSun),
            ("core_radius2", "core_radius2", units.RSun),
            ("envelope_mass1", "envelope_mass1", units.MSun),
            ("envelope_mass2", "envelope_mass2", units.MSun),
            ("envelope_radius1", "envelope_radius1", units.RSun),
            ("envelope_radius2", "envelope_radius2", units.RSun),
            ("spin1", "spin1", units.none),
            ("spin2", "spin2", units.none),
            ("epoch1", "epoch1", units.Myr),
            ("epoch2", "epoch2", units.Myr),
            ("MS_lifetime1", "MS_lifetime1", units.Myr),
            ("MS_lifetime2", "MS_lifetime2", units.Myr),
            ("age", "age", units.Myr),
            ("orbital_period", "orbital_period", units.day),
            ("eccentricity", "eccentricity", units.none),
        ),
        (
            ("end_time", 'end_time', units.Myr),
        )
    )
    
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def evolve_model(self, end_time = None):
        if end_time is None:
            self.update_time_steps()
            end_time = self.particles.time_step + self.particles.age
        self.evolve_particles_method._run(self, self.particles, end_time)
        
    def update_time_steps(self):
        self.update_time_step_method._run(self, self.particles)
    
    # Currently unused:
    def evolve_model_using_timesteps(self, end_time = None):
        
        particles_set = particles.to_set()
        while len(particles_set) > 0:
            self._evolve_particles(particles_set)
            particles_set = particles_set.select(lambda x : x < end_time, ["age"])
            print len(particles_set)
            print particles_set.age
            
                
    def _evolve_particles(self, particles):
        self.update_time_step_method._run(self, particles)
        end_times = particles.time_step + particles.age
        self.evolve_particles_method._run(self, particles, end_times)
    
    def initialize_stars(self):
        pass
    
    def initialize_module_with_current_parameters(self):
        status = self.initialize(self.parameters.metallicity.value_in(units.none),
            self.parameters.reimers_mass_loss_coefficient.value_in(units.none), 
            self.parameters.binary_enhanced_mass_loss_parameter.value_in(units.none), 
            self.parameters.helium_star_mass_loss_factor.value_in(units.none), 
            self.parameters.common_envelope_efficiency.value_in(units.none), 
            self.parameters.common_envelope_binding_energy_factor.value_in(units.none), 
            self.parameters.common_envelope_model_flag.value_in(units.none), 
            self.parameters.tidal_circularisation_flag.value_in(units.none), 
            self.parameters.white_dwarf_IFMR_flag.value_in(units.none), 
            self.parameters.white_dwarf_cooling_flag.value_in(units.none), 
            self.parameters.black_hole_kick_flag.value_in(units.none), 
            self.parameters.neutron_star_mass_flag.value_in(units.none), 
            self.parameters.maximum_neutron_star_mass.value_in(units.MSun),
            self.parameters.SN_kick_random_seed.value_in(units.none), 
            self.parameters.fractional_time_step_1.value_in(units.none), 
            self.parameters.fractional_time_step_2.value_in(units.none), 
            self.parameters.fractional_time_step_3.value_in(units.none),
            self.parameters.SN_kick_speed_dispersion.value_in(units.km / units.s),
            self.parameters.wind_velocity_factor.value_in(units.none), 
            self.parameters.wind_accretion_efficiency.value_in(units.none), 
            self.parameters.wind_accretion_factor.value_in(units.none), 
            self.parameters.nova_retained_accreted_matter_fraction.value_in(units.none), 
            self.parameters.Eddington_mass_transfer_limit_factor.value_in(units.none), 
            self.parameters.Roche_angular_momentum_factor.value_in(units.none))
        
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.initialize(self.default_path_to_inlist, self.default_path_to_MESA_data)
    
    
class MESA(MESAInterface, MESABinding):
    def __init__(self):
        MESAInterface.__init__(self)
        MESABinding.__init__(self)
    
    
        
        
        
         
        
        
