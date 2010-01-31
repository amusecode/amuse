from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *
from amuse.legacy.support.lit import LiteratureRefs

from amuse.support.data.values import Quantity
from amuse.support.data.binding import InterfaceWithParametersBinding, InterfaceWithObjectsBinding
from amuse.support.data.binding import ParticleAttributesModifier
from amuse.support.data.core import Particles, ParticlesSubset

class BSEInterface(LegacyInterface, LiteratureRefs): 
    def __init__(self):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")
        LiteratureRefs.__init__(self)

    @legacy_function   
    def initialize():
        function = LegacyFunctionSpecification()  
        function.addParameter('z_in', dtype='d', direction=function.IN)
        function.addParameter('neta_in', dtype='d', direction=function.IN)
        function.addParameter('bwind_in', dtype='d', direction=function.IN)
        function.addParameter('hewind_in', dtype='d', direction=function.IN)
        function.addParameter('alpha1_in', dtype='d', direction=function.IN)
        function.addParameter('CElambda_in', dtype='d', direction=function.IN)
        function.addParameter('ceflag_in', dtype='i', direction=function.IN)
        function.addParameter('tflag_in', dtype='i', direction=function.IN)
        function.addParameter('ifflag_in', dtype='i', direction=function.IN)
        function.addParameter('wdflag_in', dtype='i', direction=function.IN)
        function.addParameter('bhflag_in', dtype='i', direction=function.IN)
        function.addParameter('nsflag_in', dtype='i', direction=function.IN)
        function.addParameter('mxns_in', dtype='i', direction=function.IN)
        function.addParameter('pts1_in', dtype='d', direction=function.IN)
        function.addParameter('pts2_in', dtype='d', direction=function.IN)
        function.addParameter('pts3_in', dtype='d', direction=function.IN)
        function.addParameter('sigma_in', dtype='d', direction=function.IN)
        function.addParameter('beta_in', dtype='d', direction=function.IN)
        function.addParameter('xi_in', dtype='d', direction=function.IN)
        function.addParameter('acc2_in', dtype='d', direction=function.IN)
        function.addParameter('epsnov_in', dtype='d', direction=function.IN)
        function.addParameter('eddfac_in', dtype='d', direction=function.IN)
        function.addParameter('gamma_in', dtype='d', direction=function.IN)
        function.addParameter('status', dtype='i', direction=function.OUT)
        return function
        
    @legacy_function     
    def evolve():
        function = LegacyFunctionSpecification()  
        function.name = 'evolve_binary'
        function.can_handle_array = True 
        function.addParameter('type1', dtype='i', direction=function.INOUT)
        function.addParameter('type2', dtype='i', direction=function.INOUT)
        function.addParameter('initial_mass1', dtype='d', direction=function.INOUT)
        function.addParameter('initial_mass2', dtype='d', direction=function.INOUT)
        function.addParameter('mass1', dtype='d', direction=function.INOUT)
        function.addParameter('mass2', dtype='d', direction=function.INOUT)
        function.addParameter('radius1', dtype='d', direction=function.INOUT)
        function.addParameter('radius2', dtype='d', direction=function.INOUT)
        function.addParameter('luminosity1', dtype='d', direction=function.INOUT)
        function.addParameter('luminosity2', dtype='d', direction=function.INOUT)
        function.addParameter('core_mass1', dtype='d', direction=function.INOUT)
        function.addParameter('core_mass2', dtype='d', direction=function.INOUT)
        function.addParameter('core_radius1', dtype='d', direction=function.INOUT)
        function.addParameter('core_radius2', dtype='d', direction=function.INOUT)
        function.addParameter('envelope_mass1', dtype='d', direction=function.INOUT)
        function.addParameter('envelope_mass2', dtype='d', direction=function.INOUT)
        function.addParameter('envelope_radius1', dtype='d', direction=function.INOUT)
        function.addParameter('envelope_radius2', dtype='d', direction=function.INOUT)
        function.addParameter('spin1', dtype='d', direction=function.INOUT)
        function.addParameter('spin2', dtype='d', direction=function.INOUT)
        function.addParameter('epoch1', dtype='d', direction=function.INOUT)
        function.addParameter('epoch2', dtype='d', direction=function.INOUT)
        function.addParameter('MS_lifetime1', dtype='d', direction=function.INOUT)
        function.addParameter('MS_lifetime2', dtype='d', direction=function.INOUT)
        function.addParameter('age', dtype='d', direction=function.INOUT)
        function.addParameter('end_time', dtype='d', direction=function.INOUT)
        function.addParameter('orbital_period', dtype='d', direction=function.INOUT)
        function.addParameter('eccentricity', dtype='d', direction=function.INOUT)
        return function
        
    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('type1', dtype='i', direction=function.IN)
        function.addParameter('type2', dtype='i', direction=function.IN)
        function.addParameter('initial_mass1', dtype='d', direction=function.IN)
        function.addParameter('initial_mass2', dtype='d', direction=function.IN)
        function.addParameter('mass1', dtype='d', direction=function.IN)
        function.addParameter('mass2', dtype='d', direction=function.IN)
        function.addParameter('MS_lifetime1', dtype='d', direction=function.IN)
        function.addParameter('MS_lifetime2', dtype='d', direction=function.IN)
        function.addParameter('epoch1', dtype='d', direction=function.IN)
        function.addParameter('epoch2', dtype='d', direction=function.IN)
        function.addParameter('age', dtype='d', direction=function.IN)
        function.addParameter('time_step', dtype='d', direction=function.OUT)
        return function
        
    def get_time_step_for_star(self, star):
        
        current_values = {}
        current_values['kw'] = star.type.value_in(units.stellar_type)
        current_values['mass'] = star.initial_mass.value_in(units.MSun)
        current_values['mt'] = star.mass.value_in(units.MSun)
        current_values['tm'] = star.main_sequence_lifetime.value_in(units.Myr)
        current_values['age'] = star.age.value_in(units.Myr)
        current_values['epoch'] = star.epoch.value_in(units.Myr)
        
        result = self.get_time_step(**current_values)
        
        return result | units.Myr
        
        
    def evolve_particle(self, particle, time_end):
        t = particle.current_time
        if particle.type.value_in(units.none) == 15:
            return
        while t < time_end:
            t0 = t
            t  = t0 + self.get_time_step_for_star(particle)
            if t > time_end:
                t = time_end
            self.evolve_star(particle, t)
            t1 = particle.current_time
            dt = t1 - t0
            t0 = t1
            if dt.value_in(units.Myr) == 0.0:
                print t, t0, t1, dt, "BREAK BREAK BREAK!"
                return
            if particle.type.value_in(units.none) == 15:
                return
    
        
class BSEParticles(Particles):
    
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
            "type" : 1 | units.stellar_type,
            "luminosity":  0 | units.LSun,
            "core_mass": 0 | units.MSun,
            "core_radius": 0 | units.RSun,
            "envelope_mass": 0 | units.MSun,
            "envelope_radius": 0 | units.RSun,
            "epoch": 0 | units.Myr ,
            "spin": 0 | units.none ,
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
        
        super(BSEParticles, self)._set_particles(keys, all_attributes, all_values)
        
        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface.evolve_particles_method._run(self._private.code_interface, added_particles, 1e-06 | units.Myr)
    
    def _state_attributes(self):
        return ["mass", "radius"]
        
class BSEBinding(InterfaceWithParametersBinding):
    
    def __init__(self):
        InterfaceWithParametersBinding.__init__(self)
        self.particles = BSEParticles(self)
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
            0.5 | units.none
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
            1.0/8.0 | units.none
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
            3.0/2.0 | units.none
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
            ("type", "kw", units.stellar_type),
            ("initial_mass", "mass", units.MSun),
            ("mass", "mt", units.MSun),
            ("main_sequence_lifetime", "tm", units.Myr),
            ("age", "age", units.Myr),
            ("epoch", "epoch", units.Myr),
            ("time_step", "dt", units.Myr),
        )
    )
        
         
    evolve_particles_method = ParticleAttributesModifier(
        "evolve",
        (
            ("type", "kw", units.stellar_type),
            ("initial_mass", "mass", units.MSun),
            ("mass", "mt", units.MSun),
            ("radius", "r", units.RSun),
            ("luminosity", "lum", units.LSun),
            ("core_mass", "mc", units.MSun),
            ("core_radius", "rc", units.RSun),
            ("envelope_mass", "menv", units.MSun),
            ("envelope_radius", "renv", units.RSun),
            ("epoch", "epoch", units.Myr),
            ("spin", "ospin", units.none),
            ("main_sequence_lifetime", "tm", units.Myr),
            ("age", "tphys", units.Myr),
        ),
        (
            ("end_time", 'tphysf', units.Myr),
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
            self.parameters.SN_kick_speed_dispersion.value_in(units.km / units.s),
            self.parameters.white_dwarf_IFMR_flag.value_in(units.none), 
            self.parameters.white_dwarf_cooling_flag.value_in(units.none), 
            self.parameters.black_hole_kick_flag.value_in(units.none), 
            self.parameters.neutron_star_mass_flag.value_in(units.none), 
            self.parameters.maximum_neutron_star_mass.value_in(units.MSun),
            self.parameters.fractional_time_step_1.value_in(units.none), 
            self.parameters.fractional_time_step_2.value_in(units.none), 
            self.parameters.fractional_time_step_3.value_in(units.none))
        
    def initialize_module_with_default_parameters(self):
        """
        * neta is the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
        * bwind is the binary enhanced mass loss parameter (inactive for single).
        * hewind is a helium star mass loss factor (1.0 normally). 
        * sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
        *
        * ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
        * wdflag > 0 uses modified-Mestel cooling for WDs (0). 
        * bhflag > 0 allows velocity kick at BH formation (0). 
        * nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
        * mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
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
        self.initialize_module_with_current_parameters()
    
    
class BSE(BSEInterface, BSEBinding):
    
        
    def __init__(self):
        BSEInterface.__init__(self)
        BSEBinding.__init__(self)
    
    
        
        
        
         
        
        
