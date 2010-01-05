from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *

from amuse.support.data.values import Quantity
from amuse.support.data.binding import InterfaceWithParametersBinding, InterfaceWithObjectsBinding
from amuse.support.data.binding import ParticleAttributesModifier
from amuse.support.data.core import Particles

class SSEInterface(LegacyInterface): 
    def __init__(self):
        LegacyInterface.__init__(self)

    
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
        function.addParameter('mxns_in', dtype='i', direction=function.IN)
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
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = self.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
            
   
        
    def get_time_step_for_star(self, star):
        
        current_values = {}
        current_values['kw'] = star.type.value_in(units.stellar_type)
        current_values['mass'] = star.initial_mass.value_in(units.MSun)
        current_values['mt'] = star.mass.value_in(units.MSun)
        current_values['tm'] = star.main_sequence_lifetime.value_in(units.Myr)
        current_values['age'] = star.current_time.value_in(units.Myr)
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
    
        
class SSEParticles(Particles):
    
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
            "current_time": 0 | units.Myr
        }
        
        given_attributes = set(attributes)
        
        if not "initial_mass" in given_attributes:
            index_of_mass_attibute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attibute])
        
        for attribute, default_value in mapping_from_attribute_to_default_value.iteritems():
            if not attribute in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))
        
        super(SSEParticles, self)._set_particles(keys, all_attributes, all_values)
        
        self._private.code_interface.evolve_model(1e-06 | units.Myr)
    
    def _get_attributes_old(self):
        return ["mass", "radius"]
        
class SSEBinding(InterfaceWithParametersBinding):
    
    def __init__(self):
        InterfaceWithParametersBinding.__init__(self)
        self.particles = SSEParticles(self)
        
    update_time_step_method = ParticleAttributesModifier(
        "get_time_step",
        (
            ("type", "kw", units.stellar_type),
            ("initial_mass", "mass", units.MSun),
            ("mass", "mt", units.MSun),
            ("main_sequence_lifetime", "tm", units.Myr),
            ("current_time", "age", units.Myr),
            ("epoch", "epoch", units.Myr),
            ("dt", "timestep", units.Myr),
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
            ("current_time", "tphys", units.Myr),
        ),
        (
            ("end_time", 'tphysf', units.Myr),
        )
    )
    
    attribute_definitions = [
        attributes.AttributeDefinition(
            name = "type",
            setup_parameters = ["kw"],
            description = "star type",
            unit = units.stellar_type,
            default = 1 | units.stellar_type             
        ),
        attributes.AttributeDefinition(
            name = "initial_mass",
            setup_parameters = ["mass"],
            description = "initial ZAMS mass of a star",
            unit = units.MSun,
            default = None            
        ),
        attributes.AttributeDefinition(
            name = "mass",
            setup_parameters = ["mt"],
            description = "current mass of a star",
            unit = units.MSun,
            default = None         
        ),
        attributes.AttributeDefinition(
            name = "radius",
            setup_parameters = ["r"],
            description = "radius of a star",
            unit = units.RSun,
            default = None   
        ),
        attributes.AttributeDefinition(
            name = "luminosity",
            setup_parameters = ["lum"],
            description = "luminosity of a star",
            unit = units.LSun,
            default = 0 | units.LSun             
        ),
        attributes.AttributeDefinition(
            name = "core_mass",
            setup_parameters = ["mc"],
            description = "mass of the inner core of the star",
            unit = units.MSun,
            default = 0 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "core_radius",
            setup_parameters = ["rc"],
            description = "radius of the inner core of the star",
            unit = units.RSun,
            default = 0 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "envelope_mass",
            setup_parameters = ["menv"],
            description = "mass of the outer envelope of the star",
            unit = units.MSun,
            default = 0 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "envelope_radius",
            setup_parameters = ["renv"],
            description = "radius of the oute envelope of the star",
            unit = units.RSun,
            default = 0 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "epoch",
            setup_parameters = ["epoch"],
            description = "epoch of the star",
            unit = units.Myr,
            default = 0 | units.Myr             
        ),
        attributes.AttributeDefinition(
            name = "spin",
            setup_parameters = ["ospin"],
            description = "spin of the star",
            unit = units.none,
            default = 0 | units.none             
        ),
        attributes.AttributeDefinition(
            name = "main_sequence_lifetime",
            setup_parameters = ["tm"],
            description = "main sequence lifetime of the star",
            unit = units.Myr,
            default = 0 | units.Myr             
        ),
        attributes.AttributeDefinition(
            name = "current_time",
            setup_parameters = ["tphys"],
            description = "current time of the star",
            unit = units.Myr,
            default = 0 | units.Myr             
        ),
    ]
    
    
    def setup_particles(self, particles):
        self.particles = particles
        
        attributes = []
        all_units = []
        values = []
        
        for attribute_definition in self.attribute_definitions:
            if not attribute_definition.default_value is None:
                attribute_definition.for_default_fill_arguments_for_attributelist_set(
                    attributes,
                    all_units,
                    values,
                    len(particles)                    
                )
         
        list_of_values = particles._private.attribute_storage._get_values(None, ["mass"])
        
        
        attributes.append("initial_mass")
        all_units.append(units.MSun)
        values.append(list_of_values[0].value_in(units.MSun))
        
        particles._private.attribute_storage.set_values_of_particles_in_units(None, attributes, values, all_units)
        self.evolve_particles(particles, 1e-06 | units.Myr)
        
        
    def evolve_model(self, end_time):
        self.evolve_particles_method._run(self, self.particles, end_time)
        
    def update_time_steps(self):
        self.update_time_step_method._run(self, self.particles)
        
    def evolve_model_using_timesteps(self, time_end):
        self.update_time_steps()
        next_time = self.particles.time_step + self.particles.current_time
        self.evolve_model(next_time)
        
    
    def evolve_particles(self, particles, target_times):
       
        attributes = []
        all_units = []
        values = []
        keywords = []
        
        for attribute_definition in self.attribute_definitions:
            attribute_definition.for_setup_fill_arguments_for_attributelist_get(
                attributes,
                all_units,
                keywords,
            )
    
        list_of_values = particles._private.attribute_storage.get_values_of_particles_in_units(None, attributes, all_units)
        
        keyword_arguments = {}
        for keyword, values in zip(keywords, list_of_values):
            keyword_arguments[keyword] = values
        
        
        if isinstance(target_times, Quantity):
            keyword_arguments['tphysf'] = [target_times.value_in(units.Myr)] * len(particles)
        else:
            keyword_arguments['tphysf'] = [x.value_in(units.Myr) for x in target_times]
        
        
        results = self.evolve(**keyword_arguments)
        
        attributes = []
        all_units = []
        values = []
        
        
        for attribute_definition in self.attribute_definitions:
            attribute_definition.for_state_fill_arguments_for_attributelist_set(
                results,
                attributes,
                all_units,
                values,
            )
        
        time = None
        
        particles._private.attribute_storage.set_values_of_particles_in_units(None, attributes, values, all_units, time)
        
    def get_timesteps(self, particles):
        attribute_names = ["type","initial_mass","mass","main_sequence_lifetime", "current_time", "epoch"]
        keywords = ["kw", "mass", "mt", "tm", "age", "epoch"]
        all_units = [units.stellar_type, units.MSun, units.MSun, units.Myr, units.Myr, units.Myr]
                
        list_of_values = particles._private.attribute_storage.get_values_of_particles_in_units(None, attribute_names, all_units)
        
        keyword_arguments = {}
        for keyword, values in zip(keywords, list_of_values):
            keyword_arguments[keyword] = values
        
        
        results = self.get_time_step(**keyword_arguments)
        return units.Myr.new_quantity(results)
        
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
    
    
class SSE(SSEInterface, SSEBinding):
    
        
    def __init__(self):
        SSEInterface.__init__(self)
        SSEBinding.__init__(self)
    
    
        
        
        
         
        
        
