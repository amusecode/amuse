from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *

from amuse.support.data.binding import InterfaceWithParametersBinding, InterfaceWithObjectsBinding

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
        #function.can_handle_array = True 
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
        #function.can_handle_array = True  
        function.addParameter('kw', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age', dtype='d', direction=function.IN)
        function.addParameter('mt', dtype='d', direction=function.IN)
        function.addParameter('tm', dtype='d', direction=function.IN)
        function.addParameter('epoch', dtype='d', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.OUT)
        return function
    
    @property
    def prototype_star(self):
        star = Particles(1)[0]
        
        star.age = 0.0 | units.s
        star.initial_mass = 1.0 | units.MSun
        star.mass = star.initial_mass
        star.type = 1 | units.no_unit
        star.luminocity = 0.0 | units.LSun
        star.radius = 0.0 | units.RSun
        star.core_mass = 0.0 | units.MSun
        star.core_radius = 0.0 | units.RSun
        star.envelope_mass = 0.0 | units.MSun
        star.envelope_radius = 0.0 | units.RSun
        star.spin = 0.0 | units.km / units.s
        star.epoch = 0.0 | units.Myr
        star.current_time = 0.0 | units.Myr
        star.main_sequence_lifetime = 0.0 | units.Myr
        #star.nuclear_burning_lifetime = 0.0 | units.Myr
        
        return star
    
    
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
        
    def initialize_star(self, star):
        star.initial_mass = star.mass
        star.age = 0.0 | units.s
        star.initial_mass = star.mass
        star.type = 1 | units.no_unit
        star.luminocity = 0.0 | units.LSun
        star.core_mass = 0.0 | units.MSun
        star.core_radius = 0.0 | units.RSun
        star.envelope_mass = 0.0 | units.MSun
        star.envelope_radius = 0.0 | units.RSun
        star.spin = 0.0 | units.km / units.s
        star.epoch = 0.0 | units.Myr
        star.current_time = 0.0 | units.Myr
        star.main_sequence_lifetime = 0.0 | units.Myr
        star.set_default("radius", 1.0 | units.RSun)
        self.evolve_star(star, 1e-06 | units.Myr)
        
        
            
    def evolve_star(self, star, target_time):
        
        current_values = {}
        current_values['kw'] = star.type.value_in(units.none)
        current_values['mass'] = star.initial_mass.value_in(units.MSun)
        current_values['mt'] = star.mass.value_in(units.MSun)
        current_values['r'] = star.radius.value_in(units.RSun)
        current_values['lum'] = star.luminocity.value_in(units.LSun)
        current_values['mc'] = star.core_mass.value_in(units.MSun)
        current_values['rc'] = star.core_radius.value_in(units.RSun)
        current_values['menv'] = star.envelope_mass.value_in(units.MSun)
        current_values['renv'] = star.envelope_radius.value_in(units.RSun)
        current_values['ospin'] = star.spin.value_in(units.km / units.s)
        current_values['epoch'] = star.epoch.value_in(units.Myr)
        current_values['tm'] = star.main_sequence_lifetime.value_in(units.Myr)
        current_values['tphys'] = star.current_time.value_in(units.Myr)
        current_values['tphysf'] = target_time.value_in(units.Myr)
        new_values = self.evolve(**current_values)
        
        time = new_values['tphysf']| units.Myr
        
        star.initial_mass =  new_values['mass'] | units.MSun
        star.mass =  new_values['mt'] | units.MSun
        star.type =  new_values['kw']  | units.no_unit
        star.luminocity =  new_values['lum'] | units.LSun
        star.radius =  new_values['r']  | units.RSun
        star.core_mass =  new_values['mc'] | units.MSun
        star.core_radius =  new_values['rc'] | units.RSun
        star.envelope_mass =  new_values['menv'] | units.MSun
        star.envelope_radius =  new_values['renv']| units.RSun
        star.spin =  new_values['ospin'] | units.km / units.s
        star.epoch =  new_values['epoch'] | units.Myr
        star.main_sequence_lifetime =  new_values['tm'] | units.Myr
        star.current_time =  time
        star.x_time =  new_values['tphys'] | units.Myr
        
    def get_time_step_for_star(self, star):
        
        current_values = {}
        current_values['kw'] = star.type.value_in(units.none)
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
    
    def evolve_particles(self, particles, time_end):
        for x in particles:
            self.evolve_particle(x, time_end)
    
    def initialize_particles(self, particles):
        for x in particles:
            self.initialize_star(x)
        
class SSEBinding(InterfaceWithParametersBinding):
    
    def __init__(self):
        InterfaceWithParametersBinding.__init__(self)
        
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
            default = 1 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "mass",
            setup_parameters = ["mt"],
            description = "current mass of a star",
            unit = units.MSun,
            default = 1 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "radius",
            setup_parameters = ["r"],
            description = "radius of a star",
            unit = units.RSun,
            default = 1 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "luminosity",
            setup_parameters = ["lum"],
            description = "luminosity of a star",
            unit = units.RSun,
            default = 1 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "core_mass",
            setup_parameters = ["mc"],
            description = "mass of the inner core of the star",
            unit = units.MSun,
            default = 1 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "core_radius",
            setup_parameters = ["rc"],
            description = "radius of the inner core of the star",
            unit = units.RSun,
            default = 1 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "envelope_mass",
            setup_parameters = ["menv"],
            description = "mass of the outer envelope core of the star",
            unit = units.MSun,
            default = 1 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "epoch",
            setup_parameters = ["epoch"],
            description = "epoch of the star",
            unit = units.Myr,
            default = 1 | units.Myr             
        ),
        attributes.AttributeDefinition(
            name = "spin",
            setup_parameters = ["ospin"],
            description = "spin of the star",
            unit = units.none,
            default = 1 | units.none             
        ),
        attributes.AttributeDefinition(
            name = "main_sequence_lifetime",
            setup_parameters = ["tm"],
            description = "main sequence lifetime of the star",
            unit = units.Myr,
            default = 1 | units.Myr             
        ),
        attributes.AttributeDefinition(
            name = "current_time",
            setup_parameters = ["tphys"],
            description = "current time of the star",
            unit = units.Myr,
            default = 1 | units.Myr             
        ),
        attributes.AttributeDefinition(
            name = "target_time",
            setup_parameters = ["tphysf"],
            description = "target time of the star",
            unit = units.Myr,
            default = 1 | units.Myr             
        ),    
    ]
    
class SSE(SSEInterface, SSEBinding):
    
    def __init__(self):
        SSEInterface.__init__(self)
        SSEBinding.__init__(self)
        
        
