from amuse.support.units import units
from amuse.support.data.core import Particle

from amuse.legacy import *

class SSE(LegacyInterface): 
    def __init__(self):
        LegacyInterface.__init__(self)
        
    @legacy_function   
    def initialize():
        function = RemoteFunction()  
        function.id = 1
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
        function = RemoteFunction()  
        function.id = 2
        function.name = 'evolve0'
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
        function = RemoteFunction()      
        function.id = 3
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
        star = Particle(-1)
        
        star.age = 0.0 | units.s
        star.initial_mass = 1.0 | units.MSun
        star.mass = star.initial_mass.value()
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
        star.initial_mass = star.mass.value()
        star.ensure_attributes_like(self.prototype_star)
        star.x_time = 0 | units.Myr
        self.evolve_star(star, 1e-06 | units.Myr)
            
    def evolve_star(self, star, target_time):
        
        current_values = {}
        current_values['kw'] = star.type.value().number
        current_values['mass'] = star.initial_mass.to_number_in(units.MSun)
        current_values['mt'] = star.mass.to_number_in(units.MSun)
        current_values['r'] = star.radius.to_number_in(units.RSun)
        current_values['lum'] = star.luminocity.to_number_in(units.LSun)
        current_values['mc'] = star.core_mass.to_number_in(units.MSun)
        current_values['rc'] = star.core_radius.to_number_in(units.RSun)
        current_values['menv'] = star.envelope_mass.to_number_in(units.MSun)
        current_values['renv'] = star.envelope_radius.to_number_in(units.RSun)
        current_values['ospin'] = star.spin.to_number_in(units.km / units.s)
        current_values['epoch'] = star.epoch.to_number_in(units.Myr)
        current_values['tm'] = star.main_sequence_lifetime.to_number_in(units.Myr)
        current_values['tphys'] = star.current_time.to_number_in(units.Myr)
        current_values['tphysf'] = target_time.in_(units.Myr).number
        new_values = self.evolve(**current_values)
        
        time = new_values['tphysf']| units.Myr
        
        
        star.initial_mass.set_value_at_time(time, new_values['mass'] | units.MSun)
        star.mass.set_value_at_time(time, new_values['mt'] | units.MSun)
        star.type.set_value_at_time(time, new_values['kw']  | units.no_unit)
        star.luminocity.set_value_at_time(time, new_values['lum'] | units.LSun)
        star.radius.set_value_at_time(time, new_values['r']  | units.RSun)
        star.core_mass.set_value_at_time(time, new_values['mc'] | units.MSun)
        star.core_radius.set_value_at_time(time, new_values['rc'] | units.RSun)
        star.envelope_mass.set_value_at_time(time, new_values['menv'] | units.MSun)
        star.envelope_radius.set_value_at_time(time, new_values['renv']| units.RSun)
        star.spin.set_value_at_time(time, new_values['ospin'] | units.km / units.s)
        star.epoch.set_value_at_time(time, new_values['epoch'] | units.Myr)
        star.main_sequence_lifetime.set_value_at_time(time, new_values['tm'] | units.Myr)
        star.current_time.set_value_at_time(time, time)
        star.x_time.set_value_at_time(time, new_values['tphys'] | units.Myr)
        
    def get_time_step_for_star(self, star):
        
        current_values = {}
        current_values['kw'] = star.type.value().number
        current_values['mass'] = star.initial_mass.to_number_in(units.MSun)
        current_values['mt'] = star.mass.to_number_in(units.MSun)
        current_values['tm'] = star.main_sequence_lifetime.to_number_in(units.Myr)
        current_values['age'] = star.current_time.to_number_in(units.Myr)
        current_values['epoch'] = star.epoch.to_number_in(units.Myr)
        
        result = self.get_time_step(**current_values)
        
        return result | units.Myr
        
        
    def evolve_particle(self, particle, time_end):
        t = particle.current_time.value()
        if particle.type.value().number == 15:
            return
        while t < time_end:
            t0 = t
            t  = t0 + self.get_time_step_for_star(particle)
            self.evolve_star(particle, t)
            t1 = particle.current_time.value()
            dt = t1 - t0
            t0 = t1
            if dt.in_(units.Myr).number == 0.0:
                print t, t0, t1, dt, "BREAK BREAK BREAK!"
                return
            if particle.type.value().number == 15:
                return
                
    def initialize_particles(self, particles):
        for x in particles:
            self.initialize_star(x)
        
        
