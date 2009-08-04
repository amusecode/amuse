import os.path
from mpi4py import MPI
import numpy

from amuse.legacy.support import core

from amuse.legacy.support.core import RemoteFunction, legacy_global
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Star

class SSE(object): 
    def __init__(self):
        directory_of_this_module = os.path.dirname(__file__);
        full_name_of_the_worker = os.path.join(directory_of_this_module , 'muse_worker1')
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, 1)
        self.channel = core.MpiChannel(self.intercomm)
    
    def __del__(self):
        self.stop_worker()
        
    @core.legacy_function
    def stop_worker():
        function = RemoteFunction()  
        function.id = 0
        return function;
        
    @core.legacy_function   
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
        
    @core.legacy_function     
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
        
    @core.legacy_function      
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
        star = Star(-1)
        
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
        star.ensure_attributes_like(self.prototype_star)
        self.evolve_star(star, 1e-06 | units.Myr)
            
    def evolve_star(self, star, target_time):
        
        current_values = {}
        current_values['kw'] = star.type.number
        current_values['mass'] = star.initial_mass.in_(units.MSun).number
        current_values['mt'] = star.mass.in_(units.MSun).number
        current_values['r'] = star.radius.in_(units.RSun).number
        current_values['lum'] = star.luminocity.in_(units.LSun).number
        current_values['mc'] = star.core_mass.in_(units.MSun).number
        current_values['rc'] = star.core_radius.in_(units.RSun).number
        current_values['menv'] = star.envelope_mass.in_(units.MSun).number
        current_values['renv'] = star.envelope_radius.in_(units.RSun).number
        current_values['ospin'] = star.spin.in_(units.km / units.s).number
        current_values['epoch'] = star.epoch.in_(units.Myr).number
        current_values['tm'] = star.main_sequence_lifetime.in_(units.Myr).number
        current_values['tphys'] = star.current_time.in_(units.Myr).number
        current_values['tphysf'] = target_time.in_(units.Myr).number
        
        new_values = self.evolve(**current_values)
        
        star.age = 0.0 | units.s
        star.initial_mass = new_values['mass'] | units.MSun
        star.mass = new_values['mt'] | units.MSun
        star.type = new_values['kw']  | units.no_unit
        star.luminocity = new_values['lum'] | units.LSun
        star.radius = new_values['r']  | units.RSun
        star.core_mass = new_values['mc'] | units.MSun
        star.core_radius = new_values['rc'] | units.RSun
        star.envelope_mass = new_values['menv'] | units.MSun
        star.envelope_radius = new_values['renv']| units.RSun
        star.spin = new_values['ospin'] | units.km / units.s
        star.epoch = new_values['epoch'] | units.Myr
        star.current_time = new_values['tphys']| units.Myr
        star.main_sequence_lifetime = new_values['tm'] | units.Myr
        star.new_time = new_values['tphysf']| units.Myr
        
    def get_time_step_for_star(self, star):
        
        current_values = {}
        current_values['kw'] = star.type.number
        current_values['mass'] = star.initial_mass.in_(units.MSun).number
        current_values['mt'] = star.mass.in_(units.MSun).number
        current_values['tm'] = star.main_sequence_lifetime.in_(units.Myr).number
        current_values['age'] = star.current_time.in_(units.Myr).number
        current_values['epoch'] = star.epoch.in_(units.Myr).number
        
        result = self.get_time_step(**current_values)
        
        return result | units.Myr
        
        
