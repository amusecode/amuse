"""
Simultaniously evolves gas and stars in a cluster.

The stars and gas particles evolve in their respective codes (Hermite and Gadget2)
The force of the gas particles on the star particles is calculated by a seperate code (BHTree)
The force of the star particles on the gas particles is calculated by the star code (Hermite)

The bridge works in S.I. units, but for this example all times and energies are
reported in nbody units.
"""

from amuse.couple import bridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.gadget2.interface import Gadget2

from amuse.ic import plummer
from amuse.ic import gasplummer

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system

import numpy
from matplotlib import pyplot

class BridgeStarAndGasPlummerCode(object):
    """
    Calculates a cluster with stars and gas, using a fixed set of codes
    (Hermite for stars, Gadget2 for gas and BHTree for gas to star interaction).
    """

    def __init__(self,
        nstars = 10, 
        ngas = None, 
        total_mass = 1000  | units.MSun,
        gas_fraction = 0.9,
        rscale = 1.0  | units.parsec,
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = None,
        interaction_timestep = 0.01 | nbody_system.time,
        diagnostic_timestep = 0.1 | nbody_system.time
    ):
        
        if not seed is None:
            numpy.random.seed(seed)
    
        if ngas is None:
            ngas = nstars * 10
        
        self.ngas = ngas
        self.nstars = nstars
        self.total_mass = total_mass
        self.gas_fraction = gas_fraction
        self.star_fraction = 1.0 - self.gas_fraction
        self.rscale = rscale
        
        self.star_epsilon = star_smoothing_fraction * self.rscale
        self.gas_epsilon = gas_smoothing_fraction * self.rscale
        
        self.star_mass = self.star_fraction * self.total_mass
        self.gas_mass = self.gas_fraction * self.total_mass
        
        self.converter = nbody_system.nbody_to_si(self.total_mass, self.rscale)
        self.interaction_timestep_in_si = self.converter.to_si(interaction_timestep)
        self.diagnostic_timestep_in_si = self.converter.to_si(diagnostic_timestep)
        
        self.time  = 0.0 * self.interaction_timestep_in_si
        
        self.times = [] | nbody_system.time
        self.energy_at_time = [] | nbody_system.energy
        self.coreradius_at_time = [] | self.rscale.unit    
        
    def setup(self):
        self.create_codes()
        
        self.initialize_data()
        
        self.create_bridge()
        
        self.code = self.bridge_system
        
        self.print_log() 
        
        self.store_values()
   
        
    def create_codes(self):
        self.star_code = self.new_star_code_hermite()
        self.gas_code = self.new_gas_code_gadget()
        
        def new_code_to_calculate_gravity_of_gas_particles():
            result = BHTree(self.converter)
            return result
            
        calculate_gravity_code=bridge.CalculateFieldForCodes(
            new_code_to_calculate_gravity_of_gas_particles,  # the code that calculates the acceleration field
            input_codes = [self.gas_code] # the codes to calculate the acceleration field of
        )
        self.gas_to_star_codes = [calculate_gravity_code]
        self.star_to_gas_codes = [self.star_code]
    
    def initialize_data(self):
        self.particles = self.new_particles_cluster()
        self.gas_particles = self.new_gas_cluster()
        
        self.star_code.particles.add_particles(self.particles)
        self.gas_code.gas_particles.add_particles(self.gas_particles)
        
        self.channel_from_code_to_memory_for_stars = self.star_code.particles.new_channel_to(self.particles)
        self.channel_from_code_to_memory_for_gas = self.gas_code.gas_particles.new_channel_to(self.gas_particles)
        
    def create_bridge(self):
        self.bridge_system = bridge.Bridge(
            timestep = self.interaction_timestep_in_si,
            use_threading = False
        )
        
        self.bridge_system.add_system(
            self.gas_code, 
            self.star_to_gas_codes
        )
        self.bridge_system.add_system(
            self.star_code, 
            self.gas_to_star_codes
        )
    
    def stop(self):
        self.star_code.stop()
        self.gas_code.stop()
        
    def new_gas_code_gadget(self):
        result = Gadget2(self.converter)
        return result
        
    def new_star_code_hermite(self):
        result = Hermite(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        return result
        
    def new_particles_cluster(self):
        particles=plummer.new_plummer_model(self.nstars,convert_nbody=self.converter)
        particles.radius = self.star_epsilon
        particles.mass = (1.0/self.nstars) * self.star_mass
        return particles

    def new_gas_cluster(self):
        particles=gasplummer.new_plummer_gas_model(self.ngas,convert_nbody=self.converter)
        particles.h_smooth = self.gas_epsilon
        particles.mass = (1.0/self.ngas) * self.gas_mass
        return particles
        
    def evolve_model(self, time_end):

        time_end_in_si = self.converter.to_si(time_end)
        while self.time < time_end_in_si:
            self.time += self.diagnostic_timestep_in_si
            self.code.evolve_model(self.time)
            print "evolved to time:", self.converter.to_nbody(self.code.time)
            self.store_values()
            
        self.channel_from_code_to_memory_for_stars.copy()
        self.channel_from_code_to_memory_for_gas.copy()
        
        self.print_log()
        
    def print_log(self):
        time = self.converter.to_nbody(self.time).value_in(nbody_system.time)
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())
        
        print "Time          :", time
        print "Energy        :", energy
        print "Virial radius :", coreradius
        
    def store_values(self):
        
        time = self.converter.to_nbody(self.time)
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy)
        coreradius = self.star_code.particles.virial_radius()
        
        self.times.append(time)
        self.energy_at_time.append(energy)
        self.coreradius_at_time.append(coreradius)
    
if __name__ in ("__main__", "__plot__"):
    code = BridgeStarAndGasPlummerCode(
        nstars = 100,
        diagnostic_timestep = 0.5 | nbody_system.time
    )
    code.setup()
    code.evolve_model(4.0 | nbody_system.time)
    
    figure = pyplot.figure(figsize= (5,10))
    
    subplot = figure.add_subplot(2, 1, 1)
    subplot.plot(
        code.times.value_in(nbody_system.time),
        code.coreradius_at_time.value_in(units.parsec)
    )
    subplot.set_title('code radius')
    subplot.set_xlabel('time (nbody)')
    subplot.set_ylabel('radius (parsec)')
    subplot = figure.add_subplot(2, 1, 2)
    subplot.plot(
        code.times.value_in(nbody_system.time),
        (code.energy_at_time - code.energy_at_time[0])/code.energy_at_time[0],
        'g'
    )
    subplot.set_title('energy error')
    subplot.set_xlabel('time (nbody)')
    subplot.set_ylabel('(E-E0)/E')
    pyplot.show()
