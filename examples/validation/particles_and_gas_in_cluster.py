from amuse.couple import bridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.phiGRAPE.interface import PhiGRAPE

from amuse.ic import plummer
from amuse.ic import gasplummer

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system

from optparse import OptionParser

import numpy
import time
import pylab





class GasPlummerModelExternalField(object):
    """
      skeleton grav code for use in bridge.

      must have get_gravity_at_point and get_potential_at_point
    """
    def __init__(self, position = [0,0,0] | units.parsec, radius=1000.| units.parsec, total_mass=1.6e10 | units.MSun):
        self.radius = radius
        self.total_mass = total_mass
        self.gravity_constant = constants.G
        self.position = position
        self.radius_squared = (self.radius**2)

    def get_gravity_at_point(self,eps,x,y,z):
        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radii_squared=dx**2 + dy**2 + dz**2
        #radii = radii_squared**0.5
        plummer_radii_squared = radii_squared + self.radius_squared
        plummer_radii_15 = plummer_radii_squared ** 1.5
        fr=-self.gravity_constant*self.total_mass/plummer_radii_15
        ax=fr*dx
        ay=fr*dy
        az=fr*dz
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radii_squared=dx**2 + dy**2 + dz**2
        #radii = radii_squared**0.5
        
        plummer_radii = (radii_squared + self.radius_squared)**0.5
        phi=self.gravity_constant*self.total_mass/plummer_radii
        return -phi * 2 
        
    def stop(self):
        pass

    @property
    def kinetic_energy(self):
        return quantities.zero

    @property
    def potential_energy(self):
        return quantities.zero
        
    @property
    def thermal_energy(self):
        return quantities.zero



class AbstractStarAndGasPlummerCode(object):


    def __init__(self,
        nstars = 10, 
        ngas = -1, 
        endtime = 10,
        total_mass = 1000,
        gas_fraction = 0.9,
        rscale = 1.0,
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = -1,
        ntimesteps = 10,
        must_do_plot = True
    ):
        
        if seed >= 0:
            numpy.random.seed(seed)
    
        if ngas < 0:
            ngas = nstars * 10
            
        self.must_do_plot = must_do_plot
        self.line = None
        self.line2 = None
        
        self.ntimesteps = ntimesteps
        self.ngas = ngas
        self.nstars = nstars
        
        self.total_mass = total_mass | units.MSun
        self.gas_fraction = gas_fraction
        self.star_fraction = 1.0 - self.gas_fraction
        self.rscale = rscale | units.parsec
        
        self.star_epsilon = star_smoothing_fraction * self.rscale
        self.gas_epsilon = gas_smoothing_fraction * self.rscale
        
        self.star_mass = self.star_fraction * self.total_mass
        self.gas_mass = self.gas_fraction * self.total_mass
        
        self.converter = nbody_system.nbody_to_si(self.total_mass, self.rscale)
        
        self.endtime = self.converter.to_si(endtime | nbody_system.time)
        self.delta_t = self.endtime / self.ntimesteps
        
    def update_plot(self, time, code):
        
        time = self.converter.to_nbody(time).value_in(nbody_system.time), 
        sum_energy = code.kinetic_energy + code.potential_energy + code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())
        kicke =  self.converter.to_nbody(code.kick_energy).value_in(nbody_system.energy)
       
        if self.line is None:
            pylab.ion()
            pylab.subplot(1,2,1)
            self.line = pylab.plot([time], [energy])[0]
            pylab.xlim(0, self.converter.to_nbody(self.endtime).value_in(nbody_system.time))
            pylab.ylim(energy * 0.8, energy * 1.2)
            pylab.subplot(1,2,2)
            self.line2 = pylab.plot([time], [coreradius])[0]
            #self.line2 = pylab.plot([time], [kicke])[0]
            pylab.xlim(0, self.converter.to_nbody(self.endtime).value_in(nbody_system.time))
            pylab.ylim(0,3)
            #pylab.ylim(-0.1, 0.1)
        else:
            xdata = self.line.get_xdata()
            ydata = self.line.get_ydata()
            xdata = numpy.concatenate( (xdata, time) )
            ydata = numpy.concatenate( (ydata, [energy]) )
            self.line.set_xdata(xdata)
            self.line.set_ydata(ydata)
            
            
            xdata = self.line2.get_xdata()
            ydata = self.line2.get_ydata()
            xdata = numpy.concatenate( (xdata, time) )
            #ydata = numpy.concatenate( (ydata, [kicke]) )
            ydata = numpy.concatenate( (ydata, [coreradius]) )
            self.line2.set_xdata(xdata)
            self.line2.set_ydata(ydata)
            
            pylab.draw()
    
    def new_particles_cluster(self):
        particles=plummer.new_plummer_sphere(self.nstars,convert_nbody=self.converter)
        particles.radius= self.star_epsilon
        particles.mass = (1.0/self.nstars) * self.star_mass
        return particles

    def new_gas_cluster(self):
        particles=gasplummer.new_plummer_gas_model(self.ngas,convert_nbody=self.converter)
        particles.h_smooth= self.gas_epsilon
        particles.mass = (1.0/self.ngas) * self.gas_mass
        return particles
        
    def new_particles_cluster_as_gas(self):
        particles=plummer.new_plummer_sphere(self.ngas,convert_nbody=self.converter)
        particles.radius= self.gas_epsilon
        particles.mass = (1.0/self.ngas) * self.gas_mass
        return particles
        
    def stop(self):
        pass
    
    def evolve_model(self):
        
        if self.must_do_plot:
            self.update_plot(time = 0 * self.delta_t, code = self.code)
            
        for time in self.delta_t * range(1,self.ntimesteps+1):
            self.code.evolve_model(time)
            print self.converter.to_nbody(self.code.time)
            if self.must_do_plot:
                self.update_plot(time = self.bridge_system.time, code = self.code)
        
        
        


class BridgeStarAndGasPlummerCode(AbstractStarAndGasPlummerCode):


    def __init__(self,
        nstars = 10, 
        ngas = -1, 
        endtime = 10,
        total_mass = 1000,
        gas_fraction = 0.9,
        rscale = 1.0,
        star_code = 'hermite',
        gas_code = 'field', 
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = -1,
        ntimesteps = 10,
        interaction_timestep = 0.01,
        must_do_plot = True,
        gas_to_star_interaction_code = 'none',
        star_to_gas_interaction_code = 'none',
        **ignored_options
    ):
        
        AbstractStarAndGasPlummerCode.__init__(
            self,
            nstars, 
            ngas, 
            endtime,
            total_mass,
            gas_fraction,
            rscale,
            star_smoothing_fraction,
            gas_smoothing_fraction,
            seed,
            ntimesteps,
            must_do_plot
        )

        
        self.interaction_timestep = self.converter.to_si(interaction_timestep| nbody_system.time)
        
        self.create_codes(
            gas_code,
            star_code,
            gas_to_star_interaction_code,
            star_to_gas_interaction_code,
        )
        
        self.create_bridge()
        
        self.code = self.bride_system
        
        self.evolve_model()
        
        if must_do_plot:
            pylab.show()
            pylab.savefig(
                "{0}-{1}-{2}-{3}.png".format(
                    star_code,
                    gas_code,
                    nstars,
                    ngas
                )
            )
        
        self.stop()
        
        if must_do_plot:
            raw_input('Press enter...') 
        
   
        
    def create_codes(self, gas_code, star_code, gas_to_star_interaction_code, star_to_gas_interaction_code):
        self.gas_code = getattr(self, 'new_gas_code_'+gas_code)()
        self.star_code = getattr(self,'new_star_code_'+star_code)()
        self.gas_to_star_codes = getattr(self, 'new_gas_to_star_interaction_codes_'+gas_to_star_interaction_code)(self.gas_code)
        self.star_to_gas_codes = getattr(self, 'new_star_to_gas_interaction_codes_'+star_to_gas_interaction_code)(self.star_code)
        
    def create_bridge(self):
        bridge_code1 = bridge.GravityCodeInField(
            self.gas_code, self.star_to_gas_codes
        )
        bridge_code2 = bridge.GravityCodeInField(
            self.star_code, self.gas_to_star_codes
        )
        
        self.bridge_system = bridge.Bridge(
            timestep = self.interaction_timestep,
            use_threading = False
        )
        self.bridge_system.add_code(bridge_code2)
        self.bridge_system.add_code(bridge_code1)
    
    def stop(self):
        self.star_code.stop()
        self.gas_code.stop()
    
    def new_gas_to_star_interaction_codes_self(self, gas_code):
        return [gas_code]
        
    def new_star_to_gas_interaction_codes_self(self, star_code):
        return [star_code]
        
    def new_gas_to_star_interaction_codes_none(self, gas_code):
        return []
        
    def new_star_to_gas_interaction_codes_none(self, gas_code):
        return []
        
    def new_gas_to_star_interaction_codes_octgrav(self, gas_code):
        def new_octgrav():
            result = Octgrav(self.converter)
            result.parameters.epsilon_squared = self.gas_epsilon ** 2
            return result
            
        return [bridge.CalculateFieldForCodes(new_octgrav, [gas_code])]
    
    
    def new_gas_to_star_interaction_codes_bhtree(self, gas_code):
        def new_bhtree():
            result = BHTree(self.converter)
            result.parameters.epsilon_squared = self.star_epsilon ** 2
            return result
            
        return [bridge.CalculateFieldForCodes(new_bhtree, [gas_code])]\
    
    def new_gas_code_fi(self):
        result = Fi(self.converter)
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = True
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        #result.parameters.adaptive_smoothing_flag = True
        #result.parameters.epsilon_squared = self.gas_epsilon ** 2
        #result.parameters.eps_is_h_flag = False
        result.parameters.integrate_entropy_flag = False
        
        #result.parameters.self_gravity_flag = False
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_gadget(self):
        result = Gadget2(self.converter)
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_field(self):
        result = GasPlummerModelExternalField(
            radius = self.rscale,
            total_mass = self.gas_mass
        )
        return result
        
    def new_gas_code_hermite(self):
        result = Hermite(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster_as_gas())
        result.commit_particles()
        return result
        
    def new_star_code_hermite(self):
        result = Hermite(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_phigrape(self):
        result = PhiGRAPE(self.converter, mode="gpu")
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_bhtree(self):
        result = BHTree(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_bhtree(self):
        result = BHTree(self.converter)
        result.parameters.epsilon_squared = self.gas_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster_as_gas())
        result.commit_particles()
        return result



class AllInOneStarAndGasPlummerCode(AbstractStarAndGasPlummerCode):


    def __init__(self,
        nstars = 10, 
        ngas = -1, 
        endtime = 10,
        total_mass = 1000,
        gas_fraction = 0.9,
        rscale = 1.0,
        sph_code = 'fi',
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = -1,
        ntimesteps = 10,
        must_do_plot = True,
        **ignored_options
    ):
        
        AbstractStarAndGasPlummerCode.__init__(
            self,
            nstars, 
            ngas, 
            endtime,
            total_mass,
            gas_fraction,
            rscale,
            star_smoothing_fraction,
            gas_smoothing_fraction,
            seed,
            ntimesteps,
            must_do_plot
        )

        
        self.interaction_timestep = self.converter.to_si(interaction_timestep| nbody_system.time)
        
        self.create_code(sph_code)
        
        self.evolve_model()
        
        if must_do_plot:
            pylab.show()
            pylab.savefig(
                "{0}-{1}-{2}.png".format(
                    sph_code,
                    nstars,
                    ngas
                )
            )
        
        self.stop()
        
        if must_do_plot:
            raw_input('Press enter...') 
        
    def create_code(self, name):
        self.code = getattr(self, 'new_sph_code_'+name)()
        
    def stop(self):
        self.code.stop()
    
    def new_sph_code_fi(self):
        result = Fi(self.converter)
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = True
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        
        #result.parameters.adaptive_smoothing_flag = True
        #result.parameters.epsilon_squared = self.gas_epsilon ** 2
        #result.parameters.eps_is_h_flag = False
        
        result.parameters.integrate_entropy_flag = False
        result.particles.add_particles(self.new_particles_cluster())
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_sph_code_gadget(self):
        result = Gadget2(self.converter)
        result.particles.add_particles(self.new_particles_cluster())
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n", "--nstar", 
        default = 10,
        dest="nstars",
        help="number of star particles",
        type="int"
    )
    result.add_option(
        "-g", "--ngas", 
        default = -1,
        dest="ngas",
        help="number of gas particles (if -1, 10 times the number of stars)",
        type="int"
    )
    result.add_option(
        "--gas-code", 
        default = "field",
        dest="gas_code",
        help="the code modelling the gas ('fi', 'gadget', 'field')",
        type="string"
    )
    result.add_option(
        "--star-code", 
        default = "hermite",
        dest="star_code",
        help="the code modelling the particles ('hermite', 'bhtree', 'octgrav', 'phigrape')",
        type="string"
    )
    result.add_option(
        "--sph-code", 
        default = "fi",
        dest="sph_code",
        help="the code modelling the particles and the gas simultaniously",
        type="string"
    )
    result.add_option(
        "--gas-star-code", 
        default = "self",
        dest="gas_to_star_interaction_code",
        help="the code calculating the gravity field of the gas code for the star code (default is self, gas code will calculate field for star code)",
        type="string"
    )
    result.add_option(
        "--star-gas-code", 
        default = "self",
        dest="star_to_gas_interaction_code",
        help="the code calculating the gravity field of the star code for the gas code (default is self, star code will calculate field for gas code)",
        type="string"
    )
    result.add_option(
        "-m", "--total-mass", 
        default = 1000.0,
        dest="total_mass",
        help="the total mass in solar masses",
        type="float"
    )
    result.add_option(
        "--gas-fraction", 
        default = 0.9,
        dest="gas_fraction",
        help="the gas fraction between 0.0 and 1.0 (default 0.9)",
        type="float"
    )
    
    result.add_option(
        "-r", "--rscale", 
        default = 1.0,
        dest="rscale",
        help="length scale of the problem in parsec (default 1) ",
        type="float"
    )
    
    result.add_option(
        "--star_smoothing_fraction", 
        default = 0.001,
        dest="star_smoothing_fraction",
        help="smoothing length of the stars as a fraction of the length scale",
        type="float"
    )
    
    result.add_option(
        "--gas_smoothing_fraction", 
        default = 0.05,
        dest="gas_smoothing_fraction",
        help="smoothing length of the gas particles as a fraction of the length scale",
        type="float"
    )
    
    result.add_option(
        "-s", "--seed", 
        default = 0,
        dest="seed",
        help="random number seed (-1, no seed)",
        type="int"
    )
    result.add_option(
        "--interaction-timestep", 
        default = 0.01,
        dest="interaction_timestep",
        help="time between bridge interactions (0.01 nbody time)",
        type="float"
    )
    result.add_option(
        "-t", "--end-time", 
        default = 1,
        dest="endtime",
        help="end time of the simulation (in nbody time, default 1)",
        type="float"
    )
    result.add_option(
        "--ntimesteps", 
        default = 10,
        dest="ntimesteps",
        help="number of times to do reporting",
        type="int"
    )
    result.add_option(
        "--noplot", 
        dest="must_do_plot",
        default = True,
        help="do not show a plot and end as soon as possible",
        action="store_false"
    )
    result.add_option(
        "--allinoone", 
        dest="must_do_bridge",
        default = True,
        help="simulate the stars and gas with one sph code",
        action="store_false"
    )
    return result
    
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    if options.allinone:
        options = options.__dict__
        AllInOneStarAndGasPlummerCode(**options)
    else:
        options = options.__dict__
        BridgeStarAndGasPlummerCode(**options)
    
