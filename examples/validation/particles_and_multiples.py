from amuse.couple import bridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.newsmallN.interface import SmallN

from amuse.ic import plummer
from amuse.ic import gasplummer

from amuse.couple import multiples

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system

from optparse import OptionParser

import numpy
import time
try:
	import pylab
except ImportError:
	pylab = None

class MultiplesClusterCode(object):


    def __init__(self,
        nstars = 10, 
        endtime = 10,
        total_mass = 1000,
        rscale = 1.0,
        interaction_radius = -1.0,
        star_code = 'hermite',
        star_smoothing_fraction = 0.0,
        seed = -1,
        ntimesteps = 10,
        must_do_plot = True,
        **ignored_options
    ):
        if seed >= 0:
            numpy.random.seed(seed)
        
        if interaction_radius < 0.0:
            self.interaction_radius = 0.01 | nbody_system.length
        else:
            self.interaction_radius = interaction_radius | nbody_system.length
                
        self.must_do_plot = must_do_plot
        self.line = None
        self.line2 = None
        
        self.ntimesteps = ntimesteps
        self.nstars = nstars
        
        self.total_mass = total_mass | nbody_system.mass
        self.rscale = rscale |  nbody_system.length
        
        self.star_epsilon = star_smoothing_fraction * self.rscale
        
        self.star_mass = self.total_mass
        
        self.endtime = endtime | nbody_system.time
        self.delta_t = self.endtime / self.ntimesteps
        
        self.create_code(star_code)
        
        self.code = multiples.Multiples(self.star_code, self.new_smalln)
        
        time = 0
        sum_energy = self.code.kinetic_energy + self.code.potential_energy
        energy0 = sum_energy.value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())

        print "Time          :", time
        print "Energy        :", energy0
        print "Virial radius :", coreradius
        
        self.evolve_model()
        
        if must_do_plot:
            pylab.show()
            pylab.savefig(
                "multiples-{0}-{1}.png".format(
                    star_code,
                    nstars
                )
            )
        
        
        time = self.code.model_time.value_in(nbody_system.time)
        sum_energy = self.code.kinetic_energy + self.code.potential_energy - self.code.multiples_energy_correction
        energy = sum_energy.value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())
        
        print "Time          :", time
        print "Energy        :", energy
        print "Delta E       :", (energy-energy0)/energy0
        print "Virial radius :", coreradius
        
        self.stop()
        
        if must_do_plot:
            raw_input('Press enter...') 
        
   
    def update_plot(self, time, code):
        
        time = time.value_in(nbody_system.time), 
        sum_energy = code.kinetic_energy + code.potential_energy - self.code.multiples_energy_correction
        energy = sum_energy.value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())
       
        if self.line is None:
            pylab.ion()
            pylab.subplot(1,2,1)
            self.line = pylab.plot([time], [energy])[0]
            pylab.xlim(0, self.endtime.value_in(nbody_system.time))
            pylab.ylim(energy * 0.8, energy * 1.2)
            pylab.subplot(1,2,2)
            self.line2 = pylab.plot([time], [coreradius])[0]
            #self.line2 = pylab.plot([time], [kicke])[0]
            pylab.xlim(0, self.endtime.value_in(nbody_system.time))
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
        particles=plummer.new_plummer_model(self.nstars)
        particles.scale_to_standard()
        particles.radius= self.interaction_radius
        return particles
        
    def evolve_model(self):
        
        if self.must_do_plot:
            self.update_plot(time = 0 * self.delta_t, code = self.code)
            
        for time in self.delta_t * range(1,self.ntimesteps+1):
            self.code.evolve_model(time)
            print self.code.model_time
            if self.must_do_plot:
                self.update_plot(time = self.code.model_time, code = self.code)
            self.code.print_trees_summary()
            
        
    def create_code(self, star_code):
        self.star_code = getattr(self,'new_star_code_'+star_code)()
        
    
    def stop(self):
        self.star_code.stop()
    
    def new_smalln(self):
        result = SmallN()
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 50000
        return result
        
        
    def new_star_code_fi(self):
        result = Fi()
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = False
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
        
    def new_star_code_hermite(self):
        result = Hermite()
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_phigrape(self):
        result = PhiGRAPE(mode="gpu")
        result.parameters.initialize_gpu_once = 1
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_ph4(self):
        result = ph4(mode="gpu")
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_bhtree(self):
        result = BHTree()
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_octgrav(self):
        result = Octgrav()
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
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
        "--star-code", 
        default = "hermite",
        dest="star_code",
        help="the code modelling the particles ('hermite', 'bhtree', 'octgrav', 'phigrape')",
        type="string"
    )
    result.add_option(
        "-m", "--total-mass", 
        default = 1,
        dest="total_mass",
        help="the total mass in nbody units",
        type="float"
    )
    result.add_option(
        "-r", "--rscale", 
        default = 1.0,
        dest="rscale",
        help="length scale of the problem in nbody units (default 1) ",
        type="float"
    )
    result.add_option(
        "--star_smoothing_fraction", 
        default = 0.0,
        dest="star_smoothing_fraction",
        help="smoothing length of the stars as a fraction of the length scale (default 0)",
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
        "--interaction-radius", 
        default = 0.01,
        dest="interaction_radius",
        help="radius of all stars, collision detection will depend on this",
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
    return result
    
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    options = options.__dict__
    MultiplesClusterCode(**options)
    
