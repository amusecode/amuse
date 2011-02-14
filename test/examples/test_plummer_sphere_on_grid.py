"""
In this script we simulate a plummer sphere on a grid
"""

from amuse.support.core import late
from amuse.support.data.values import VectorQuantity
from amuse.support.data.core import Grid, Particles
from amuse.support import io
from amuse.support.io import text
from amuse.support.units.generic_unit_system import *
from amuse.support.units.generic_unit_converter import *
from amuse.support.units import constants
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.ext import cloud
from amuse.community.athena.interface import Athena, AthenaInterface
from amuse.community.capreole.interface import Capreole
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.octgrav.interface import Octgrav

import sys

try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False
    
    
from numpy import sqrt, arange, searchsorted, tanh, pi
from optparse import OptionParser


#Grid.add_global_vector_attribute("position", ["x","y","z"])


class HydroGridAndNbody(object):
    
    def __init__(self, gridcode, nbodycode):
        self.gridcode = gridcode
        self.nbodycode = nbodycode
    
    def setup_positions_for_potential_grid(self):
        self.x = self.gridcode.potential_grid.x.flatten()
        self.y = self.gridcode.potential_grid.y.flatten()
        self.z = self.gridcode.potential_grid.z.flatten()
        self.eps =  self.x.aszeros()
        
    def setup_particles_in_nbodycode(self):
        staggered_grid_shape =  numpy.asarray(self.gridcode.grid.shape) + 2
        corner0 = self.gridcode.grid[0][0][0].position
        corner1 = self.gridcode.grid[-1][-1][-1].position
        
        delta = self.gridcode.grid[1][1][1].position - corner0
        print delta.prod()
        self.volume = delta.prod()
        staggered_corner0 = corner0 - delta
        staggered_corner1 = corner1
        self.staggered_grid = Grid.create(staggered_grid_shape, staggered_corner1 - staggered_corner0 + delta)
        #self.staggered_grid.x += staggered_corner0.x
        #self.staggered_grid.y += staggered_corner0.x
        #self.staggered_grid.z += staggered_corner0.x
        
        #self.staggered_grid.p000 = 0.0 | mass
        #self.staggered_grid.p100 = 0.0 | mass
        #self.staggered_grid.p010 = 0.0 | mass
        #self.staggered_grid.p110 = 0.0 | mass
        #self.staggered_grid.p000 = 0.0 | mass
        #self.staggered_grid.p100 = 0.0 | mass
        #self.staggered_grid.p011 = 0.0 | mass
        #self.staggered_grid.p111 = 0.0 | mass
        #self.staggered_grid.p001 = 0.0 | mass
        #self.staggered_grid.p101 = 0.0 | mass
        
        self.normal_grid = Grid.create(self.gridcode.grid.shape, corner1 - corner0 + delta)
        particles = Particles(self.normal_grid.size)
        particles.mass = 0.0 | mass
        particles.x = self.grid.x.flatten()
        particles.y = self.grid.y.flatten()
        particles.z = self.grid.z.flatten()
        particles.vx = 0.0 | speed
        particles.vy = 0.0 | speed
        particles.vz = 0.0 | speed
        particles.radius = 0.0 | length
        
        self.nbodycode.particles.add_particles(particles)
        self.from_model_to_nbody = particles.new_channel_to(self.nbodycode.particles)
        self.nbodycode.commit_particles()
        self.particles = particles
        
    @property
    def grid(self):
        return self.gridcode.grid
        
    @property
    def parameters(self):
        return self.gridcode.parameters
    
    def commit_parameters(self):
        self.gridcode.commit_parameters()
        self.setup_positions_for_potential_grid()
        self.setup_particles_in_nbodycode()
        
    def initialize_grid(self):
        self.gridcode.initialize_grid()
        self.position = self.gridcode.grid.position
        
        
        
    def evolve(self, time):
        #masses = self.gridcode.grid.rho * self.volume * 1.0 / 8.0
        #self.staggered_grid.p000[1:,1:,1:] = masses
        #self.staggered_grid.p100[:-1,1:,1:] = masses
        # self.staggered_grid.p010[1:,:-1,1:] = masses
        #self.staggered_grid.p110[:-1,:-1,1:] = masses
        #self.staggered_grid.p001[1:,1:,:-1] = masses
        #self.staggered_grid.p101[:-1,1:,:-1] = masses
        #self.staggered_grid.p011[1:,:-1,:-1] = masses
        #self.staggered_grid.p111[:-1,:-1,:-1] = masses
        
        #self.staggered_grid.masses = (
        #    self.staggered_grid.p000 + 
        #    self.staggered_grid.p100 + 
        #    self.staggered_grid.p010 +
        #    self.staggered_grid.p110 +
        #    self.staggered_grid.p001 +
        #    self.staggered_grid.p101 +
        #    self.staggered_grid.p011 +
        #    self.staggered_grid.p111
        #)
        
        #self.particles.mass = self.staggered_grid.masses.flatten()
        self.particles.mass = (self.gridcode.grid.rho * self.volume).flatten()
        self.from_model_to_nbody.copy_attribute('mass')
        
        correction = (length ** 3) / (mass * (time ** 2))
        print correction, nbody_system.G
        
        print "getting potential energy"
        potential = self.nbodycode.get_potential_at_point(
            self.eps,
            self.x,
            self.y,
            self.z
            
        )
        print potential.shape, self.gridcode.potential_grid.shape, self.grid.shape
        self.staggered_grid.rho  = 0.0 | density
        self.staggered_grid[1:-1,1:-1,1:-1].rho = self.gridcode.grid.rho
        correction = ((self.staggered_grid.rho * self.volume) / numpy.sqrt(self.nbodycode.parameters.epsilon_squared)) * nbody_system.G
        print correction.flatten().shape, potential.shape
        potential = potential + correction.flatten()
        print "got potential enery"
        potential = potential.reshape(self.gridcode.potential_grid.shape)
        self.gridcode.potential_grid.potential = potential 
        
        print "evolve hydro", time
        self.gridcode.evolve(time)
        print "end evolve hydro"


class HydroGridAndNbodyWithAccelerationTransfer(object):
    
    def __init__(self, gridcode, nbodycode):
        self.gridcode = gridcode
        self.nbodycode = nbodycode
        
        self.gridcode.grid.add_vector_attribute("acceleration_field", ["fx","fy","fz"])
    
    def setup_positions_for_potential_grid(self):
        self.x = self.gridcode.acceleration_grid.x.flatten()
        self.y = self.gridcode.acceleration_grid.y.flatten()
        self.z = self.gridcode.acceleration_grid.z.flatten()
        self.eps =  self.x.aszeros()
        
    def setup_particles_in_nbodycode(self):
        staggered_grid_shape =  numpy.asarray(self.gridcode.grid.shape) + 1
        corner0 = self.gridcode.grid[0][0][0].position
        corner1 = self.gridcode.grid[-1][-1][-1].position
        
        delta = self.gridcode.grid[1][1][1].position - corner0
        print delta.prod()
        self.volume = delta.prod()
        staggered_corner0 = corner0 - delta
        staggered_corner1 = corner1
        self.staggered_grid = Grid.create(staggered_grid_shape, staggered_corner1 - staggered_corner0 + delta)
        self.staggered_grid.x += staggered_corner0.x
        self.staggered_grid.y += staggered_corner0.x
        self.staggered_grid.z += staggered_corner0.x
        
        self.staggered_grid.p000 = 0.0 | mass
        self.staggered_grid.p100 = 0.0 | mass
        self.staggered_grid.p010 = 0.0 | mass
        self.staggered_grid.p110 = 0.0 | mass
        self.staggered_grid.p000 = 0.0 | mass
        self.staggered_grid.p100 = 0.0 | mass
        self.staggered_grid.p011 = 0.0 | mass
        self.staggered_grid.p111 = 0.0 | mass
        self.staggered_grid.p001 = 0.0 | mass
        self.staggered_grid.p101 = 0.0 | mass
        
        self.normal_grid = Grid.create(self.gridcode.grid.shape, corner1 - corner0 + delta)
        particles = Particles(self.staggered_grid.size)
        particles.mass = 0.0 | mass
        particles.x = self.normal_grid.x.flatten()
        particles.y = self.normal_grid.y.flatten()
        particles.z = self.normal_grid.z.flatten()
        particles.vx = 0.0 | speed
        particles.vy = 0.0 | speed
        particles.vz = 0.0 | speed
        particles.radius = 0.0 | length
        
        self.nbodycode.particles.add_particles(particles)
        self.from_model_to_nbody = particles.new_channel_to(self.nbodycode.particles)
        self.nbodycode.commit_particles()
        self.particles = particles
        
    @property
    def grid(self):
        return self.gridcode.grid
        
    @property
    def parameters(self):
        return self.gridcode.parameters
    
    def commit_parameters(self):
        self.gridcode.commit_parameters()
        self.setup_positions_for_potential_grid()
        self.setup_particles_in_nbodycode()
        
    def initialize_grid(self):
        self.gridcode.initialize_grid()
        self.position = self.gridcode.grid.position
        
        
        
    def evolve(self, time):
        masses = self.gridcode.grid.rho * self.volume * 1.0 / 8.0
        self.staggered_grid.p000[1:,1:,1:] = masses
        self.staggered_grid.p100[:-1,1:,1:] = masses
        self.staggered_grid.p010[1:,:-1,1:] = masses
        self.staggered_grid.p110[:-1,:-1,1:] = masses
        self.staggered_grid.p001[1:,1:,:-1] = masses
        self.staggered_grid.p101[:-1,1:,:-1] = masses
        self.staggered_grid.p011[1:,:-1,:-1] = masses
        self.staggered_grid.p111[:-1,:-1,:-1] = masses
        
        self.staggered_grid.masses = (
            self.staggered_grid.p000 + 
            self.staggered_grid.p100 + 
            self.staggered_grid.p010 +
            self.staggered_grid.p110 +
            self.staggered_grid.p001 +
            self.staggered_grid.p101 +
            self.staggered_grid.p011 +
            self.staggered_grid.p111
        )
        
        #self.particles.mass = self.staggered_grid.masses.flatten()
        self.particles.mass = self.gridcode.grid.rho * self.volume
        self.from_model_to_nbody.copy_attribute('mass')
        print "getting acceleration field"
        acc_x, acc_y, acc_z = self.nbodycode.get_gravity_at_point(
            self.eps,
            self.x,
            self.y,
            self.z
        )
        print "got acceleration field"
        acc_x  =  acc_x.reshape(self.grid.shape)
        acc_y  =  acc_y.reshape(self.grid.shape)
        acc_z  =  acc_z.reshape(self.grid.shape)
        self.gridcode.acceleration_grid._set_values(None, ["fx","fy","fz"], [acc_x, acc_y, acc_z]) 
        
        print "evolve hydro", time
        self.gridcode.evolve(time)
        print "end evolve hydro"
        
class CalculateSolutionIn3D(object):
    size = 10.0 | length
    
    number_of_workers = 1
    number_of_grid_points = 25
    gamma = 5.0/3.0
    rho_medium = 0 | density
    rho_sphere = 0.5 | density
    center = [1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0] * size
    radius = 1.0 | length
    total_mass = 1.0 | mass
    name_of_the_code = "athena_selfgravity"
    convert_generic_units = ConvertBetweenGenericAndSiUnits(1.0 | units.kpc, 1.0e10 | units.MSun, constants.G)

    def __init__(self, **keyword_arguments):
        for x in keyword_arguments:
            print x, keyword_arguments[x]
            setattr(self, x, keyword_arguments[x])
            
        self.dimensions_of_mesh = (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            self.number_of_grid_points
        )
        
        self.instance = None
        
    def new_instance_of_code(self):
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_athena_selfgravity_code(self):
        result=Athena(mode = AthenaInterface.MODE_SELF_GRAVITY, redirection="none", number_of_workers=1)
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.3
        result.set_four_pi_G( 4.0 * pi ) # G == 1
        result.set_grav_mean_rho(self.rho_mean.value_in(density))
        return result
        
    def new_instance_of_athena_code(self):
        result=Athena(redirection="none", number_of_workers=1)
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.3
        return result
        
    def new_instance_of_athena_and_nbody_code(self):
        gridcode=Athena(redirection="none", number_of_workers=1)
        gridcode.initialize_code()
        gridcode.parameters.gamma = self.gamma
        gridcode.parameters.courant_number=0.3
        
        gridcode.set_has_external_gravitational_potential(1)
        
        nbodycode = PhiGRAPE(mode="gpu")
        #nbodycode = Octgrav(mode="gpu")
        nbodycode.initialize_code()
        nbodycode.parameters.epsilon_squared = (self.size / (10000.0 * self.number_of_grid_points)) ** 2
        nbodycode.commit_parameters()
        
        result = HydroGridAndNbody(gridcode, nbodycode)
        return result
        
    def new_instance_of_capreole_and_nbody_code(self):
        gridcode=Capreole(number_of_workers=self.number_of_workers)
        gridcode.initialize_code()
        nbodycode = PhiGRAPE(mode="gpu")
        #nbodycode = Octgrav(mode="gpu")
        nbodycode.initialize_code()
        nbodycode.parameters.epsilon_squared = (self.size / (10000.0 * self.number_of_grid_points)) ** 2
        nbodycode.commit_parameters()
        
        result = HydroGridAndNbodyWithAccelerationTransfer(gridcode, nbodycode)
        return result
        
    def set_parameters(self, instance):
        
        instance.parameters.mesh_size = self.dimensions_of_mesh
        
        instance.parameters.length_x = self.size
        instance.parameters.length_y = self.size
        instance.parameters.length_z = self.size
        
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
    
    def new_grid(self):
        
        density = mass / length**3
        
        density = density
        momentum =  speed * density
        energy =  mass / (time**2 * length)
        
        grid = Grid.create(self.dimensions_of_mesh, (10.0, 10.0, 10.0) | length)
        
        grid.rho =  self.rho_medium
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy

        return grid
    
    def initialize_grid_with_plummer_sphere(self, grid):
        scaled_radius = self.radius #/ 1.695
        
        
        radii = (grid.position - self.center).lengths()
        selected_radii  = radii[radii < self.radius ]
        print "number of cells in cloud (number of cells in grid)", len(selected_radii), grid.shape, grid.size
        
        self.rho_sphere = ((0.75 * self.total_mass /  (pi * (scaled_radius ** 3))))
        grid.rho = self.rho_medium + (self.rho_sphere * ((1 + (radii ** 2) / (scaled_radius ** 2))**(-5.0/2.0)))
        
        internal_energy = (0.25 | potential) * ( 1.0 | length / mass) * self.total_mass / scaled_radius
        # 1.0 * grid.rho * 
        grid.energy =  grid.rho * (internal_energy/((1.0+(radii/scaled_radius)**2)**(1.0/2.0)))
        
    def setup_code(self):
        print "setup code"
        
        self.grid = self.new_grid()
        self.initialize_grid_with_plummer_sphere(self.grid)
        print "Mean density", self.grid.rho.flatten().mean()
        self.rho_mean = self.grid.rho.flatten().mean()
        
        self.instance=self.new_instance_of_code()
        self.set_parameters(self.instance)
        
        
        
        self.from_model_to_code = self.grid.new_channel_to(self.instance.grid)
        self.from_code_to_model = self.instance.grid.new_channel_to(self.grid)
        
        self.from_code_to_model.copy_attributes(("x","y","z",))
        self.from_model_to_code.copy()
        self.instance.initialize_grid()
        
    def get_solution_at_time(self, time):
        if self.instance is None:
            self.setup_code()
        
        print "start evolve"
        self.instance.evolve(time)
        
        print "copying results"
        self.from_code_to_model.copy()
        
        #print "max,min",  max(self.grid.rhovx.flatten()),  min(self.grid.rhovx.flatten())
        #print "max,min",  max(self.grid.rhovy.flatten()),  min(self.grid.rhovy.flatten())
        #print "max,min",  max(self.grid.rhovz.flatten()),  min(self.grid.rhovz.flatten())
        #print "max,min",  max(self.grid.energy.flatten()),  min(self.grid.energy.flatten())
        #print "max,min",  max(self.grid.rho.flatten()),  min(self.grid.rho.flatten())

        return self.grid

def store_attributes(x, rho, rhovx, energy, filename):
    output = text.TableFormattedText(filename = filename)
    output.quantities = (x, rho, rhovx, energy)
    output.attribute_names = ("x", "rho", "rhovx", "energy")
    output.store()


def store_attributes_of_line(grid, yindex = 0, zindex = 0, **options):
    store_attributes(
        grid.x[...,yindex,zindex],
        grid.rho[...,yindex,zindex],
        grid.rhovx[...,yindex,zindex],
        grid.energy[...,yindex,zindex],
        filename = "plummer_sphere_{name_of_the_code}_{number_of_grid_points}_{number_of_workers}_{time}.csv".format(**options)
    )
    

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n",
        "--mesh-size", 
        dest="number_of_grid_points",
        type="int",
        default = 25,
        help="number of grid cells in the x, y and z direction"
    )
    
    result.add_option(
        "-w",
        "--workers", 
        dest="number_of_workers",
        type="int",
        default = 1,
        help="number of parallel workers to start"
    )
    result.add_option(
        "-c",
        "--code",
        dest="name_of_the_code",
        default="athena_selfgravity",
        help="name of the code to use"
    )
    return result
    
def main(**options):
    center = options['number_of_grid_points'] / 2
    
    print "calculating shock using code"
    model = CalculateSolutionIn3D(**options)
    for t in range(50):
        grid = model.get_solution_at_time( t * 0.1 | time)
        print "saving data"
        store_attributes_of_line(grid, yindex = center, zindex = center, time = t, **options)
    


if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)
