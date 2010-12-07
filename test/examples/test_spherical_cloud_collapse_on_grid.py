"""
In this script we simulate the evrard cload collapse in 3d.
"""

from amuse.support.core import late
from amuse.support.data.values import VectorQuantity
from amuse.support.data.core import Grid
from amuse.support import io
from amuse.support.io import text
from amuse.support.units.generic_unit_system import *
from amuse.support.units.generic_unit_converter import *
from amuse.support.units import constants
from amuse.support.units import units
from amuse.ext import cloud
from amuse.legacy.athena.interface import Athena, AthenaInterface
from amuse.legacy.capreole.interface import Capreole

try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False
    
    
from numpy import sqrt, arange, searchsorted, tanh, pi
from optparse import OptionParser


Grid.add_global_vector_attribute("position", ["x","y","z"])


class CalculateSolutionIn3D(object):
    number_of_workers = 1
    number_of_grid_points = 25
    gamma = 5.0/3.0
    name_of_the_code = "athena_selfgravity"
    convert_generic_units = ConvertBetweenGenericAndSiUnits(1.0 | units.kpc, 1.0e10 | units.MSun, constants.G)

    def __init__(self, **keyword_arguments):
        for x in keyword_arguments:
            print x, keyword_arguments[x]
            setattr(self, x, keyword_arguments[x])
            
        self.dimensions_of_mesh = (
            self.number_of_grid_points * 3, 
            self.number_of_grid_points, 
            self.number_of_grid_points
        )
        
        self.instance = None
        
    def new_instance_of_code(self):
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_athena_selfgravity_code(self):
        #unit_converter = self.convert_generic_units,
        result=Athena(mode = AthenaInterface.MODE_SELF_GRAVITY, redirection="none", number_of_workers=1)
        result.initialize_code()
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.3
        result.set_four_pi_G( 4 * pi )
        result.set_grav_mean_rho(1e-04)
        return result
        
    def new_instance_of_capreole_code(self):
        result=Capreole(number_of_workers=self.number_of_workers)
        result.initialize_code()
        return result
        
    def set_parameters(self, instance):
        
        instance.parameters.mesh_size = self.dimensions_of_mesh
        
        instance.parameters.length_x = 1 | length
        instance.parameters.length_y = 1 | length
        instance.parameters.length_z = 1 | length
        
        instance.x_boundary_conditions = ("periodic","periodic")
        instance.y_boundary_conditions = ("periodic","periodic")
        instance.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
    
    def new_grid(self):
        
        density = mass / length**3
        
        density = density
        momentum =  speed * density
        energy =  mass / (time**2 * length)
        
        grid = Grid(*self.dimensions_of_mesh)
        
        grid.rho =  1e-04 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        
    
        return grid
    
    def initialize_grid_with_sphere(self, grid):
        density = mass / length**3
        energy =  mass / (time**2 * length)
        
        center = [0.5, 0.5, 0.5] | length
        radii = (grid.position - center).lengths()
        rho0 = 0 | density
        rho1 = 1.0 | density
        radius = 0.1 | length
        total_mass = 1.0 | mass
        sig0 = 1.0
        rho_mean = 3.0 *  total_mass  / (4.0 * pi * (radius ** 3))
        rho_mean = 0.01 | density
        p0 = 0.1 | energy
        selection = radii < radius
        #grid.rho[selection] = (1.0 | density) * ((1.0 | length) / radii[selection]) 
        #grid.rho[selection] = grid.rho[selection] / grid.rho.sum().number
        
        grid.rho[selection] = rho_mean
        #print  "SUM:", grid.rho.sum()
        #print  "SUM:", grid.rho[selection]
        #grid.energy =  (self.gamma - 1.0) * grid.rho  * (0.05 | length ** 2 / time ** 2)
        #print  4 * pi / (self.gamma ** 2) * rho_mean ** 2 * radius ** 2  * (0.05 | ( length**3 / (mass * time**2)))
        grid.energy[selection] =  (4 * pi / (self.gamma ** 2)) * rho_mean ** 2 * radius ** 2  * (0.10 | ( length**3 / (mass * time**2)))
        #print grid.energy[selection]
    
    def setup_code(self):
        print "setup code"
        
        self.instance=self.new_instance_of_code()
        self.set_parameters(self.instance)
        
        
        self.grid = self.new_grid()
        
        self.from_model_to_code = self.grid.new_channel_to(self.instance.grid)
        self.from_code_to_model = self.instance.grid.new_channel_to(self.grid)
        
        self.from_code_to_model.copy_attributes(("x","y","z",))
        self.initialize_grid_with_sphere(self.grid)
        self.from_model_to_code.copy()
        
        self.instance.initialize_grid()
        
    def get_solution_at_time(self, time):
        if self.instance is None:
            self.setup_code()
        
        print "start evolve"
        self.instance.evolve(time)
        
        print "copying results"
        self.from_code_to_model.copy()

        return self.grid

def store_attributes(x, rho, rhovx, energy, filename):
    output = text.CsvFileText(filename = filename)
    output.quantities = (x, rho, rhovx, energy)
    output.attribute_names = ("x", "rho", "rhovx", "energy")
    output.store()


def store_attributes_of_line(grid, yindex = 0, zindex = 0, **options):
    store_attributes(
        grid.x[...,yindex,zindex],
        grid.rho[...,yindex,zindex],
        grid.rhovx[...,yindex,zindex],
        grid.energy[...,yindex,zindex],
        filename = "spherical_collapse_{name_of_the_code}_{number_of_grid_points}_{number_of_workers}_{time}.csv".format(**options)
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

def xtest_sperical_collapse():
    exact = CalculateExactSolutionIn1D()
    x, rho, p, u = exact.get_solution_at_time(0.12 | time)
    
    model = CalculateSolutionIn3D()
    model.name_of_the_code = "athena"
    model.dimensions_of_mesh = (500,1,1)
    
    grid = model.get_solution_at_time(0.12 | time)
    model_x = grid.x[...,0,0]
    density = grid.rho[...,0,0]
    
    index_in_model = searchsorted(model_x.value_in(length), 0.56)
    index_in_exact = searchsorted(x.value_in(length), 0.56)
    
    #store_attributes_of_line(grid, name_of_the_code = "athena-test", number_of_grid_points = 500, number_of_workers = 1)
    
    #assert abs((rho[index_in_exact] - density[index_in_model])/ density[index_in_model]) < 1.e-3 |units.none
    
    
def main(**options):
    #print "calculating shock using exact solution"
    #exact = CalculateExactSolutionIn1D()
    #x, rho, p, u = exact.get_solution_at_time(0.12 | time)
    
    center = options['number_of_grid_points'] / 2
    print center
    
    print "calculating shock using code"
    model = CalculateSolutionIn3D(**options)
    for t in range(50):
        grid = model.get_solution_at_time(t * 0.1 | time)
        print "saving data"
        store_attributes_of_line(grid, yindex = center, zindex = center, time = t, **options)
    


if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)
