"""
Runs the Kelvin-Helmholtz Instability problem in two dimensions with Athena.
"""
import numpy
from matplotlib import pyplot
#from amuse.lab import *
from amuse.community.athena.interface import Athena
from amuse.community.capreole.interface import Capreole
from amuse.units.generic_unit_system import *
from amuse.datamodel import Grid

GAMMA = 1.4
DIMENSIONS_OF_MESH = (200,200,1)
PERTUBATION_AMPLITUDE = 0.1 | speed

def new_instance_of_athena_code(number_of_workers=1):
    result=Athena(number_of_workers = number_of_workers)
    result.parameters.gamma = GAMMA
    result.parameters.courant_number=0.8
    set_parameters(result)
    return result

def new_instance_of_capreole_code(number_of_workers=1):
    result=Capreole(number_of_workers=number_of_workers)
    result.initialize_code()
        
    result.parameters.x_boundary_conditions = ("reflective","reflective")
    result.parameters.y_boundary_conditions = ("outflow","outflow")
    result.parameters.z_boundary_conditions = ("reflective","reflective")
    set_parameters(result)
    return result

def new_instance_of_mpiamrvac_code(number_of_workers=1):
    from amuse.community.mpiamrvac.interface import MpiAmrVac
    #result=MpiAmrVac(number_of_workers=number_of_workers) #, redirection="none")
    result=MpiAmrVac(mode="2d", number_of_workers=number_of_workers)
    result.set_parameters_filename(result.default_parameters_filename)
    result.initialize_code()
    result.parameters.maximum_number_of_grid_levels = 3
#    result.parameters.spatial_discretization_method = 'tvdmu'
#    result.parameters.predictor_step_discretization_method = 'tvdmu'
#    result.parameters.entropy_type = 'powell'
#    result.parameters.courant_number = 0.8
        
#    result.parameters.x_boundary_conditions = ("cont","cont")
#    result.parameters.y_boundary_conditions = ("cont","cont")
#    result.parameters.z_boundary_conditions = ("cont","cont")

    result.parameters.length_x = 1 | length
    result.parameters.length_y = 1 | length
    result.parameters.length_z = 1 | length
    
    result.parameters.x_boundary_conditions = ("periodic","periodic")
    result.parameters.y_boundary_conditions = ("periodic","periodic")
    result.parameters.z_boundary_conditions = ("periodic","periodic")
    
    return result

def set_parameters(instance):
    instance.parameters.mesh_size = DIMENSIONS_OF_MESH
    
    instance.parameters.length_x = 1 | length
    instance.parameters.length_y = 1 | length
    instance.parameters.length_z = 1 | length
    
    instance.parameters.x_boundary_conditions = ("periodic","periodic")
    instance.parameters.y_boundary_conditions = ("periodic","periodic")
    instance.parameters.z_boundary_conditions = ("periodic","periodic")
    
    
def new_grid():
    grid = Grid.create(DIMENSIONS_OF_MESH, [1,1,1] | length)
    clear_grid(grid)
    return grid
    
def clear_grid(grid):
    density = mass / length**3
    momentum =  speed * density
    energy =  mass / (time**2 * length)

    grid.rho =  0.0 | density
    grid.rhovx = 0.0 | momentum
    grid.rhovy = 0.0 | momentum
    grid.rhovz = 0.0 | momentum
    grid.energy = 0.0 | energy

    return grid
    
def initialize_grid(grid):        
    vx = 0.5 | speed
    p = 2.5 | (mass / (length * time**2))
    
    halfway = DIMENSIONS_OF_MESH[0]/2 - 1
    
    outerregion = numpy.logical_or(grid.y <= 0.25 | length, grid.y >= 0.75 | length)
    innerregion = numpy.logical_and(grid.y > 0.25 | length, grid.y < 0.75 | length)
    
    grid[outerregion].rho = 1  | density
    grid[outerregion].rhovx =  vx * grid[outerregion].rho
    
    grid[innerregion].rho = 2.0  | density
    grid[innerregion].rhovx = -vx * grid[innerregion].rho
    
    grid.energy = p / (GAMMA - 1)
        
def pertubate_grid(grid):
    grid.rhovx += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    grid.rhovy += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    
    grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
def simulate_kelvin_helmholtz_instability(end_time):
#    instance=new_instance_of_athena_code()
#    instance=new_instance_of_capreole_code()    
    instance=new_instance_of_mpiamrvac_code()    
#    set_parameters(instance)
    
    print("setup grid")
    for x in instance.itergrids():
        inmem = x.copy()
        
        clear_grid(inmem)
        initialize_grid(inmem)
        pertubate_grid(inmem)
        
        from_model_to_code = inmem.new_channel_to(x)
        from_model_to_code.copy()
    
    print("start evolve")
    dt = end_time / 10.0
    t = dt
    while t < end_time:
        instance.evolve_model(t)
        
        print(("time : ", t))
        t += dt
    
    print("copying results")
    result = []
    for x in instance.itergrids():
        result.append(x.copy())

    print("terminating code")
    instance.stop()

    return result
    
def plot_grid(grid):
    rho = grid.rho[...,0].value_in(density)
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(rho, origin = 'lower', cmap='jet')
    figure.savefig('kelvin_helmholtz.png')
    pyplot.show()
    
if __name__ == "__main__":
    grids = simulate_kelvin_helmholtz_instability(1.0 | time)
    plot_grid(grids[0])
