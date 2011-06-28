"""
In this script we simulate Kelvin-Helmholtz Instability in two dimentsions.
"""
import numpy
from matplotlib import pyplot
    

from amuse.support.data.core import Grid
from amuse.support.units.generic_unit_system import *
from amuse.community.athena.interface import Athena

density = mass / length**3
momentum =  speed * density
energy =  mass / (time**2 * length)
magnetic_field = mass / current / time**2
vector_potential = magnetic_field * length

PI = numpy.pi

GAMMA = 5.0/3.0
DIMENSIONS_OF_MESH = (64,64,1)
DIMENSIONS_OF_A_MESH = (65,65,1)
B0 = 1.0/numpy.sqrt(4.0*PI) | magnetic_field
D0 = 25.0/(36.0*PI)   | density
V0 = 1.0              | speed
P0 = 5.0/(12.0*PI)    | energy

def new_instance_of_hydro_code(number_of_workers=1):
    result=Athena(number_of_workers = number_of_workers)
    result.initialize_code()
    result.parameters.gamma = GAMMA
    result.parameters.courant_number=0.8
    return result

def set_parameters(instance):
    instance.parameters.mesh_size = DIMENSIONS_OF_MESH

    instance.parameters.length_x = 1 | length
    instance.parameters.length_y = 1 | length
    instance.parameters.length_z = 1 | length

    instance.parameters.x_boundary_conditions = ("periodic","periodic")
    instance.parameters.y_boundary_conditions = ("periodic","periodic")
    instance.parameters.z_boundary_conditions = ("periodic","periodic")

    result = instance.commit_parameters()

def new_grid():
    grid = Grid.create(DIMENSIONS_OF_MESH, [1,1,1] | length)
    self.clear_grid(grid)
    return grid

def clear_grid(grid):
    grid.rho =  0.0 | density
    grid.rhovx = 0.0 | momentum
    grid.rhovy = 0.0 | momentum
    grid.rhovz = 0.0 | momentum
    grid.energy = 0.0 | energy
    grid.B1i    = 0.0 | magnetic_field
    grid.B2i    = 0.0 | magnetic_field
    grid.B3i    = 0.0 | magnetic_field

    return grid

def initialize_grid(grid):        

    temp = [1,1,1] | length
    l    = 1.0 | length;
    A_grid = Grid.create(DIMENSIONS_OF_A_MESH, temp + grid.cellsize()*[1,1,0])
    A_grid.position -= A_grid.cellsize()*[0.5,0.5,0.0]
    A_grid.Az = (
        B0*l/(4.0*PI)*(4.0*PI/l*A_grid.x).cos() + 
        B0*l/(2.0*PI)*(2.0*PI/l*A_grid.y).cos() ) 

    assert grid.cellsize().x == A_grid.cellsize().x
    assert grid.cellsize().y == A_grid.cellsize().y
#    assert grid.dim[0]      == A_grid.dim[0]

    grid.B1i =  (A_grid.Az[:-1,1:,...] - A_grid.Az[:-1,:-1,...])/A_grid.cellsize().y;
    grid.B2i = -(A_grid.Az[1:,:-1,...] - A_grid.Az[:-1,:-1,...])/A_grid.cellsize().x;
    grid.B3i =  0.0 | magnetic_field;

    grid.B1c = 0.0 | magnetic_field;
    grid.B2c = 0.0 | magnetic_field;
    grid.B3c = 0.0 | magnetic_field;
    grid.divB = 0.0 | (magnetic_field / length)

    grid.B1c[0:-1,...,...] = (grid.B1i[:-1,...,...] + grid.B1i[1:,...,...])*0.5
    grid.B2c[...,0:-1,...] = (grid.B2i[...,:-1,...] + grid.B2i[...,1:,...])*0.5
    grid.B3c = grid.B3i

    grid.B1c[-1,...,...] = grid.B1c[0,...,...]
    grid.B2c[...,-1,...] = grid.B2c[...,0,...]


    mu0 = 1.0 | (mass*length/time**2/current**2)

    grid.divB[0:-1,0:-1,...] = (
        (grid.B1i[:-1,0:-1,...] - grid.B1i[1:,0:-1,...])/grid.cellsize().x +
        (grid.B2i[0:-1,:-1,...] - grid.B2i[0:-1,1:,...])/grid.cellsize().y
        )

    assert (abs(grid.divB) < 1.0e-10 | grid.divB.unit).all()

    assert grid.B1i.unit == magnetic_field

    grid.rho   =  D0 
    grid.rhovx = -D0*V0*(2.0*PI/l*grid.y).sin()
    grid.rhovy =  D0*V0*(2.0*PI/l*grid.x).sin()
    grid.rhovz =  0.0 | momentum



    print grid.rhovx.sum()
    print grid.rhovy.sum()
    print grid.rhovz.sum()

    grid.energy = (
      P0 / (GAMMA - 1) + 
      0.5 * (grid.rhovx**2 + grid.rhovy**2 + grid.rhovz**2)/grid.rho +
      0.5 * (grid.B1c**2   + grid.B2c**2   + grid.B3c**2)/mu0
      )



def simulate_orszag_tang_problem(end_time):
    instance=new_instance_of_hydro_code()
    set_parameters(instance)

    print "setup grid"
    for x in instance.itergrids():
        inmem = x.copy_to_memory()

        clear_grid(inmem)
        initialize_grid(inmem)

        from_model_to_code = inmem.new_channel_to(x)
        from_model_to_code.copy()
    instance.initialize_grid()

    print "start evolve"
    dt = end_time / 10
    t = dt
    while t < end_time:
        instance.evolve_model(t)

        print "time : ", t
        t += dt

    print "copying results"
    result = []
    for x in instance.itergrids():
      result.append(x.copy_to_memory())

    print "terminating code"
    instance.stop()

    return result

def plot_grid(grid):
    rho = grid.rho[...,...,0].value_in(density)
    print rho
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(rho, origin = 'lower')
    figure.savefig('orszag_tang.png')
    pyplot.show()

if __name__ in ("__main__", "__plot__"):
    grids = simulate_orszag_tang_problem(0.2 | time)
    plot_grid(grids[0])
