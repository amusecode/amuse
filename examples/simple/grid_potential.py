from amuse.lab import *

from matplotlib import pyplot

def hydro_grid_in_potential_well(mass = 1 | units.MSun, length = 1 | units.AU):
    converter = nbody_system.nbody_to_si(mass, length)
    
    instance=Capreole(converter)
    instance.initialize_code()
    instance.parameters.nx = 20
    instance.parameters.ny = 20
    instance.parameters.nz = 20
    instance.parameters.length_x = length
    instance.parameters.length_y = length
    instance.parameters.length_z =  length
    instance.parameters.x_boundary_conditions = ("outflow","outflow")
    instance.parameters.y_boundary_conditions = ("outflow","outflow")
    instance.parameters.z_boundary_conditions = ("outflow","outflow")
    
    instance.stopping_conditions.number_of_steps_detection.enable()
    #instance.set_has_external_gravitational_potential(1)
    #instance.set_courant_friedrichs_lewy_number(0.3)
    instance.commit_parameters()
    
    grid_in_memory = instance.grid.copy()
    grid_in_memory.rho  = 1e-2 | units.MSun / units.AU**3
    
    channel = grid_in_memory.new_channel_to(instance.grid)
    channel.copy()
    #print instance.grid.energy[0]
    #instance.energy = 0.001 | units.m**-1 * units.s**-2 * units.kg

    instance.initialize_grid()
    particle = Particle(
        mass = mass,
        position = length * [0.5, 0.5, 0.5],
        velocity = [0.0, 0.0, 0.0] | units.kms
    )
    
    gravity = Hermite(converter)
    gravity.particles.add_particle(particle)
    
    
    potential = gravity.get_potential_at_point(
        0 * instance.acceleration_grid.x.flatten(),
        instance.acceleration_grid.x.flatten(),
        instance.acceleration_grid.y.flatten(),
        instance.acceleration_grid.z.flatten()
    )
    fx,fy,fz = gravity.get_gravity_at_point(
        0 * instance.acceleration_grid.x.flatten(),
        instance.acceleration_grid.x.flatten(),
        instance.acceleration_grid.y.flatten(),
        instance.acceleration_grid.z.flatten()
    )
    print fx.min()
    print fx.max()
    #potential =  potential.reshape(instance.potential_grid.shape)
    #print potential.max()
    #print potential.min()
    
    grid_in_memory = instance.acceleration_grid.copy()
    grid_in_memory.fx = fx.reshape(instance.acceleration_grid.shape) 
    grid_in_memory.fy = fy.reshape(instance.acceleration_grid.shape) 
    grid_in_memory.fz = fz.reshape(instance.acceleration_grid.shape) 
    channel = grid_in_memory.new_channel_to(instance.acceleration_grid)
    channel.copy()
    
    instance.evolve_model(1e-8 | units.s)
    
    plot_grid(instance.grid[2])
    
def plot_grid(grid):
    rho = grid.rho[...,...,0].value_in(units.MSun / units.AU**3)
    print rho
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(rho, origin = 'lower')
    #figure.savefig('orszag_tang.png')
    pyplot.show()
    
    
if __name__ == '__main__':
    hydro_grid_in_potential_well()
        
