from amuse.lab import *

from matplotlib import pyplot

def hydro_grid_in_potential_well(mass = 1 | units.MSun, length = 1 | units.AU):
    converter = nbody_system.nbody_to_si(mass, length)
    
    instance=Athena(converter)
    instance.initialize_code()
    instance.parameters.nx = 20
    instance.parameters.ny = 20
    instance.parameters.nz = 1
    instance.parameters.length_x = length
    instance.parameters.length_y = length
    instance.parameters.length_z = length
    instance.parameters.x_boundary_conditions = ("outflow","outflow")
    instance.parameters.y_boundary_conditions = ("outflow","outflow")
    instance.parameters.z_boundary_conditions = ("outflow","outflow")
    
    instance.stopping_conditions.number_of_steps_detection.enable()
    instance.set_has_external_gravitational_potential(1)
    #instance.set_courant_friedrichs_lewy_number(0.3)
    instance.commit_parameters()
    
    grid_in_memory = instance.grid.copy()
    grid_in_memory.rho  = 1e-6 | units.MSun / units.AU**3
    
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
    print grid_in_memory.shape
    dx = (grid_in_memory.x[1][0][0] - grid_in_memory.x[0][0][0]).as_quantity_in(units.AU)
    print  dx**2 
    gravity.parameters.epsilon_squared = (10*dx)**2 
    gravity.particles.add_particle(particle)
    
    
    potential = gravity.get_potential_at_point(
        0 * instance.potential_grid.x.flatten(),
        instance.potential_grid.x.flatten(),
        instance.potential_grid.y.flatten(),
        instance.potential_grid.z.flatten()
    )
    
    potential =  potential.reshape(instance.potential_grid.x.shape)
    print potential.max()
    print potential.min()
    instance.potential_grid.potential =  potential
    
    instance.evolve_model(1 | units.s)
    
    rho = instance.grid.rho[...,...,0].value_in(units.MSun / units.AU**3)
    #rho = potential[...,...,0].value_in(potential.unit)
    plot_grid(rho)
    
def plot_grid(x):
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(x, origin = 'lower')
    #figure.savefig('orszag_tang.png')
    pyplot.show()
    
    
if __name__ == '__main__':
    hydro_grid_in_potential_well()
        
