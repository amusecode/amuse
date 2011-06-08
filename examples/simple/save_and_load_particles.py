from amuse.lab import *
from amuse.plot import *

number_of_particles = 1000
end_time = 1000 | units.yr

nbody_converter = nbody_system.nbody_to_si(number_of_particles | units.MSun, 1 | units.parsec)


def evolve_cluster_and_save_particles_halfway(particles, end_time):
    dynamics_code = Hermite(nbody_converter)
    dynamics_code.particles.add_particles(particles)
    print "Start evolve"
    dynamics_code.evolve_model(0.5 * end_time)
    print "Writing particles to file at time", (0.5 * end_time).as_quantity_in(units.yr)
    write_set_to_file(dynamics_code.particles, "save_and_load_particles.csv","csv")
    print "Continue evolution of original particles to time", end_time
    dynamics_code.evolve_model(end_time)
    positions_evolved_particles = dynamics_code.particles.position.as_quantity_in(units.parsec)
    dynamics_code.stop()
    return positions_evolved_particles

def evolve_saved_particles(end_time):
    loaded_particles = read_set_from_file("save_and_load_particles.csv","csv")
    dynamics_code = Hermite(nbody_converter)
    dynamics_code.parameters.time = 0.5 * end_time
    dynamics_code.particles.add_particles(loaded_particles)
    print "Continue evolution of saved particles to time", end_time
    dynamics_code.evolve_model(end_time)
    positions_evolved_loaded_particles = dynamics_code.particles.position.as_quantity_in(units.parsec)
    dynamics_code.stop()
    return positions_evolved_loaded_particles

def plot_compare_original_and_loaded(original, loaded):
    print "Making a plot to compare the results"
    figure = native_plot.figure(figsize = (6, 10))
    subplot1 = figure.add_subplot(2, 1, 1, aspect='equal', adjustable = 'datalim')
    scatter(original.x, original.y)
    ylabel(r'y$_\mathrm{original}$')
    
    subplot2 = figure.add_subplot(2, 1, 2, sharex = subplot1)
    scatter(original.x, loaded.x - original.x)
    xlabel(r'x$_{\mathrm{original}}$')
    ylabel(r'x$_{\mathrm{loaded}}$ - x$_{\mathrm{original}}$')
    native_plot.subplots_adjust(hspace = 0.05, left = 0.2)
    native_plot.setp(subplot1.get_xticklabels(), visible=False)
    native_plot.show()
    
if __name__ in ('__main__', '__plot__'):
    print "Creating a Plummer model consisting of", number_of_particles, "equal-mass particles."
    particles = new_plummer_sphere(number_of_particles, nbody_converter)
    particles_after_evolve = evolve_cluster_and_save_particles_halfway(particles, end_time)
    saved_particles_after_evolve = evolve_saved_particles(end_time)
    plot_compare_original_and_loaded(particles_after_evolve, saved_particles_after_evolve)
