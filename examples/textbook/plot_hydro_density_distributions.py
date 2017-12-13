import numpy
from matplotlib import pyplot 
from amuse.lab import *
from amuse.plot import sph_particles_plot

from prepare_figure import figure_frame
from distinct_colours import get_distinct

###BOOKLISTSTART1###
def stellar_model(N, M, t=0.0|units.Myr):
    star = Particle(mass=M)
    stellar_evolution = EVtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    print "Evolving", star.mass, "to t=", t.in_(units.Myr)
    stellar_evolution.evolve_model(t)
    print "Stellar type:", stellar_evolution.particles.stellar_type.number
    print "Creating SPH particles from the (1D) stellar evolution model"
    sph_particles = convert_stellar_model_to_SPH(se_star, N).gas_particles
    stellar_evolution.stop()
    return sph_particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def plot_ZAMS_stellar_model(N, M):
    sph_particles = stellar_model(N, M)
    x_label = "x [R$_\odot$]"
    y_label = "y [R$_\odot$]"
    figure, ax = figure_frame(x_label, y_label, xsize=12, ysize=12)
    sph_particles_plot(sph_particles, min_size = 500,
                       max_size = 500, alpha=0.01,
                       view=(-2, 2, -2, 2)|units.RSun)
    ax.set_facecolor('white')
    
    save_file = 'stellar_2MSunZAMS_projected.png'
    pyplot.savefig(save_file)
    print 'Saved figure in file', save_file, '\n'
    pyplot.show()
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
from amuse.ext.molecular_cloud import molecular_cloud
def GMC_model(N, M, R):
    converter = nbody_system.nbody_to_si(M, R)
    sph_particles = molecular_cloud(targetN=N, convert_nbody=converter).result
    sph = Fi(converter)
    sph.gas_particles.add_particle(sph_particles)
    sph.evolve_model(1|units.day)
    ch = sph.gas_particles.new_channel_to(sph_particles)
    ch.copy()
    sph.stop()
    return sph_particles
###BOOKLISTSTOP3###

def plot_GMC_model(N, M, R):
    sph_particles = GMC_model(N, M, R)
    x_label = "x [pc]"
    y_label = "y [pc]"
    figure, ax = figure_frame(x_label, y_label, xsize=12, ysize=12)
    sph_particles_plot(sph_particles, min_size = 500,
                       max_size = 500, alpha=0.01,
                       view=(-20, 20, -20, 20)|units.parsec)
    ax.set_facecolor('white')

    file = 'molecular_cloud_projected.png'
    pyplot.savefig(file)
    print 'Saved figure in file', file


if __name__ in ("__main__","__plot__"):
    N = 10000
    M = 10000. | units.MSun
    R = 10. | units.parsec

    M = 2 | units.MSun
    print ''
    plot_ZAMS_stellar_model(N, M)
    plot_GMC_model(N, M, R)
    print ''
