import numpy
from matplotlib import pyplot 
from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.plot import sph_particles_plot, native_plot

from prepare_figure import figure_frame
from distinct_colours import get_distinct

def create_molecular_cloud(N, Mcloud, Rcloud, t_end):

    converter = nbody_system.nbody_to_si(Mcloud,Rcloud)
    parts = new_plummer_gas_model(N, convert_nbody=converter)

    sph=Fi(converter)
    sph.gas_particles.add_particle(parts)

    sph.evolve_model(t_end)
    ch = sph.gas_particles.new_channel_to(parts)
    ch.copy()
    sph.stop()
    return parts

def convert_star_to_hydro_model(M, t_end):
    star =  Particle(mass=M)
    
    stellar_evolution = EVtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    
    print "Evolving", star.mass, "star with", stellar_evolution.__class__.__name__, "up to", t_end.in_(units.Myr)
    stellar_evolution.evolve_model(t_end)
    
    print "Creating SPH particles from the (1D) stellar evolution model"
    sph_particles = convert_stellar_model_to_SPH(
        se_star, 
        10000
    ).gas_particles
    stellar_evolution.stop()
    return sph_particles

def stellar_model(N, M, t=0.0|units.Myr):
    star =  Particle(mass=M)
    stellar_evolution = EVtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    print "Evolving", star.mass, "star with", stellar_evolution.__class__.__name__, "to t=", t.in_(units.Myr)
    stellar_evolution.evolve_model(t)
    print "Stellar type:", stellar_evolution.particles.stellar_type
    
    print "Creating SPH particles from the (1D) stellar evolution model"
    sph_particles = convert_stellar_model_to_SPH(se_star, N).gas_particles
    stellar_evolution.stop()
    return sph_particles

def plot_ZAMS_stellar_model(N, M):
    sph_particles = stellar_model(N, M)
    figure = pyplot.figure(figsize=(12, 12))
    plot = figure.add_subplot(1,1,1)
    pyplot.rcParams.update({'font.size': 30})
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.minorticks_on() # switch on the minor ticks
    lim = 2
    sph_particles_plot(sph_particles, min_size = 500, max_size = 500, alpha=0.01,
                       view=(-lim, lim, -lim, lim)|units.RSun)
    ax.set_facecolor('white')
    pyplot.xlabel("x [R$_\odot$]")
    pyplot.ylabel("y [R$_\odot$]")
    #native_plot.savefig("stellar_2MSunZAMS_projected")
    pyplot.savefig("stellar_2MSunZAMS_projected")

def plot_stellar_model(N, M, t):
    sph_particles = stellar_model(N, M, t)
    x_label = "x [R$_\odot$]"
    y_label = "y [R$_\odot$]"
    figure, ax = figure_frame(x_label, y_label, xsize=12, ysize=12)
    sph_particles_plot(sph_particles, min_size = 500, max_size = 500, alpha=0.01, view=(-5, 5, -5, 5)|units.RSun)
#    sph_particles_plot(sph_particles, min_size = 500, max_size = 500, alpha=0.01, view=(-2, 2, -2, 2)|units.RSun)
    ax.set_facecolor('white')
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    native_plot.savefig("stellar_2MSun_projected")
    #pyplot.savefig("stellar_2MSun_projected")

def plummer_model(N, M, R):
    converter = nbody_system.nbody_to_si(M, R)
    parts = new_plummer_gas_model(N, convert_nbody=converter)
    sph=Fi(converter)
    sph.gas_particles.add_particle(parts)
    sph.evolve_model(1|units.day)
    ch = sph.gas_particles.new_channel_to(parts)
    ch.copy()
    sph.stop()
    return parts

def GMC_model(N, M, R):
    converter = nbody_system.nbody_to_si(M, R)
    gmc=molecular_cloud(targetN=N, convert_nbody=converter)
    parts= gmc.result
    sph=Fi(converter)
    sph.gas_particles.add_particle(parts)
    sph.evolve_model(1|units.day)
    ch = sph.gas_particles.new_channel_to(parts)
    ch.copy()
    sph.stop()
    return parts

def plot_plummer_model(N, M, R):
    sph_particles = plummer_model(N, M, R)
    x_label = "x [pc]"
    y_label = "y [pc]"
    figure, ax = figure_frame(x_label, y_label, xsize=12, ysize=12)
    sph_particles_plot(sph_particles, min_size = 500, max_size = 500, alpha=0.01, view=(-50, 50, -50, 50)|units.parsec)
    ax.set_facecolor('white')
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    native_plot.savefig("plummer_projected")
    
def plot_GMC_model(N, M, R):
    sph_particles = GMC_model(N, M, R)

    x_label = "x [pc]"
    y_label = "y [pc]"
    figure, ax = figure_frame(x_label, y_label, xsize=12, ysize=12)
    lim = 20
    sph_particles_plot(sph_particles, min_size = 500, max_size = 500, alpha=0.01,
                       view=(-lim, lim, -lim, lim)|units.parsec)
    ax.set_facecolor('white')
    #native_plot.savefig("molecular_cloud_projected")
    pyplot.savefig("molecular_cloud_projected")


if __name__ in ("__main__","__plot__"):
    N = 10000
    M = 10000. | units.MSun
    R = 10. | units.parsec
#    plot_plummer_model(N, M, R)
#    plot_GMC_model(N, M, R)

    M = 2 | units.MSun
    plot_ZAMS_stellar_model(N, M)
#    plot_stellar_model(N, M, t=1.3|units.Gyr)
