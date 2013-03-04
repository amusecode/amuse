"""
Make plots of the stellar structure of a star of given mass and age
"""

import numpy
from matplotlib import pyplot
from amuse.plot import plot, semilogx, semilogy, loglog, xlabel, ylabel

from amuse.units import units
from amuse.community.evtwin.interface import EVtwin
from amuse.datamodel import Particle

def structure_from_star(mass, age):
    stellar_evolution = EVtwin()
    star = stellar_evolution.particles.add_particle(Particle(mass=mass))
    star.evolve_for(age)
    
    radius_profile = star.get_radius_profile()
    density_profile = star.get_density_profile()
    if hasattr(star, "get_mass_profile"):
        mass_profile = star.get_mass_profile()* star.mass
    else:
        radii_cubed = radius_profile**3
        radii_cubed.prepend(0|units.m**3)
        mass_profile = (4.0/3.0 * numpy.pi) * density_profile * (radii_cubed[1:] - radii_cubed[:-1])
        print "Derived mass profile from density and radius."
    
    return dict(
        radius = radius_profile.as_quantity_in(units.RSun),
        density = density_profile,
        mass = mass_profile,
        temperature = star.get_temperature_profile(),
        pressure = star.get_pressure_profile(),
        composition = star.get_chemical_abundance_profiles()
    )

def temperature_density_plot(data, mass, age):
    figure = pyplot.figure(figsize = (8, 10))
    pyplot.subplot(2, 1, 1)
    ax = pyplot.gca()
    plotT = semilogy(data["radius"], data["temperature"], 'r-', label = r'$T(r)$')
    xlabel('Radius')
    ylabel('Temperature')
    ax.twinx()
    plotrho = semilogy(data["radius"], data["density"], 'g-', label = r'$\rho(r)$')
    plots = plotT + plotrho
    labels = [one_plot.get_label() for one_plot in plots]
    ax.legend(plots, labels, loc=3)
    ylabel('Density')
    pyplot.subplot(2, 1, 2)
    semilogy(data["radius"], data["composition"][0], label = "H")
    semilogy(data["radius"], data["composition"][1] + data["composition"][2], label = "He")
    semilogy(data["radius"], data["composition"][3], label = "C")
    semilogy(data["radius"], data["composition"][4], label = "N")
    semilogy(data["radius"], data["composition"][5], label = "O")
    pyplot.ylim(0.0, 1.0)
    xlabel('Radius')
    ylabel('Mass fraction')
    pyplot.legend()
    pyplot.suptitle('Structure of a {0} star at {1}'.format(mass, age))
    pyplot.show()   
    
if __name__ in ('__main__', '__plot__'):
    mass = 1.0 | units.MSun
    age = 5.0 | units.Gyr
    
    data = structure_from_star(mass, age)
    temperature_density_plot(data, mass, age)
    
