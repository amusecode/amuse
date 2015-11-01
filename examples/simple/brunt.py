from matplotlib import pyplot
from amuse.plot import semilogy, xlabel, ylabel

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse.datamodel import Particle

def brunt_vaisala_frequency_squared_profile(mass, age):
    stellar_evolution = MESA()
    star = stellar_evolution.particles.add_particle(Particle(mass=mass))
    star.evolve_for(age)
    radius_profile = star.get_radius_profile()
    brunt_profile = star.get_brunt_vaisala_frequency_squared_profile()
    stellar_evolution.stop()
    return radius_profile.as_quantity_in(units.RSun), brunt_profile

def make_plot(radius_profile, brunt_profile, mass, age):
    figure = pyplot.figure()
    semilogy(radius_profile, -brunt_profile, 'g-', label = r'convective, $N^2$ < 0')
    semilogy(radius_profile, brunt_profile, 'r-', label = r'radiative, $N^2$ > 0')
    xlabel('Radius')
    ylabel(r'$\|N^2\|$')
    pyplot.title('Brunt-Vaisala frequency squared of a {0} star at {1}'.format(mass, age))
    pyplot.legend(loc=3)
    pyplot.show()   
    
if __name__ == '__main__':
    mass = 1.0 | units.MSun
    age = 1.0 | units.Gyr
    
    radius_profile, brunt_profile = brunt_vaisala_frequency_squared_profile(mass, age)
    make_plot(radius_profile, brunt_profile, mass, age)
