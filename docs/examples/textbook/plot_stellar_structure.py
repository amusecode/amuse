from matplotlib import pyplot
import numpy
from amuse.lab import *
from amuse.plot import plot, scatter

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def stellar_density_profile_at_time(mass, time):
    stellar = MESA()
    stellar.parameters.metallicity = 0.02
    star = stellar.particles.add_particle(Particle(mass=mass))

    stellar.evolve_model(time)

    rho = star.get_density_profile(star.get_number_of_zones())
    L = star.get_luminosity_profile(star.get_number_of_zones())
    m = star.get_cumulative_mass_profile(star.get_number_of_zones()) * mass
    stellar.stop()

    return m, L, rho

if __name__ == "__main__":
    
    x_label = "L [LSun]"
    y_label = "density [$g/cm^{3}$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True,
                          xsize=12, ysize=10)

    mass = [1, 1, 5, 5, 10, 10] | units.MSun
    age = [1, 10000, 1, 90, 1, 20] | units.Myr
    c = get_distinct(3)
    color = []
    ls = []
    for ci in c:
        color.append(ci)
        color.append(ci)
        ls.append("-")
        ls.append("--")
    symbol = ['v', 'v', 'o', 'o', '^', '^']
    i = 0
    for imass in range(len(mass)):
        dm = 0.2*mass[imass]
        time = age[imass]
        m, L, rho = stellar_density_profile_at_time(mass[imass], time)
        pyplot.plot(L.value_in(units.LSun),
                    rho.value_in(units.g/units.cm**3),
                    lw=4, c=color[imass], label='$10M_\odot$', ls=ls[imass])
        mlim = dm
        for mi in range(len(m)):
            if m[mi]>mlim:
                mlim += dm
                print mi, len(L), len(rho)
                pyplot.scatter(L[mi].value_in(units.LSun),
                               rho[mi].value_in(units.g/units.cm**3),
                               c=color[imass], s=150, marker=symbol[imass],
                               lw=0)
        pyplot.scatter(L[-1].value_in(units.LSun),
                       rho[-1].value_in(units.g/units.cm**3),
                       c=color[imass], s=150, marker=symbol[imass], lw=0)

    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    save_file = 'fig_1_5_10_MSun_stellar_core_luminosity'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
