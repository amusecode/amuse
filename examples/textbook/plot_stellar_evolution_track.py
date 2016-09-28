from matplotlib import pyplot
import numpy
from amuse.lab import *

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def get_stellar_track(mass):

    stars = Particles(1)
    stars[0].mass = mass
#    stellar = SeBa()
    stellar = SSE()
    stellar.particles.add_particles(stars)
    attributes = ["temperature", "luminosity", "stellar_type"]
    to_framework = stellar.particles.new_channel_to(stars,
                                                    attributes=attributes,
                                                    target_names=attributes)
    t = [] | units.Myr
    T = [] | units.K
    L = [] | units.LSun
    stp = [] | units.stellar_type
    #dt = 0.1|units.Myr
    Helium_White_Dwarf = 10 | units.stellar_type
    while stellar:
        #dt = stellar.stars[0].time_step
        #stellar.evolve_model(stellar.model_time + dt)
        stellar.evolve_model()
        to_framework.copy()
        if stars[0].stellar_type>=Helium_White_Dwarf:
            stellar.stop()
            stellar=False
        else:
            t.append(stellar.model_time)
            T.append(stars[0].temperature)
            L.append(stars[0].luminosity)
            stp.append(stars[0].stellar_type)
            
    return t, T, L, stp
        
if __name__ == "__main__":

    x_label = "T [$K$]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True, xsize=12, ysize=10)
#    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=10)
    pyplot.xlim(1e+6, 1e+3)
    pyplot.ylim(1.e-1, 3e+5)

    color = get_distinct(10)

    t, T, L, stp = get_stellar_track(1|units.MSun)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=4, c=color[stp[i].value_in(units.stellar_type)], label='$10M_\odot$')

    t, T, L, stp = get_stellar_track(5|units.MSun)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=6, c=color[stp[i].value_in(units.stellar_type)], label='$10M_\odot$')

    t, T, L, stp = get_stellar_track(20|units.MSun)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=4, c=color[stp[i].value_in(units.stellar_type)], label='$10M_\odot$')
        

#    pyplot.text(20., 1.1e+7,'$1M_\odot$')
#    pyplot.text(100, 6.8e+7,'$10M_\odot$')
#    pyplot.text(48000., 2.4e+9,'$100M_\odot$')

#    pyplot.show()
    pyplot.savefig("fig_stellar_evolution_track")
