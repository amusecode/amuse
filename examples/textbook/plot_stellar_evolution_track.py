from matplotlib import pyplot
import numpy
from amuse.lab import *

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def get_stellar_track(mass):

    stars = Particles(1)
    stars[0].mass = mass
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

sigma = constants.Stefan_hyphen_Boltzmann_constant
def stellar_radius(L, T):
    return numpy.sqrt(L/(4*numpy.pi*sigma*T**4))

def stellar_luminosity(R, T):
    return R**2 * 4*numpy.pi*sigma*T**4

def get_color_based_on_stellar_type(istp):
    color = get_distinct(5)
    if istp.value_in(units.stellar_type)<2:
        c = color[0]
    elif istp.value_in(units.stellar_type)<6:
        c = color[3]
    else:
        c = color[4]
    return c

def stellar_temperature(R, L):
    return (L/(4*numpy.pi*sigma*R**2))**(1./4.)

if __name__ == "__main__":

    x_label = "T [$K$]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True, xsize=12, ysize=10)
    pyplot.xlim(1e+6, 1e+3)
    pyplot.ylim(1.e-1, 3e+5)

    L = [1, 1.e+5] | units.LSun
    T = stellar_temperature(0.01|units.RSun, L)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    pyplot.text(9e+5, 6.e+4, '$0.1 R_\odot$', rotation=-60)
    L = [1, 1.e+5] | units.LSun
    T = stellar_temperature(1|units.RSun, L)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    pyplot.text(1.8e+4, 80, '$1 R_\odot$', rotation=-60)
    L = [1, 1.e+5] | units.LSun
    T = stellar_temperature(100|units.RSun, L)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    pyplot.text(3.e+3, 600, '$100 R_\odot$', rotation=-60)
    L = [1, 1.e+5] | units.LSun
    T = stellar_temperature(10000|units.RSun, L)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    
    t, T, L, stp = get_stellar_track(1|units.MSun)
    #pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        c = get_color_based_on_stellar_type(stp[i])
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=4, c=c, label='$1M_\odot$')
#        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
#                    lw=4, c=color[stp[i].value_in(units.stellar_type)], label='$1M_\odot$')

    t, T, L, stp = get_stellar_track(5|units.MSun)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        c = get_color_based_on_stellar_type(stp[i])
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=6, c=c, label='$5M_\odot$')
#        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
#                    lw=6, c=color[stp[i].value_in(units.stellar_type)], label='$5M_\odot$')

    t, T, L, stp = get_stellar_track(20|units.MSun)
    pyplot.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c='k')
    for i in range(len(T)-2):
        c = get_color_based_on_stellar_type(stp[i])
        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
                    lw=4, c=c, label='$20M_\odot$')
#        pyplot.plot(T[i:i+2].value_in(units.K), L[i:i+2].value_in(units.LSun),
#                    lw=4, c=color[stp[i].value_in(units.stellar_type)], label='$20M_\odot$')


#    pyplot.legend(loc="lower right")
    
    Ttext = 700000
    Ltext = 0.1
#    stellar_type = ["", "main sequence", "Hertzsprung gap", "first giant branch", "core helium burning", "first asymptotic giant brach", "second AGB"]
    stellar_type = ["", "main sequence", "giant", "remnant"]
    dL = 1.8
    i = 0
    for sti in range(len(stellar_type)):
        if sti<=1:
            c = get_color_based_on_stellar_type(1|units.stellar_type)
        elif sti == 2:
            c = get_color_based_on_stellar_type(5|units.stellar_type)
        else:
            c = get_color_based_on_stellar_type(10|units.stellar_type)
        pyplot.text(Ttext, dL**sti*Ltext, stellar_type[sti], color=c, fontsize=24)
        
    #pyplot.show()
    pyplot.savefig("fig_stellar_evolution_track")

