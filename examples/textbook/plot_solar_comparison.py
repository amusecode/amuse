import numpy
from amuse.lab import *
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from amuse.plot import scatter

from prepare_figure import single_frame
from distinct_colours import get_distinct

from amuse.community.seba.interface import SeBa
    
def main(t_end, mass, z, Tstar, Lstar):

    stellar_evolution_codes = [SeBa(), SSE(), MESA(), EVtwin()]
    label = ["SeBa", "SSE", "MESA", "EVtwin"]
    marker = ["o", "v", "<", ">"]

    x_label = "$(T-T_\odot)/T_\odot)$"
    y_label = "$(L-L_\odot)/L_\odot)$"
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=10)
    pyplot.xlim(-0.006, 0.004)
    pyplot.ylim(-0.1, 0.1)    
    color = get_distinct(6)
    pyplot.scatter([0], [0], marker="o", c=color[3], label="Sun", s=200, lw=0)

    for si in range(len(stellar_evolution_codes)):
        stellar = stellar_evolution_codes[si]
        stellar.parameters.metallicity = z

        star = Particles(1)
        star.mass = mass
        t_end = 6000.0 | units.Myr
        stellar.particles.add_particles(star)
        attributes = ["temperature", "luminosity","age"]
        to_framework = stellar.particles.new_channel_to(star,
                                                    attributes=attributes,
                                                    target_names=attributes)
        t = [] | units.Myr
        T = []
        L = []
        min_dist_sun = 10000.0
        current_time = 3000.0 | units.Myr

        print(label[si])
        #dt = 50 | units.Myr
        #time = 4000 | units.Myr
        while stellar:
            print(stellar.model_time.value_in(units.Myr)) 
            current_time = current_time + stellar.particles[0].time_step
            stellar.evolve_model(current_time)

            to_framework.copy()

            if star[0].age >= t_end:
                stellar.stop()
                stellar = False
            else:
                L.append((star[0].luminosity - Lstar)/Lstar)
                T.append((star[0].temperature - Tstar)/Tstar)
                t.append(star[0].age)

                deltaL = numpy.abs((star[0].luminosity - Lstar)/Lstar)
                deltaT = numpy.abs((star[0].temperature - Tstar)/Tstar)       

                dist = numpy.sqrt(deltaL*deltaL + deltaT*deltaT)
                
                if min_dist_sun > dist:

                    min_dist_sun = dist

                    L_sim_sun = (star[0].luminosity - Lstar)/Lstar
                    T_sim_sun = (star[0].temperature - Tstar)/Tstar

                    eta = star[0].age
        print(eta)
        if si==3: 
            pyplot.plot(T, L,ls='-', marker=marker[si], color=color[5], markersize=10)
            pyplot.scatter(T_sim_sun, L_sim_sun, marker=marker[si],
                           color=color[5], label=label[si], s=300, lw=1)
        else:
            pyplot.plot(T, L,ls='-', marker=marker[si], color=color[si], markersize=10)
            pyplot.scatter(T_sim_sun, L_sim_sun, marker=marker[si],
                           color=color[si], label=label[si], s=300, lw=1)

    pyplot.legend(scatterpoints=1, loc='best')

    save_file = 'fig_SunComparison.png'
    pyplot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    pyplot.show()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-T", unit=units.K,
                      dest="Tstar", type="float",default = 5778 |units.K,
                      help="stellar temperature [%defailt]")
    result.add_option("-L", unit=units.LSun,
                      dest="Lstar", type="float",default = 1 |units.LSun,
                      help="stellar luminosity [%defailt]")
    result.add_option("-m", unit=units.MSun,
                      dest="mass", type="float",default = 1.0 |units.MSun,
                      help="stellar mass [%defailt]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 4.6 |units.Gyr,
                      help="end time of the simulation [%defailt]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%defailt]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

