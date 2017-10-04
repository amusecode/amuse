"""
   Evolve a population of N stars.
   Initial mass function between Mmin and Mmax with stellar evolution
   for metallicity z.
"""
import sys
import numpy
from amuse.lab import *
from matplotlib import pyplot

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def _get_stellar_temperature_and_luminosity(stars, C, z=0.02,
                                            t_end=100|units.Myr, write=False):

    if C.find("SeBa") >= 0 :
        stellar = SeBa()
        filename = "Stellar_SeBa.h5"
    if C.find("SSE") >= 0 :
        stellar = SSE()
        filename = "Stellar_SSE.h5"
    if C.find("MESA") >= 0 :
        stellar = MESA()
        filename = "Stellar_MESA.h5"
    if C.find("EVtwin") >= 0 :
        stellar = SeBa()
        filename = "Stellar_EVtwin.h5"

    stellar.parameters.metallicity = z
    stellar.particles.add_particles(stars)

    for si in stellar.particles:
        try:
            si.evolve_for(t_end)
        except:
            print "Failed to evolve star: m=", si.mass.in_(units.MSun)
    if write:
        write_set_to_file(stellar.particles, filename, 'hdf5')
    stellar.stop()

def get_stellar_temperature_and_luminosity(stars, C, z=0.02,
                                           t_end=100|units.Myr, write=False):

    if C.find("SeBa") >= 0 :
        filename = "Stellar_SeBa.h5"
        stellar = SeBa()
        stellar.parameters.metallicity = z
        stellar.particles.add_particles(stars)
        channel_to_framework = stellar.particles.new_channel_to(stars)
        stellar.evolve_model(t_end)
        channel_to_framework.copy_attributes(["radius", "temperature",
                                              "luminosity"])
        stellar.stop()
    else:
        if C.find("MESA") >= 0 :
            filename = "Stellar_MESA.h5"
        if C.find("EVtwin") >= 0 :
            filename = "Stellar_EVtwin.h5"

        for si in stars:
            stellar = MESA()
            stellar.parameters.metallicity = z
            stellar.particles.add_particle(si)
            channel_to_framework = stellar.particles.new_channel_to(stars)
            try:
                stellar.evolve_model(t_end)
                channel_to_framework.copy_attributes(["radius", "temperature",
                                                      "luminosity"])
                print "Successfully evolved star: m=", si.mass.in_(units.MSun)
            except:
                print "Failed to evolve star: m=", si.mass.in_(units.MSun)
            stellar.stop()
    if write:
        write_set_to_file(stars, filename, 'hdf5')
        
def plot_HRD(filename, color):
        stars = read_set_from_file(filename, 'hdf5')
        T = stars.temperature.value_in(units.K)
        L = stars.luminosity.value_in(units.LSun)
        R = stars.radius.value_in(units.RSun)
        
        R = 80*numpy.sqrt(R)
        pyplot.scatter(T, L, c=color, lw=0, s=R)

def main(N, t_end, z, C="SeBa", plot=False):
    if "SeBa" not in C:
        x_label = "T [K]"
        y_label = "L [L$_\odot$]"
        figure = single_frame(x_label, y_label, logx=True, logy=True,
                              xsize=14, ysize=10)
        color = get_distinct(4)
        pyplot.xlim(1.e+5, 1.e+3)
        pyplot.ylim(1.e-4, 1.e+4)
        filename = "Stellar_"+"SeBa"+".h5"
        plot_HRD(filename, color[0])
        filename = "Stellar_"+C+".h5"
        plot_HRD(filename, color[1])

        save_file = 'HRD_N3000at4500Myr.png'
        pyplot.savefig(save_file)
        print '\nSaved figure in file', save_file,'\n'
        
    elif not plot:
        numpy.random.seed(1)
        masses = new_salpeter_mass_distribution(N)
        stars = Particles(mass=masses)
        get_stellar_temperature_and_luminosity(stars, C=C, z=z,
                                               t_end=t_end, write=True)
    else:
        x_label = "T [K]"
        y_label = "L [L$_\odot$]"
        figure = single_frame(x_label, y_label, logx=True, logy=True,
                              xsize=14, ysize=10)
        color = get_distinct(4)
        pyplot.xlim(1.e+5, 1.e+3)
        pyplot.ylim(1.e-4, 1.e+4)
        filename = "Stellar_"+C+".h5"
        plot_HRD(filename, color[0])

        save_file = 'HRD_N3000at4500Myr.png'
        pyplot.savefig(save_file)
        print '\nSaved figure in file', save_file,'\n'
        pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-C", dest="C", default = "SeBa",
                      help="stellar evolution code [SeBa]")
    result.add_option("-N", dest="N", type="int",default = 3000,
                      help="number of stars [10]")
    result.add_option("-t", dest="t_end", unit=units.Myr,
                      type="float", default = 4500.0|units.Myr,
                      help="end time of the simulation [100] Myr")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    result.add_option("-p", dest="plot", action="store_true", default = False,
                      help="plot")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

                  
