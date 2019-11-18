"""
   Evolve a population of N stars.
   initial mass function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import sys
import numpy
from amuse.lab import *
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel

from amuse.community.seba.interface import SeBa
    
def main(N=10, t_end=10|units.Myr, dt=1|units.Myr, filename="stellar.hdf5", 
         Mmin=1.0|units.MSun, Mmax= 100|units.MSun, z=0.02, C="SeBa"):

    if C.find("SeBa")>=0:
        print "SeBa"
        stellar = SeBa()
    elif C.find("SSE")>=0:
        print "SSE"
        stellar = SSE()
    elif C.find("MESA")>=0:
        stellar = MESA()
    elif C.find("EVtwin")>=0:
        stellar = EVtwin()
    else:
        print "No stellar model specified."
        return
    stellar.parameters.metallicity = z

    mZAMS = new_salpeter_mass_distribution(N, Mmin, Mmax)
    mZAMS = mZAMS.sorted()
    bodies = Particles(mass=mZAMS)
    stellar.particles.add_particles(bodies)
    stellar.commit_particles()

    write_set_to_file(stellar.particles, filename, 'hdf5')

    Mtot_init = stellar.particles.mass.sum()
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt

        stellar.evolve_model(time)
        write_set_to_file(stellar.particles, filename, 'hdf5')

        Mtot = stellar.particles.mass.sum()

        print "T=", time, 
        print "M=", Mtot, "dM[SE]=", Mtot/Mtot_init #, stellar.particles[0].stellar_type

    T = [] 
    L = [] 
    r = []
    import math
    for si in stellar.particles:
        T.append(math.log10(si.temperature.value_in(units.K)))
        L.append(math.log10(si.luminosity.value_in(units.LSun)))
        r.append(1.0+10*si.radius.value_in(units.RSun))
    stellar.stop()
    pyplot.xlim(4.5, 3.3)
    pyplot.ylim(-4.0, 3.0)
    scatter(T, L, s=r)
    pyplot.xlabel("$log_{10}(T/K)$")
    pyplot.ylabel("$log_{10}(L/L_\odot)$")
    pyplot.show()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-C", dest="C", default = "SeBa",
                      help="stellar evolution code [SeBa]")
    result.add_option("-d", unit=units.Myr,
                      dest="dt", type="float", default = 100.0 |units.Myr,
                      help="diagnostics time step [%defailt]")
    result.add_option("-f", dest="filename", default = "stellar.hdf5",
                      help="output filename [stellar.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [%defailt]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mmax", type="float",default = 100|units.MSun,
                      help="maximal stellar mass [%defailt]")
    result.add_option("-m", unit=units.MSun,
                      dest="Mmin", type="float",default = 1.0 |units.MSun,
                      help="minimal stellar mass [%defailt]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 100.0 |units.Myr,
                      help="end time of the simulation [%defailt]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%defailt]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

