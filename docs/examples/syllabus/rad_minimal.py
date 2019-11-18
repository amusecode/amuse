"""
   Minimalistic routine for running a radiative transfer code.
"""
from amuse.lab import *

def main(N=1000, Lstar=100|units.LSun, boxsize=10|units.parsec, 
         rho=1.0 | (units.amu/units.cm**3), t_end=0.1 |units.Myr):

    internal_energy = (9. |units.kms)**2

    source=Particle()
    source.position = (0, 0, 0) |units.parsec
    source.flux = Lstar/(20. | units.eV)
    source.rho = rho
    source.xion = 0.0
    source.u = internal_energy

    from amuse.ext.molecular_cloud import ism_cube
    ism = ism_cube(N, boxsize/2., rho, internal_energy).result
    ism.rho = rho
    ism.flux = 0. | units.s**-1
    ism.xion = source.xion

    radiative = SimpleX()
    radiative.parameters.box_size=1.001*boxsize    
    radiative.parameters.timestep=0.001 | units.Myr

    radiative.particles.add_particle(source)
    radiative.particles.add_particles(ism)

    radiative.evolve_model(t_end)
    print "min ionization:", radiative.particles.xion.min()
    print "average Xion:", radiative.particles.xion.mean()
    print "max ionization:", radiative.particles.xion.max()
    radiative.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 1000,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", default = 0.1|units.Myr,
                      help="radiation time [%default]")
    result.add_option("-L", unit=units.LSun,
                      dest="Lstar", default = 100|units.LSun,
                      help="luminosity of ionizing source [%default]")
    result.add_option("-p", unit=units.amu/units.cm**3, 
                      dest="rho", default = 1|units.amu/units.cm**3,
                      help="interstellar density [%default] amu/cm^3")
    result.add_option("-d", unit=units.parsec,
                      dest="boxsize", default = 100|units.parsec,
                      help="size of the density box [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

