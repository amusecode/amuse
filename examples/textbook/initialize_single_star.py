from amuse.lab import *
from matplotlib import pyplot
from amuse.plot import plot, xlabel, ylabel

###BOOKLISTSTART###
def get_density_profile(code=EVtwin, M=2.0|units.MSun, z=0.02):
    stellar = code()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=M))
    radius = stellar.particles[0].get_radius_profile()
    rho    = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho
###BOOKLISTSTOP###

def initialize_single_star(M, z):
    r, rho = get_density_profile(EVtwin, M, z)
    plot(r.in_(units.RSun), rho, label="EVtwin")
    r, rho = get_density_profile(MESA, M, z)
    plot(r.in_(units.RSun), rho, label="MESA")
    pyplot.xlabel("$R$ [$R_\odot$]")
    pyplot.ylabel("density [$g/cm^3$]")
    pyplot.semilogy()
    
    save_file = 'initialize_single_star.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit= units.MSun,
                      dest="M", type="float",default = 1.0 | units.MSun,
                      help="stellar mass [1.0] %unit")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    initialize_single_star(**o.__dict__)
