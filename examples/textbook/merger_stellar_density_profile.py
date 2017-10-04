"""
   Script to initialize a star and print its structure
"""
import numpy
from amuse.lab import *
from matplotlib import pyplot
from amuse.plot import plot, xlabel, ylabel

from prepare_figure import figure_frame
from distinct_colours import get_distinct

def merge_two_stars(Mprim, Msec, t_coll):
    bodies = Particles(mass=[Mprim.value_in(units.MSun),
                             Msec.value_in(units.MSun)] |units.MSun)
        
    stellar = MESA()
    primary = stellar.particles.add_particles(bodies[0].as_set())
    secondary = stellar.particles.add_particles(bodies[1].as_set())

    stellar.evolve_model(t_coll)

    print "Pre merger:\n", stellar.particles
    stellar.merge_colliding(primary.copy(), secondary.copy(),
                            MakeMeAMassiveStar,
                            dict(), dict(target_n_shells_mixing = 2000),
                            return_merge_products=["se"])
    print "Post merger:\n", stellar.particles

    radius = stellar.particles[0].get_radius_profile()
    rho    = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho

def get_density_profile(code=MESA, M=1.0|units.MSun, z=0.02):
    stellar = code()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=M))
    print "Nzones=", stellar.particles.get_number_of_zones()
    radius = stellar.particles[0].get_radius_profile()
    rho    = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho

def main(M, z, output_filename):
    numpy.random.seed(31415)
    x_label = "$R$ [R$_\odot$]"
    y_label = "$\\rho$ [g/cm$^{3}$]"
    fig, ax = figure_frame(x_label, y_label, xsize=12, ysize=8)
    cols = get_distinct(3)

    r, rho = get_density_profile(EVtwin, M, z)
    pyplot.plot(r.value_in(units.RSun),
                rho.value_in(units.g/units.cm**3),
                label="EVtwin", c=cols[0])
    r, rho = get_density_profile(MESA, M, z)
    pyplot.plot(r.value_in(units.RSun),
                rho.value_in(units.g/units.cm**3),
                label="MESA", c=cols[1])

    # Run the merger code.
    
    r, rho = merge_two_stars(0.5*M, 0.5*M, 1|units.yr)
    
    pyplot.plot(r.value_in(units.RSun),
                rho.value_in(units.g/units.cm**3),
                label="MESA", c=cols[2])
    pyplot.semilogy()
    
    if output_filename is not None:
        pyplot.savefig(output_filename)
        print '\nSaved figure in file', output_filename,'\n'
    else:
        output_filename = 'merger_stellar_density_profile.png'
        pyplot.savefig(output_filename)
        print '\nSaved figure in file', output_filename,'\n'
        pyplot.show()
   
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit= units.MSun,
                      dest="M", type="float",default = 2.0 | units.MSun,
                      help="stellar mass [1.0] %unit")
    result.add_option("-o", 
                      dest="output_filename",
                      default=None,
                      help="output filename [%default]")
    result.add_option("-z", dest="z", type="float", default=0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
