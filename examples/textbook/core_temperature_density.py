from amuse.lab import *

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
HeWhiteDwarf = 10 | units.stellar_type
    
def stellar_remnant(stellar):
    remnant = True
    if stellar.particles[0].stellar_type<HeWhiteDwarf or stellar.particles[0].stellar_type>11|units.stellar_type:
        remnant = False
    return remnant

def stellar_core_temperature_and_density(M, z):
    stellar = MESA()
    stellar.parameters.metallicity = z
    star = stellar.particles.add_particle(Particle(mass=M))

    while not stellar_remnant(stellar):
        stellar.evolve_model()
        T_core = star.get_temperature_profile(star.get_number_of_zones())[0]
        density_core = star.get_density_profile(star.get_number_of_zones())[0]

        T_surface = star.get_temperature_profile(star.get_number_of_zones())[-1]
        density_surface = star.get_density_profile(star.get_number_of_zones())[-1]

        print(star.age, T_surface, density_surface, T_core, density_core)

    stellar.stop()
    
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
    stellar_core_temperature_and_density(**o.__dict__)
