from amuse.units import units
from amuse.datamodel import Particle
from amuse.community.mesa.interface import MESA
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel

def main():
    temperatures_original = [] | units.K
    luminosities_original = [] | units.LSun
    temperatures_helium = [] | units.K
    luminosities_helium = [] | units.LSun
    
    star =  Particle()
    star.mass = 10.0 | units.MSun
    stop_radius = 100 | units.RSun
    
    stellar_evolution = MESA()
    se_star = stellar_evolution.particles.add_particle(star)
    
    print "Evolving a", star.mass, "star with", stellar_evolution.__class__.__name__, 
    print "until its radius exceeds", stop_radius
    while (se_star.radius < stop_radius):
        se_star.evolve_one_step()
        temperatures_original.append(se_star.temperature)
        luminosities_original.append(se_star.luminosity)
    
    number_of_zones = se_star.get_number_of_zones()
    composition     = se_star.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
    index = (composition[0] > 1.0e-9).nonzero()[0][0] # first index where H fraction > 1.0e-9
    
    print "Creating helium star, from the inner", index, "(out of", str(number_of_zones)+") shells."
    helium_star = stellar_evolution.new_particle_from_model(dict(
        mass = (se_star.get_cumulative_mass_profile(number_of_zones = number_of_zones) * se_star.mass)[:index],
        radius = se_star.get_radius_profile(number_of_zones = number_of_zones)[:index],
        rho    = se_star.get_density_profile(number_of_zones = number_of_zones)[:index],
        temperature = se_star.get_temperature_profile(number_of_zones = number_of_zones)[:index],
        luminosity  = se_star.get_luminosity_profile(number_of_zones = number_of_zones)[:index],
        X_H  = composition[0][:index],
        X_He = composition[1][:index] + composition[2][:index],
        X_C  = composition[3][:index],
        X_N  = composition[4][:index],
        X_O  = composition[5][:index],
        X_Ne = composition[6][:index],
        X_Mg = composition[7][:index],
        X_Si = composition[7][:index]*0.0,
        X_Fe = composition[7][:index]*0.0), 0.0 | units.Myr)
    
    print "\nStar properties before helium star evolution:\n", stellar_evolution.particles
    for i in range(1000):
        helium_star.evolve_one_step()
        temperatures_helium.append(helium_star.temperature)
        luminosities_helium.append(helium_star.luminosity)
    print "\nStar properties after helium star evolution:\n", stellar_evolution.particles
    
    stellar_evolution.stop()
    return temperatures_original, luminosities_original, temperatures_helium, luminosities_helium

def plot_tracks(temperatures_original, luminosities_original, temperatures_helium, luminosities_helium):
    pyplot.figure(figsize = (8, 6))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    loglog(temperatures_original, luminosities_original, label='progenitor')
    loglog(temperatures_helium, luminosities_helium, label='helium star')
    xlabel('Effective Temperature')
    ylabel('Luminosity')
    pyplot.xlim(pyplot.xlim()[::-1])
    pyplot.ylim(1.0,1.0e5)
    pyplot.legend(loc=3)
    pyplot.show()


if __name__ in ('__main__', '__plot__'):        
    temperatures_original, luminosities_original, temperatures_helium, luminosities_helium = main()
    plot_tracks(temperatures_original, luminosities_original, temperatures_helium, luminosities_helium)
