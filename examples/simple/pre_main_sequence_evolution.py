"""
Theoretical pre-main-sequance evolutionary tracks for stars of various masses.
After Iben, ApJ 141, 993, 1965
"""

from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel

def simulate_evolution_tracks():
    stellar_evolution = MESA()
    stellar_evolution.parameters.AGB_wind_scheme = 0
    stellar_evolution.parameters.RGB_wind_scheme = 0

    masses = [0.5, 1.0, 1.25, 1.5, 2.25, 3.0, 5.0, 9.0, 15.0] | units.MSun
    stars = datamodel.Particles(len(masses))
    stars.mass = masses
    stars = stellar_evolution.pre_ms_stars.add_particles(stars)
    
    luminosity_at_time = []
    temperature_at_time = []
    time = []
    
    for star in stars:
        one_luminosity_at_time = [] | units.LSun
        one_temperature_at_time = [] | units.K
        one_time = [] | units.yr
        while star.stellar_type == 17 | units.stellar_type:
            one_luminosity_at_time.append(star.luminosity)
            one_temperature_at_time.append(star.temperature)
            one_time.append(star.age)
            star.evolve_one_step()
            print star.stellar_type, star.age, star.mass, star.luminosity, star.radius
        luminosity_at_time.append(one_luminosity_at_time)
        temperature_at_time.append(one_temperature_at_time)
        time.append(one_time)
        
    stellar_evolution.stop()
    
    return temperature_at_time, luminosity_at_time, time
    
def plot_track(temperature_at_time, luminosity_at_time):
    pyplot.figure(figsize = (6, 8))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    for temperature, luminosity in zip(temperature_at_time, luminosity_at_time):
        loglog(temperature[5:], luminosity[5:], marker="s")
    xlabel('Effective Temperature')
    ylabel('Luminosity')
    pyplot.xlim(10**4.6, 10**3.5)
    pyplot.ylim(1.0e-2,1.e5)
    pyplot.show()
    

if __name__ in ('__main__', '__plot__'):        
    temperatures, luminosities, time = simulate_evolution_tracks()
    plot_track(temperatures, luminosities)
    
