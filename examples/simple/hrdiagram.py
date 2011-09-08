"""
Show how to generate Hertzsprung-Russell diagram
"""

import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel

from amuse.units import units
from amuse.community.sse.interface import SSE
from amuse import datamodel
from amuse.rfi.core import is_mpd_running
end_time = 2 | units.Gyr
stellar_mass =  2.0 | units.MSun

def simulate_evolution_tracks():
    stellar_evolution = SSE()
    stellar_evolution.commit_parameters()

    star = datamodel.Particle()
    star.mass = stellar_mass

    star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.commit_particles()
    
    luminosity_at_time = [] | units.LSun
    temperature_at_time = [] | units.K
    
    is_evolving = True
    while is_evolving and star.age < end_time:
        luminosity_at_time.append(star.luminosity)
        temperature_at_time.append(star.temperature)
        previous_age = star.age
        stellar_evolution.evolve_model()
        is_evolving = (star.age != previous_age)
        
    stellar_evolution.particles.remove_particle(star)

    stellar_evolution.stop()
    
    return temperature_at_time, luminosity_at_time
    
def plot_track(temperature_at_time, luminosity_at_time):
    pyplot.figure(figsize = (8, 6))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    loglog(temperature_at_time, luminosity_at_time)
    xlabel('Effective Temperature')
    ylabel('Luminosity')
    pyplot.xlim(pyplot.xlim()[::-1])
    pyplot.ylim(.1,1.e4)
    pyplot.show()
    

if __name__ in ('__main__', '__plot__'):        
    temperatures, luminosities = simulate_evolution_tracks()
    plot_track(temperatures, luminosities)
