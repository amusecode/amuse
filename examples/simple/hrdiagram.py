import sys
import numpy 
import os
import warnings
from matplotlib import pyplot

from amuse.support.units import units
from amuse.support.data import core
from amuse.community.sse.interface import SSE
from amuse.support.codes.core import is_mpd_running

def stellar_remnant_state(star):
    return 10 <= star.stellar_type.value_in(units.stellar_type) and \
        star.stellar_type.value_in(units.stellar_type) < 16

def simulate_evolution_tracks():
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_current_parameters()
    star = core.Particle()
    star.mass = 2.0 | units.MSun

    star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.initialize_stars()
    
    luminosity_at_time     = [] | units.LSun
    temperature_at_time     = [] | units.K
    stellar_type_at_time = [] | units.stellar_type
        
    stopped_evolving = False
    while not stellar_remnant_state(star) and not stopped_evolving:
        luminosity_at_time.append(star.luminosity)
        temperature_at_time.append(star.temperature)
        stellar_type_at_time.append(star.stellar_type)
        previous_age = star.age
        try:
            stellar_evolution.evolve_model()
            stopped_evolving = (star.age == previous_age) # Check whether the age has stopped increasing
        except Exception as ex:
            stopped_evolving = True
        
    stellar_evolution.particles.remove_particle(star)

    stellar_evolution.stop()
    plot_track(temperature_at_time, luminosity_at_time)
    
def plot_track(temperature_at_time, luminosity_at_time):
    pyplot.figure(figsize = (7, 8))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    pyplot.xlabel('Effective Temperature (K)')
    pyplot.ylabel('Luminosity (solar luminosity)')
    
    x_values = temperature_at_time.value_in(units.K)
    y_values = luminosity_at_time.value_in(units.LSun)
    pyplot.loglog(x_values, y_values)
    pyplot.show()
    

if __name__ in ['__main__', '__plot__']:
    if not is_mpd_running():
        print "There is no mpd server running. Please do 'mpd &' first."
        sys.exit()
        
    simulate_evolution_tracks()
