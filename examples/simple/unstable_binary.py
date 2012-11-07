"""
Evolves two stars dynamically (hermit, nbody code) each star will 
lose mass during the evolution (evtwin, stellar evolution code)

We start with two stars, one 10.0 and one 1.0 solar mass star. These
stars start orbiting with a stable kepler orbit. 
After 2 orbital periods the stars will begin to lose mass and the binary
will become unstable.
"""

from amuse.plot import scatter, xlabel, ylabel, plot
from matplotlib import pyplot
from math import pi
from amuse.units.optparse import OptionParser

from amuse.units import units
from amuse.units import constants
from amuse.units.nbody_system import nbody_to_si
from amuse.community.evtwin.interface import EVtwin
from amuse.community.sse.interface import SSE
from amuse.community.hermite0.interface import Hermite

from amuse.datamodel import Particles


def set_up_initial_conditions(orbital_period, kinetic_to_potential_ratio):
    print "Setting up initial conditions"
    stars =  Particles(2)
    stars.mass = [10.0, 1.0] | units.MSun
    stars.radius = 0 | units.RSun
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    
    print "Binary with masses: "+str(stars.mass)+", and orbital period: ", orbital_period
    semimajor_axis = ((constants.G * stars.total_mass() * (orbital_period / (2 * pi))**2.0)**(1.0/3.0))
    separation = 2 * semimajor_axis * (1 - kinetic_to_potential_ratio)
    print "Initial separation:", separation.as_quantity_in(units.AU)
    relative_velocity = ( (kinetic_to_potential_ratio / (1.0 - kinetic_to_potential_ratio)) * 
        constants.G * stars.total_mass() / semimajor_axis).sqrt()
    print "Initial relative velocity:", relative_velocity.as_quantity_in(units.km / units.s)
    
    stars[0].x = separation
    stars[0].vy = relative_velocity
    stars.move_to_center()
    return stars

def set_up_stellar_evolution_code(stars):
    stellar_evolution = EVtwin()
    stellar_evolution.initialize_code()
    # if you run with mesa, you can play with the wind efficiency
    # stellar_evolution.parameters.RGB_wind_scheme = 1
    # stellar_evolution.parameters.reimers_wind_efficiency = 1.0e6 # ridiculous, but instructive
    stellar_evolution.particles.add_particles(stars)
    return stellar_evolution
    
def set_up_gravitational_dynamics_code(stars):
    convert_nbody = nbody_to_si(11.0 | units.MSun, 10.0 | units.AU)
    gravitational_dynamics = Hermite(convert_nbody)
    gravitational_dynamics.parameters.epsilon_squared = 0.0 | units.AU ** 2
    view_on_the_primary = gravitational_dynamics.particles.add_particle(stars[0])
    gravitational_dynamics.particles.add_particle(stars[1])
    return gravitational_dynamics, view_on_the_primary
    

def simulate_binary_evolution(binary, orbital_period, t_offset_stars, t_end):
    distance = [] | units.AU
    mass = [] | units.MSun
    time = [] | units.yr
    
    stellar_evolution = set_up_stellar_evolution_code(binary)
    gravitational_dynamics, primary = set_up_gravitational_dynamics_code(binary)
    from_se_to_gd = stellar_evolution.particles.new_channel_to(gravitational_dynamics.particles)
    
    current_time = 0.0 * t_end
    
        
    print "Evolving with stellar wind"
    while current_time < t_end:
        current_time += orbital_period / 10
        gravitational_dynamics.evolve_model(current_time)
        stellar_evolution.evolve_model(current_time + t_offset_stars)
        from_se_to_gd.copy_attributes(['mass'])
        separation = (gravitational_dynamics.particles[0].position - gravitational_dynamics.particles[1].position).length()
        distance.append(separation)
        mass.append(primary.mass)
        time.append(current_time)
        print "System evolved to time: ", current_time, 
        print ", primary mass:", primary.mass.as_quantity_in(units.MSun),
        print ", separation:", separation.as_quantity_in(units.AU)
    
    print "Evolution done"
    return distance, mass, time

def orbit_plot(distance, mass, time):
    figure = pyplot.figure(figsize = (6, 10), dpi = 100)
    subplot = figure.add_subplot(2, 1, 1)
    plot(time, distance)
    xlabel('t')
    ylabel('separation')
    pyplot.margins(0.05)
    subplot = figure.add_subplot(2, 1, 2)
    plot(time, ((mass - mass[0]) / mass[0]) * 100.0)
    xlabel('t')
    ylabel('mass')
    pyplot.margins(0.05)
    pyplot.show()

def main(
        orbital_period = 1000.0 | units.yr, 
        kinetic_to_potential_ratio = 0.8, 
        periods = 10,
        age = 10 | units.Myr
    ):
    
    t_offset_stars = age
    t_end = periods * orbital_period
    
    binary = set_up_initial_conditions(orbital_period, kinetic_to_potential_ratio)
    distance, mass, time = simulate_binary_evolution(binary, orbital_period, t_offset_stars, t_end)
    orbit_plot(distance, mass, time)
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-o", "--orbitalperiod", 
        default = 1000 | units.yr,
        dest="orbital_period",
        help="initial orbital period of the binary (in years)",
        type="float",
        unit=units.yr
    )
    
    result.add_option(
        "-k", "--kpratio", 
        default = 0.8,
        dest="kinetic_to_potential_ratio",
        help="kinetec to potential energy ratio, values less than 1.0 correspond to bound systems",
        type="float"
    )
    result.add_option(
        "--periods", 
        default = 10,
        dest="periods",
        help="number of orbital periods to evolve the binary",
        type="int"
    )
    result.add_option(
        "--age", 
        default = 10 | units.Myr,
        dest="age",
        help="initial age of the stars to start the simulation with",
        type="float",
        unit=units.Myr
    )
    
    
    return result
    
    
if __name__ == "__plot__":
    main(1000 | units.yr, 0.8, 10, 10 | units.Myr)
    
if __name__ == "__main__":
    options, args = new_option_parser().parse_args()
    main(**options.__dict__)
