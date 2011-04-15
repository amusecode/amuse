from amuse.plot import scatter, xlabel, ylabel, plot
from matplotlib import pyplot
from math import pi

from amuse.support.data.core import Particles
from amuse.support.units import units, constants
from amuse.support.units.nbody_system import nbody_to_si
from amuse.community.mesa.interface import MESA
from amuse.community.hermite0.interface import Hermite


orbital_period = 1000.0 | units.yr # initial orbital period, will increase when the binary loses mass
kinetic_to_potential_ratio = 0.8 # (K/|U|), values < 1.0 correspond to bound systems
t_start_evolution_with_wind = 3.0 * orbital_period
t_end = 20.0 * orbital_period


def set_up_initial_conditions():
    print "Setting up initial conditions"
    stars =  Particles(2)
    stars.mass = [10.0, 1.0] | units.MSun
    stars.radius = 0 | units.RSun
    stars.position = [[0.0]*3]*2 | units.AU
    stars.velocity = [[0.0]*3]*2 | units.km / units.s
    
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
    stellar_evolution = MESA()
    stellar_evolution.initialize_code()
    stellar_evolution.parameters.RGB_wind_scheme = 1
    stellar_evolution.parameters.reimers_wind_efficiency = 1.0e6 # ridiculous, but instructive
    stellar_evolution.particles.add_particles(stars)
    return stellar_evolution
    
def set_up_gravitational_dynamics_code(stars):
    convert_nbody = nbody_to_si(11.0 | units.MSun, 10.0 | units.AU)
    gravitational_dynamics = Hermite(convert_nbody)
    gravitational_dynamics.parameters.epsilon_squared = 0.0 | units.AU ** 2
    view_on_the_primary = gravitational_dynamics.particles.add_particle(stars[0])
    gravitational_dynamics.particles.add_particle(stars[1])
    return gravitational_dynamics, view_on_the_primary
    
def simulate_binary_evolution(binary):
    distance = [] | units.AU
    mass = [] | units.MSun
    time = [] | units.yr
    
    stellar_evolution = set_up_stellar_evolution_code(binary)
    gravitational_dynamics, primary = set_up_gravitational_dynamics_code(binary)
    from_se_to_gd = stellar_evolution.particles.new_channel_to(gravitational_dynamics.particles)
    
    current_time = 0.0 | units.yr
    print "Evolving without stellar wind"
    while current_time < t_start_evolution_with_wind:
        current_time += orbital_period / 100
        gravitational_dynamics.evolve_model(current_time)
        distance.append((gravitational_dynamics.particles[0].position - gravitational_dynamics.particles[1].position).length())
        mass.append(primary.mass)
        time.append(current_time)
    
    print "Evolving with stellar wind"
    while current_time < t_end:
        current_time += orbital_period / 100
        gravitational_dynamics.evolve_model(current_time)
        stellar_evolution.evolve_model(current_time - t_start_evolution_with_wind)
        from_se_to_gd.copy_attributes(['mass'])
        distance.append((gravitational_dynamics.particles[0].position - gravitational_dynamics.particles[1].position).length())
        mass.append(primary.mass)
        time.append(current_time)
    
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
    plot(time, mass)
    xlabel('t')
    ylabel('mass')
    pyplot.margins(0.05)
    pyplot.show()

if __name__ in ["__main__", "__plot__"]:
    binary = set_up_initial_conditions()
    distance, mass, time = simulate_binary_evolution(binary)
    orbit_plot(distance, mass, time)
