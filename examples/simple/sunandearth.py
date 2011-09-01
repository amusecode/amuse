from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.units import values

from amuse.community.hermite0.interface import Hermite

from matplotlib import pyplot
def new_system_of_sun_and_earth():
    stars = core.Particles(2)
    
    sun = stars[0]
    sun.mass = 1.0 | units.MSun
    sun.position = (0.0,0.0,0.0) | units.m
    sun.velocity = (0.0,0.0,0.0) | (units.m / units.s)
    sun.radius = 1.0 | units.RSun

    earth = stars[1]
    earth.mass = 5.9736e24 | units.kg
    earth.radius = 6371.0 | units.km
    earth.position = (149.5e6,0.0,0.0) | units.km
    earth.velocity = (0.0,29800,0.0) | (units.m / units.s) 
    
    return stars

def simulate_system_until(particles, end_time):
    convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

    instance = Hermite(convert_nbody)
    instance.parameters.epsilon_squared = 0.0 | units.AU**2
    instance.particles.add_particles(particles)
    
    
    t0 = 0 | units.day
    dt = 10 | units.day
    t = t0
    earth = instance.particles[1]
    
    x_values = values.AdaptingVectorQuantity()
    y_values = values.AdaptingVectorQuantity()
    
    while t < end_time:
        instance.evolve_model(t)
        
        x_values.append(earth.x)
        y_values.append(earth.y)
        
        t += dt
    
    instance.stop()
    
    return x_values, y_values
    
    
def plot_track(x,y):

    figure = pyplot.figure(figsize=(5,5))
    plot = figure.add_subplot(1,1,1)
    
    x_values_in_AU = x.value_in(units.AU)
    y_values_in_AU = y.value_in(units.AU)
    
    plot.plot(x_values_in_AU,y_values_in_AU, color = "b")
    
    plot.set_xlim(-1.5, 1.5)
    plot.set_ylim(-1.5, 1.5)
    
    plot.set_xlabel('x (AU)')
    plot.set_ylabel('y (AU)')
               
    pyplot.show()

if __name__ in ('__main__','__plot__'):
    particles = new_system_of_sun_and_earth()
    x,y = simulate_system_until(particles, 20 | units.yr)
    plot_track(x,y)
    
