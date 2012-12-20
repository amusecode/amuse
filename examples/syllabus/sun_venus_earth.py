"""
   Example AMUSE script for evolving Venus and Earth around the Sun.
   Syllabus List. 2.1 and 2.2
"""
from amuse.lab import Particles, units 

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.position = (0.0,0.0,0.0) | units.m
    sun.velocity = (0.0,0.0,0.0) | (units.m/units.s)
    sun.radius = 1.0 | units.RSun

    venus = particles[1]
    venus.mass = 0.0025642 | units.MJupiter
    venus.radius = 3026.0 | units.km
    venus.position = (0.6335, 0.3499, -0.03179) | units.AU
    venus.velocity = (-17.0420, 30.5055, 1.4004) | units.kms

    earth = particles[2]
    earth.mass = 5.9736e24 | units.kg
    earth.radius = 6371.0 | units.km
    earth.position = (0.2421, -0.9875, -0.00004) | units.AU
    earth.velocity = (28.4468, 6.98125, 0.0002) | units.kms

    particles.move_to_center()
    return particles

def integrate_solar_system(particles, end_time):
    from amuse.lab import Huayno, nbody_system
    from amuse.units import quantities
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[1].position.length())

    gravity = Huayno(convert_nbody)
    gravity.particles.add_particles(particles)
    venus = gravity.particles[1]
    earth = gravity.particles[2]
    
    x_earth = quantities.AdaptingVectorQuantity()
    y_earth = quantities.AdaptingVectorQuantity()
    x_venus = quantities.AdaptingVectorQuantity()
    y_venus = quantities.AdaptingVectorQuantity()
    
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (10 | units.day))
        x_earth.append(earth.x)
        y_earth.append(earth.y)
        x_venus.append(venus.x)
        y_venus.append(venus.y)
    gravity.stop()
    return x_earth, y_earth, x_venus, y_venus
    
def plot_track(xe,ye,xv,yv):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(5,5))
    plot = figure.add_subplot(1,1,1)

    plot.scatter([0.0], [0.0], color='y')
    plot.plot(xe.value_in(units.AU), ye.value_in(units.AU), color = "b")
    plot.plot(xv.value_in(units.AU), yv.value_in(units.AU), color = "r")
    plot.set_xlim(-1.5, 1.5)
    plot.set_ylim(-1.5, 1.5)
    plot.set_xlabel('x (AU)')
    plot.set_ylabel('y (AU)')
    pyplot.show()

if __name__ in ('__main__','__plot__'):
    particles = sun_venus_and_earth()
    xe,ye, xv,yv = integrate_solar_system(particles, 2 | units.yr)
    plot_track(xe,ye,xv,yv)
    
