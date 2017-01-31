from amuse.lab import Particles, units 

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (0.0,0.0,0.0) | units.m
    sun.velocity = (0.0,0.0,0.0) | (units.m/units.s)

    venus = particles[1]
    venus.mass = 0.0025642 | units.MJupiter
    venus.radius = 3026.0 | units.km
    venus.position = (0.6335, 0.3499, -0.03179) | units.AU
    venus.velocity = (-17.0420, 30.5055, 1.4004) | units.kms

    earth = particles[2]
    earth.mass = 1.0 | units.MEarth
    earth.radius = 1.0 | units.REarth
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
    
    x_earth = [] | units.AU
    y_earth = [] | units.AU
    x_venus = [] | units.AU
    y_venus = [] | units.AU

    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (1 | units.day))
        x_earth.append(earth.x)
        y_earth.append(earth.y)
        x_venus.append(venus.x)
        y_venus.append(venus.y)
    gravity.stop()
    return x_earth, y_earth, x_venus, y_venus
    
def plot_track(xe,ye,xv,yv, output_filename):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(10, 10))
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)

    x_label = 'x [au]'
    y_label = 'y [au]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    plot.scatter([0.0], [0.0], color='y', lw=8)
    plot.plot(xe.value_in(units.AU), ye.value_in(units.AU), color = 'b')
    plot.plot(xv.value_in(units.AU), yv.value_in(units.AU), color = 'r')
    plot.set_xlim(-1.3, 1.3)
    plot.set_ylim(-1.3, 1.3)
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-o", 
                      dest="output_filename", default ="SunVenusEarth",
                      help="output filename [%default]")
    return result
    

if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = sun_venus_and_earth()
    xe,ye, xv,yv = integrate_solar_system(particles, 2 | units.yr)
    plot_track(xe,ye,xv,yv,o.output_filename)
    
