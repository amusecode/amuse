###BOOKLISTSTART1###
from amuse.lab import Particles, units 

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (855251, -804836, -3186) |units.km
    sun.velocity = (7.893, 11.894, 0.20642) | (units.m/units.s)

    venus = particles[1]
    venus.mass = 0.0025642 | units.MJupiter
    venus.radius = 3026.0 | units.km
    venus.position = (-0.3767, 0.60159, 0.03930) | units.AU
    venus.velocity = (-29.7725, -18.849, 0.795) | units.kms

    earth = particles[2]
    earth.mass = 1.0 | units.MEarth
    earth.radius = 1.0 | units.REarth
    earth.position = (-0.98561, 0.0762, -7.847e-5) | units.AU  
    earth.velocity = (-2.927, -29.803, -0.0005327) | units.kms

    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.lab import Huayno, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(),
                                             particles[1].position.length())

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
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(xe,ye,xv,yv, output_filename):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(10, 10))
    pyplot.rcParams.update({'font.size': 30})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() 
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

    save_file = 'sun_venus_earth.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
###BOOKLISTSTOP3###

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
    plot_track(xe, ye, xv, yv, o.output_filename)
    
