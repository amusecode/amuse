from amuse.lab import * 

from amuse.ext.solarsystem import solar_system_in_time
from prepare_figure import single_frame, figure_frame, \
    			   set_tickmarks, get_distinct

def integrate_solar_system(particles, end_time):

    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(),
                                             particles[1].position.length())

    gravity = Huayno(convert_nbody)
    gravity.particles.add_particles(particles)

    SunEarth = Particles()
    SunEarth.add_particle(particles[0])
    SunEarth.add_particle(particles[3])
    channel_from_to_SunEarth = gravity.particles.new_channel_to(SunEarth)

    from matplotlib import pyplot, rc
    x_label = "$a-a_0$ [AU]"
    y_label = "eccentricty"
    figure = single_frame(x_label, y_label, xsize=14, ysize=10)

    from amuse.plot import scatter
    from matplotlib import pyplot, rc
    
    kep = Kepler(convert_nbody)
    kep.initialize_from_particles(SunEarth)
    a, e = kep.get_elements()
    a0 = a
    color = get_distinct(1)
    scatter((a-a0).in_(units.AU), e, c=color[0], lw=0)
    pyplot.text((a-a0).value_in(units.AU)-0.0003, e,
            "{0:.0f}".format(0.001*gravity.model_time.value_in(units.yr))+"kyr")

    dt = 100 | units.yr
    t_diag = 10000 |units.yr
    while gravity.model_time < end_time:

        kep.initialize_from_particles(SunEarth)
        a, e = kep.get_elements()
        scatter((a-a0).in_(units.AU), e, c=color[0], lw=0, s=80)

        time = gravity.model_time
        if time > t_diag:
            t_diag += 10000 |units.yr
            fmt = "{0:.0f}".format(0.001*time.value_in(units.yr))+" kyr"
            if time < 100|units.yr:
                pyplot.text(-0.0006, e, fmt)
            elif time < 15000|units.yr:
                pyplot.text(0.0003, e, fmt)
            elif time < 25000|units.yr:
                pyplot.text(0.0017, e, fmt)
            elif time < 35000|units.yr:
                pyplot.text(0.0031, e, fmt)
            elif time < 45000|units.yr:
                pyplot.text(0.0028, e, fmt)
            elif time < 55000|units.yr:
                pyplot.text(0.0021, e, fmt)
            elif time < 65000|units.yr:
                pyplot.text(0.0017, e+0.002, fmt)
            elif time < 75000|units.yr:
                pyplot.text(0.0014, e-0.002, fmt)
            else:
                pass

        gravity.evolve_model(gravity.model_time + dt)
        channel_from_to_SunEarth.copy()
    gravity.stop()

    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.xlim(-0.001, 0.004)

    file = 'EarthOrbitVariation.png'
    pyplot.savefig(file)
    print '\nSaved figure in file', file,'\n'
    pyplot.show()

    return
    
if __name__ in ('__main__','__plot__'):
    particles = solar_system_in_time(time_JD=2474649.5|units.day)
    integrate_solar_system(particles, 84000 | units.yr)
    
