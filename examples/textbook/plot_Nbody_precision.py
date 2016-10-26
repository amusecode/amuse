from amuse.lab import * 
import numpy

from amuse.ext.solarsystem import new_solar_system

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def energy_error_of_integrated_Nbody_system(code, particles, end_time, precision):

    #gravity = Huayno()
    #gravity = ph4()
    gravity = code()
    gravity.parameters.timestep_parameter = precision
    gravity.particles.add_particles(particles)
    channel_from_to_framework = gravity.particles.new_channel_to(particles)

    E0 = gravity.particles.potential_energy(G=nbody_system.G) 
    E0 += gravity.particles.kinetic_energy()
    gravity.evolve_model(end_time)
    channel_from_to_framework.copy()
    Et = gravity.particles.potential_energy(G=nbody_system.G) + gravity.particles.kinetic_energy()
    gravity.stop()
    de = (Et-E0)/E0
    return de

def get_dE(code, precision, t_end):
    dE = []
    for pri in precision:
        dEi = energy_error_of_integrated_Nbody_system(code, particles, t_end, pri)
        dE.append(abs(dEi))
        print "Integrated with eps=", pri, "dE/E=", dEi
    return dE
    
if __name__ in ('__main__','__plot__'):

    numpy.random.seed(31415)

    particles = new_plummer_model(100)
#    precision = numpy.exp(numpy.arange(2, -10, -1))
    precision = numpy.exp(numpy.arange(2, -8, -1))
    print precision

    from matplotlib import pyplot, rc
    x_label = "time step"
    y_label = "$|E(t)-E(0)|/E(0)$"
    figure = single_frame(x_label, y_label, logx=True, logy=True, xsize=12, ysize=10)
#    ax = pyplot.gca()
#    ax.set_xlim(0.0001, 10)
#    ax.set_ylim(1.e-16, 10)

    """
    from matplotlib import pyplot, rc
    figure = pyplot.figure(figsize=(10,10))
    font = {'size' : 20}
    rc('font', **font)
    plot = figure.add_subplot(1,1,1)
    pyplot.loglog()
    """

    cols = get_distinct(2)

    t_end = 1.0| nbody_system.time
    code = ph4
    dE = get_dE(code, precision, t_end)
    pyplot.scatter(precision, dE, c=cols[0], lw=0, s=200, marker='o')

    t_end = 1.0| nbody_system.time
    code = BHTree
    dE = get_dE(code, precision, t_end)
    pyplot.scatter(precision, dE, c=cols[1], lw=0, s=200, marker='^')

    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    pyplot.plot([4, 0.01], [0.1, (0.01)**5], c=cols[0], lw=4)
    pyplot.plot([4, 0.01], [0.02, (0.01)**2], c=cols[1], lw=4)
#    pyplot.show()

    pyplot.savefig("precision_N100t1")

    
