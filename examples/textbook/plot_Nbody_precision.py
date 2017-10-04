from amuse.lab import * 
import numpy
from distinct_colours import get_distinct
from matplotlib import pyplot

def energy_error_of_integrated_Nbody_system(code, particles,
                                            end_time, precision):

    gravity = code(number_of_workers=4)
    gravity.parameters.timestep_parameter = precision
    #gravity.parameters.timestep = precision | nbody_system.time

    gravity.particles.add_particles(particles)
    channel_from_to_framework = gravity.particles.new_channel_to(particles)

    E0 = gravity.particles.potential_energy(G=nbody_system.G) 
    E0 += gravity.particles.kinetic_energy()
    gravity.evolve_model(end_time)
    channel_from_to_framework.copy()
    Et = gravity.particles.potential_energy(G=nbody_system.G) \
                            + gravity.particles.kinetic_energy()
    gravity.stop()

    de = (Et-E0)/E0
    return de

def get_dE(code, precision, t_end):
    dE = []
    for pri in precision:
        dEi = energy_error_of_integrated_Nbody_system(code, particles,
                                                      t_end, pri)
        dE.append(abs(dEi))
        print "integrated with precision=", pri, "dE/E=", dEi
    return dE
    
if __name__ in ('__main__','__plot__'):

    numpy.random.seed(31415)

    particles = new_plummer_model(1000)
    precision = 10.**numpy.linspace(0., -3., 10)

    t_end = 1.0| nbody_system.time
    cols = get_distinct(2)
    
    print 'ph4'
    code = ph4
    dE = get_dE(code, precision, t_end)
    pyplot.scatter(precision, dE, c=cols[0], lw=0, s=50, marker='o')

    print 'BHTree'
    code = BHTree
    dE = get_dE(code, precision, t_end)
    pyplot.scatter(precision, dE, c=cols[1], lw=0, s=50, marker='^')

    t0 = 0.8
    t1 = 0.02
    ep = 4.e-5
    eb = 0.07
    pyplot.plot([t0, t1], [ep, ep*(t1/t0)**4], c=cols[0], lw=2)
    pyplot.plot([t0, t1], [eb, eb*(t1/t0)**2], c=cols[1], lw=2)

    pyplot.xlabel('time step parameter')
    pyplot.xlim(1.e-4, 3.)
    pyplot.xscale('log')
    pyplot.ylabel('$|E(t)-E(0)|/|E(0)|$')
    pyplot.ylim(1.e-15, 10.)
    pyplot.yscale('log')

    save_file = 'precision_N100t1.png'
    pyplot.savefig(save_file)
    print "\nOutput saved in", save_file
    pyplot.show()

    
