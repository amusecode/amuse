"""
   Simple routine for running a gravity code
"""
from amuse.lab import *
from matplotlib import pyplot
from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def virial_ratio_evolution(code, bodies, Q_init, t_end):
    dt = 0.06125 | t_end.unit
    bodies.scale_to_standard(virial_ratio=Q_init)
    bodies.radius = 0 | nbody_system.length
    gravity = code()
    gravity.particles.add_particles(bodies)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)

    Etot_prev = Etot_init = gravity.kinetic_energy + gravity.potential_energy
    time = [0.0] | t_end.unit
    Q = [Q_init]
    while time[-1] < t_end:
        time.append(time[-1]+dt)
        gravity.evolve_model(time[-1])
        channel_from_gravity_to_framework.copy()
        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        Q.append(-1*Ekin/Epot)
        print "T=", time[-1], "Q= ", Q[-1],
        print "M=", bodies.mass.sum(), "E= ", Etot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot
    gravity.stop()
    return time, Q

def main(N, t_end):
    t_end = t_end | nbody_system.time
    Q_init = 0.2
    particles = new_plummer_model(N)
    codes = [ph4, Huayno, BHTree]
    cols = get_distinct(3)
    ci = 0
    x_label = "time [N-body units]"
    y_label = "virial ratio $Q$"
    figure = single_frame(x_label, y_label, xsize=14, ysize=10)
    ax1 = pyplot.gca()
    ax1.set_xlim(0, t_end.value_in(t_end.unit))
    ax1.set_ylim(0, 0.65)
    pyplot.plot([0, t_end.value_in(t_end.unit)], [0.5, 0.5],
                lw=1, ls='--', c='k')
    for code in codes:
        time, Q = virial_ratio_evolution(code, particles, Q_init, t_end)
        pyplot.plot(time.value_in(t_end.unit), Q, c=cols[ci])
        ci+=1

    save_file = 'gravity_to_virial.png'
    pyplot.savefig(save_file)
    print "\nOutput saved in", save_file, '\n'
    pyplot.show()
    
def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [10]")
    result.add_option("-t", dest="t_end", type="float", default = 2,
                      help="end time of the simulation [1] N-body units")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

