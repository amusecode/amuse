import numpy
import time as _time
from matplotlib import pyplot

from amuse.units import units
from amuse.units import nbody_system
from amuse.community.hermite0.interface import Hermite
from amuse.community.ph4.interface import ph4
from amuse.community.huayno.interface import Huayno
try:
    from amuse.community.PyNbody.interface import PyNbody
except ImportError:
    pass
import logging

from amuse.ic.plummer import new_plummer_model
#logging.basicConfig(level=logging.DEBUG)

#smoothing_length = 0.0 | nbody_system.length ** 2

numpy.random.seed(123456)

def print_log(time, gravity, particles, total_energy_at_t0):
    kinetic_energy = gravity.kinetic_energy
    potential_energy = gravity.potential_energy
    total_energy_at_this_time = kinetic_energy + potential_energy
    de_e0 = (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
    print "time                    : " , time
    print "energy error            : " , (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0

    
def simulate_small_cluster(
        number_of_stars = 1000, 
        end_time = 40 | nbody_system.time,
        number_of_workers = 1
    ):
    particles = new_plummer_model(number_of_stars)
    particles.scale_to_standard()

    p0 = particles.copy()
    p1 = particles.copy()
    p2 = particles.copy()
    p3 = particles.copy()

    for meth, particles, tstep in ((ph4, p0, 0.2), (Huayno, p1, 0.03), (PyNbody, p2, 0.03), (Hermite, p3, None)):
        gravity = meth()
        print "\nworker code:", gravity.__class__.__name__

        gravity.parameters.epsilon_squared = ((4.0/len(particles))**1) | nbody_system.length ** 2
        if tstep:
            gravity.parameters.timestep_parameter = tstep

        gravity.particles.add_particles(particles)
        
        from_model_to_gravity = particles.new_channel_to(gravity.particles)
        from_gravity_to_model = gravity.particles.new_channel_to(particles)
        
        time = 0.0 | end_time.unit
        total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy
        
        positions_at_different_times = []
        positions_at_different_times.append(particles.position)
        times = []
        times.append(time)

        fname = "log_"+gravity.__class__.__name__
        total_energy_at_this_time = gravity.kinetic_energy + gravity.potential_energy
        de_e0 = (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
        with open(fname, 'w') as fobj:
            fobj.write("{0} {1}\n".format(time.value_in(nbody_system.time), de_e0))

        print "evolving the model until t = " + str(end_time)
        start = _time.time()
        while time < end_time:
            time +=  end_time / 9.0
            
            gravity.evolve_model(time)
            from_gravity_to_model.copy()

#            print gravity.get_time()
            
            positions_at_different_times.append(particles.position)
            times.append(time)
            print_log(time, gravity, particles, total_energy_at_t0)
            
            total_energy_at_this_time = gravity.kinetic_energy + gravity.potential_energy
            de_e0 = (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
            with open(fname, 'a') as fobj:
                fobj.write("{0} {1}\n".format(time.value_in(nbody_system.time), de_e0))

        gravity.stop()

        print "{} CPU time: {} seconds".format(gravity.__class__.__name__, _time.time() - start)

#        plot_positions(times, positions_at_different_times)
    
    return times, positions_at_different_times


def adjust_spines(ax,spines, ticks):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_ticks(ticks)
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_ticks(ticks)

    else:
        ax.xaxis.set_ticks([])


def plot_positions(times, positions_at_different_times):
    figure = pyplot.figure()
    plot_matrix_size = numpy.ceil(numpy.sqrt(len(positions_at_different_times))).astype(int)
    number_of_rows = len(positions_at_different_times) / plot_matrix_size
    figure.subplots_adjust(wspace = 0.15, hspace = 0.15)
    
    for index, (time, positions) in enumerate(zip(times, positions_at_different_times)):
        subplot = figure.add_subplot(plot_matrix_size, plot_matrix_size, index + 1)
        
        subplot.scatter(
            positions[...,0].value_in(nbody_system.length),
            positions[...,1].value_in(nbody_system.length),
            s = 1,
            edgecolors = 'red',
            facecolors = 'red'
        )
        
        subplot.set_xlim(-4.0,4.0)
        subplot.set_ylim(-4.0,4.0)
        
        title = 'time = {0:.2f}'.format(time.value_in(nbody_system.time))
        
        subplot.set_title(title)#, fontsize=12)
        spines = []
        if index % plot_matrix_size == 0:
            spines.append('left')

        if index >= ((number_of_rows - 1)*plot_matrix_size):
            spines.append('bottom')
            

        adjust_spines(subplot, spines,numpy.arange(-4.0,4.1, 1.0))
        
        if index % plot_matrix_size == 0:
            subplot.set_ylabel('y')
            
        if index >= ((number_of_rows - 1)*plot_matrix_size):
            subplot.set_xlabel('x')
            
    pyplot.show()
    
if __name__ in ('__main__', '__plot__'):
    times, positions_at_different_times = simulate_small_cluster(
        512,
        9.0 | nbody_system.time
    )
#    plot_positions(times, positions_at_different_times)
