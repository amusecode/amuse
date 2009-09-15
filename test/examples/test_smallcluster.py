import random
import sys

import numpy 
from matplotlib import pyplot

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.muse_dynamics_mpi import Hermite
from amuse.legacy.bhtree.muse_dynamics_mpi import BHTree
from amuse.legacy.sse.muse_stellar_mpi import SSE

from amuse.experiments.plummer import MakePlummerModel


def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    initial_mass = 6
    convert_nbody = nbody_system.nbody_to_si(number_of_stars * initial_mass | units.MSun, 1 | units.lightyear)
    
    particles = MakePlummerModel(number_of_stars, convert_nbody).result;
            
    gravity = BHTree(convert_nbody)
    gravity.setup_module()
    gravity.dt_dia = 10000
    
    sse = SSE()
    sse.initialize_module_with_default_parameters() 
    
    gravity.set_eps2(0.3 | units.lightyear ** 2)
    
    print "initializing the particles"
    sse.initialize_particles(particles)
    gravity.add_particles(particles)
    
    print "setting age of stars"
    for x in particles:
        start_time = random.randint(40, 115) | units.Myr
        sse.evolve_particle(x, start_time )
        x.start_time = x.current_time.value()
        
    time = 1 | units.Myr
    
    
    print "evolving the model until t = " + str(end_time)
    while time < end_time:
        gravity.evolve_model(time)
        
        gravity.update_particles(particles)
        
        for x in particles:
            sse.evolve_particle(x, x.start_time.value() + time)
        
        print "evolved model to t = " + str(time)
        
        gravity.update_attributes(particles.mass)
        time += 2 |units.Myr
    
   
    print "plotting the data"
    figure = pyplot.figure(figsize = (40,40))
    plots = map(lambda x : figure.add_subplot(5,5,x), range(1,int(int(end_time.number) / 2 + 2)))
    
    for x in particles:
        for index, (time,position) in enumerate(x.position.values):
            x_point = position.in_(units.lightyear).number[0]
            y_point = position.in_(units.lightyear).number[1]
            t, mass = x.mass.get_value_at_time(x.start_time.value() + time)
            
            color = 'b'
            if mass > 4.999 | units.MSun:
                color = 'b'
            elif mass > 4 | units.MSun:
                color = 'y'
            elif mass > 2 | units.MSun:
                color = 'm'
            elif mass > 1 | units.MSun:
                color = 'r'
            else:
                color = 'c'
                
            plots[index].plot([x_point],[y_point], color + 'o')
        
    for plot in plots:
        plot.set_xlim(-10.0,10.0)
        plot.set_ylim(-10.0,10.0)

    figure.savefig(name_of_the_figure)
    
    del gravity
    del sse
    
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    simulate_small_cluster(5, 4 | units.Myr)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | units.Myr, sys.argv[3])
