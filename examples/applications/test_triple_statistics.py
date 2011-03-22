"""
In this script we simulate triple evolution.

See papers:
http://arxiv.org/abs/1004.2506
and
http://arxiv.org/abs/1101.0399
"""


import numpy
import sys
import time
from amuse.lab import *
from amuse.support.io import text

from multiprocessing import Process,Queue

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


from mpi4py import MPI

class SimulateTripleSystemUntilDecay(object):
    gravitational_constant = nbody_system.G
    distance_relative_to_inner_binary_axis = 50
    
    t0 = 0 | nbody_system.time
    tmax = 7000 | nbody_system.time
    dt = 10 | nbody_system.time
        
    def __init__(self, gravity_code, particles, rank = 0, index = 0):
        self.particles = particles
        self.gravity_code = gravity_code
        self.gravity_code.particles.add_particles(self.particles)
        self.from_code_to_model_channel = self.gravity_code.particles.new_channel_to(self.particles)
        self.rank = rank
        self.index = index
        
    def is_hyperbolic(self, inner_binary, outer_particle):
        mass_inner = inner_binary.mass.sum()
        velocity_cm = inner_binary.center_of_mass_velocity()
        position_cm = inner_binary.center_of_mass()
        distance_to_cm = (outer_particle.position - position_cm).length()
        
        energy = (
            ((mass_inner * (velocity_cm.length_squared())) / 2.0) +
            ((outer_particle.mass * (outer_particle.velocity.length_squared())) / 2.0) + 
            - self.gravitational_constant * mass_inner * outer_particle.mass /  distance_to_cm
        ) # formula (3) in Orlov(2010)
        
        return energy > 0 | nbody_system.energy

    def is_outside_escape_distance(self, inner_binary, outer_particle):
        position_cm = inner_binary.center_of_mass()
        distance_to_cm = (outer_particle.position - position_cm).length()
        distance_between_inner_particles = (inner_binary[0].position - inner_binary[1].position).length()
        return (distance_to_cm / distance_between_inner_particles) > self.distance_relative_to_inner_binary_axis

    def has_triple_an_escaper(self, particles):
        breakup_of_indices = (
            ((0,1),2),
            ((0,2),1),
            ((1,2),0)
        )
        for binary_indices, outer_particle_index in breakup_of_indices:
            binary = particles.get_all_particles_at(*binary_indices)
            outer = particles[outer_particle_index]
            is_possible_escape = self.is_outside_escape_distance(
                binary,
                outer
            )
            if is_possible_escape:
                is_on_hyperbolic_trajectory = self.is_hyperbolic(
                    binary,
                    outer
                )
                if is_on_hyperbolic_trajectory:
                    return True
                    
        return False
        
    def evolve_model_until_escape(self):
        self.time = self.t0
        has_escaper = False
        
        self.setup_point_positions()
        
        while self.time < self.tmax and not has_escaper:
            self.gravity_code.evolve_model(self.time)
            self.time += self.dt
            self.from_code_to_model_channel.copy()
            has_escaper = self.has_triple_an_escaper(self.particles)
            self.store_point_positions()
        
        self.make_plot()
        
        return self.time
    
    def setup_point_positions(self):
        self.x_positions = [AdaptingVectorQuantity() for x in range(len(self.particles))]
        self.y_positions = [AdaptingVectorQuantity() for x in range(len(self.particles))]
        
    def store_point_positions(self):
        for index in range(len(self.particles)):
            self.x_positions[index].append(self.particles[index].x)
            self.y_positions[index].append(self.particles[index].y)
    
    def make_plot(self):
        figure = pyplot.figure()
        plot = figure.add_subplot(1,1,1)
        
        for index in range(len(self.particles)):
            plot.plot(
                self.x_positions[index].value_in(nbody_system.length),
                self.y_positions[index].value_in(nbody_system.length)
            )
            
        figure.savefig("triple_{0:03}_{1:03}.png".format(self.rank, self.index))
            
    def stop(self):
        self.gravity_code.stop()

def new_initial_system():
    n = 10.0
    xi = 0.001
    result = Particles(3)
    mass =  numpy.random.rand(3)
    mass = mass / mass.sum()
    result.mass = mass * (1.0 | nbody_system.mass)
    a = (numpy.random.rand() * n) + (n+xi)
    b = numpy.random.rand() * 2*n
    c = (numpy.random.rand() * 2*n) + xi
    result[0].position = [0,0,0] | nbody_system.length
    result[1].position = [a,0,0] | nbody_system.length
    result[2].position = [b,c,0] | nbody_system.length
    result.radius = 0 | nbody_system.length
    result.velocity = [0,0,0] | nbody_system.speed
    return result

def calculate_escape_time(index, number_of_systems):
    numpy.random.seed()
    
    #from amuse.community.smallN.muse_dynamics_mpi import SmallN
    print "start of subprocess", index
    x0 = AdaptingVectorQuantity()
    y0 = AdaptingVectorQuantity()
    x1 = AdaptingVectorQuantity()
    y1 = AdaptingVectorQuantity()
    x2 = AdaptingVectorQuantity()
    y2 = AdaptingVectorQuantity()
    m0 = AdaptingVectorQuantity()
    m1 = AdaptingVectorQuantity()
    m2 = AdaptingVectorQuantity()
    tends = AdaptingVectorQuantity()
    
    try:
        for x in range(number_of_systems):
            gravity = Hermite()
            particles = new_initial_system()
            code = SimulateTripleSystemUntilDecay(gravity, particles, index, x)
            tend = code.evolve_model_until_escape()
            print index, x, tend
            code.stop()
            x0.append(particles[0].x)
            y0.append(particles[0].y)
            x1.append(particles[1].x)
            y1.append(particles[1].y)
            x2.append(particles[2].x)
            y2.append(particles[2].y)
            m0.append(particles[0].mass)
            m1.append(particles[1].mass)
            m2.append(particles[2].mass)
            tends.append(tend)
    except KeyboardInterrupt:
        pass
    
    print "end of subprocess", index
    
    return x0,y0,x1,y1,x2,y2,m0,m1,m2,tends
    
if __name__ == '__main__':
        
    result = AdaptingVectorQuantity()
    
    rank = MPI.COMM_WORLD.Get_rank()
    
    total_number_of_systems = int(sys.argv[1])
    data = calculate_escape_time(rank, total_number_of_systems)
   
    
    
    output = text.CsvFileText(filename = "triple_data_{0:03}.csv".format(rank))
    output.quantities = data
    output.attribute_names = ("x0","y0","x1","y1","x2","y2","m0","m1","m2","tend")
    output.store()
    
    times = AdaptingVectorQuantity()
    counts = AdaptingVectorQuantity()
    
    
    print "doing statistics"
    sorted_times = data[-1].sorted()
    t = 0 | nbody_system.time
    tend = 1000 | nbody_system.time
    dt = 1 | nbody_system.time
    while t < tend:
        i = 0
        while i < len(sorted_times):
            tn = sorted_times[i]
            i += 1
            if tn < t:
                continue
            else:
                break
                
        times.append(t)
        counts.append((len(sorted_times) - i) | units.none)
        t += dt
        
    output = text.CsvFileText(filename = "triple_statistics_{0:03}.csv".format(rank))
    output.quantities = (times,counts)
    output.attribute_names = ("t","n(t)")
    output.store()
    
    print "done"
"""

set logscale xy
"""
