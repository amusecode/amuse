"""
Runs several tests to determine how fast common operations
are in AMUSE.

to profile (in the amuse root directory):

./amuse.sh -m cProfile -s cumulative test/reports/speed_report.py 10000 speed_copy_to_set > profile.txt

"""





from amuse.lab import *
import traceback

import subprocess
import os
import numpy
import time
import sys
import signal


from mpi4py import MPI

from amuse.datamodel import ParticlesSuperset
class TimeoutException(Exception):
    pass
    
class SkipException(Exception):
    pass
    
class RunSpeedTests(object):
    
    @late
    def report_lines(self):
        return []
                
    @late
    def header_line(self):
        return ('action', 'duration (seconds)')
        
    @late
    def row_formatters(self):
        return (self.method_to_action, '{0:0.3f}'.format)
    
    @late
    def maximum_number_of_seconds_to_allow_per_test(self):
        return 1200
        
    def method_to_action(self, x):
        name = x.__name__
        name = name[len(self.method_prefix):]
        name = name.replace('_d_', '.')
        name = name.replace('_', ' ')
        return name
        
    def __init__(self, total_number_of_points, name_of_the_method = None):
        self.total_number_of_points = total_number_of_points
        self.name_of_the_method = name_of_the_method
    
    def handle_timeout(self, signum, frame):
        self.t1=-1
        self.t0=0
        raise TimeoutException("Test did not finish in allocated time frame")

        
    def start_measurement(self):
        self.t0 = time.time()
        signal.setitimer(signal.ITIMER_REAL, self.maximum_number_of_seconds_to_allow_per_test)
        signal.signal(signal.SIGALRM, self.handle_timeout)

        
    def end_measurement(self):
        self.t1 = time.time()
        signal.setitimer(signal.ITIMER_REAL, 0)
        
    def run(self):
        
        self.total_time = 0.0
    
        for x in self.names_of_testing_methods():
            if not self.name_of_the_method is None:
                if x != self.name_of_the_method:
                    continue
            method = getattr(self, x)
            print >> sys.stderr, self.row_formatters[0](method), '...',
            try:        
                method()
            except TimeoutException, ex:
                print >> sys.stderr, "timed out,", ex
                continue
            except SkipException, ex:
                print >> sys.stderr, "skipped,", ex
                continue
            except Exception, ex:
                print ex
                traceback.print_exc()
                self.t1=-1
                self.t0=0
                pass
            delta_time = self.t1-self.t0
            self.total_time = self.total_time + delta_time
            print >> sys.stderr, self.row_formatters[1](delta_time)
            self.report_lines.append((method,delta_time))
        
        lines = []
        
        lines.append(':run date:')
        lines.append('    {0}'.format(time.asctime()))
        lines.append('')
        lines.append(':number of points:')
        lines.append('    {0}'.format(self.total_number_of_points))
        lines.append('')
        lines.append(self.grid_table_row_separator_line())
        lines.append(self.grid_table_row_line(self.header_line))
        lines.append(self.grid_table_row_separator_line('='))
        for x in self.report_lines_as_strings[:-1]:
            lines.append(self.grid_table_row_line(x))
        lines.append(self.grid_table_row_line(self.report_lines_as_strings[-1]))
        lines.append(self.grid_table_row_separator_line('='))
        lines.append(self.grid_table_row_line(['total time', self.row_formatters[1](self.total_time) ]))
        lines.append(self.grid_table_row_separator_line('='))
        
        for x in lines:
            print x
            
    def names_of_testing_methods(self):
        for x in dir(type(self)):
            if x.startswith(self.method_prefix):
                yield x
    
    @late
    def method_prefix(self):
        return "speed_"
        
    @late
    def report_lines_as_strings(self):
        return list(self.iter_report_lines_as_strings())
    
    def iter_report_lines_as_strings(self):
        for line in self.report_lines:
            yield map(lambda x, formatter : formatter(x), line, self.row_formatters)
        
    def speed_make_plummer_sphere(self):
        """plummer sphere"""
        self.start_measurement()
        new_plummer_model(self.total_number_of_points)
        self.end_measurement()
        
    def speed_make_salpeter_mass_distribution(self):
        """plummer sphere"""
        self.start_measurement()
        new_salpeter_mass_distribution(self.total_number_of_points)
        self.end_measurement()
        
    def speed_make_salpeter_mass_distribution_nbody(self):
        """plummer sphere"""
        self.start_measurement()
        new_salpeter_mass_distribution_nbody(self.total_number_of_points)
        self.end_measurement()
        
    def speed_scale_plummer_sphere(self):
        input = new_plummer_model(self.total_number_of_points)
        if self.total_number_of_points > 20000:
            raise SkipException("too many points")
        self.start_measurement()
        input.scale_to_standard()
        self.end_measurement()
        
        
    def speed_calculate_potential_energy(self):
        
        if self.total_number_of_points > 20000:
            raise SkipException("too many points")
    
        input = new_plummer_model(self.total_number_of_points)
        self.start_measurement()
        input.potential_energy()
        self.end_measurement()
        
        
    def speed_calculate_kinetic_energy(self):
        input = new_plummer_model(self.total_number_of_points)
        self.start_measurement()
        input.kinetic_energy()
        self.end_measurement()
        
    
    def speed_start_and_stop_BHTree_code(self):
        self.start_measurement()
        code = BHTree()
        code.stop()
        self.end_measurement()
    
    def speed_add_particles_to_code(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
        
    
    def speed_add_particles_to_code_SI(self):
        """plummer sphere"""
        converter = nbody.nbody_to_si(1 | units.parsec, self.total_number_of_points | units.MSun)
        code = BHTree(converter)
        particles = new_plummer_model(self.total_number_of_points, converter)
        particles.radius = 0| units.RSun
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
    
    
    def speed_copy_particles_from_code(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        code.particles.add_particles(particles)
        channel = code.particles.new_channel_to(particles)
        self.start_measurement()
        channel.copy()
        self.end_measurement()
        code.stop()
    
    def speed_evolve_code_0_d_001_time_in_BHTree(self):
        """plummer sphere"""
    
        if self.total_number_of_points > 10000:
            raise SkipException("too many points")
    
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        code.particles.add_particles(particles)
        self.start_measurement()
        code.evolve_model(0.001 | nbody.time)
        self.end_measurement()
        code.stop()
    
    def speed_evolve_code_0_d_001_time_in_Hermite(self):
        """plummer sphere"""
    
        if self.total_number_of_points > 5000:
            raise SkipException("too many points")
    
        code = Hermite()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        code.particles.add_particles(particles)
        self.start_measurement()
        code.evolve_model(0.001 | nbody.time)
        self.end_measurement()
        code.stop()
        

    @late
    def maximum_column_widths(self):
        maximums = [len(x) for x in self.header_line]
        for row in self.report_lines_as_strings:
            maximums = [max(len(x), y) for x,y in zip(row, maximums)]
        return maximums
    
    @late
    def column_widths(self):
        return [x + 2 for x in self.maximum_column_widths]
        
    def grid_table_row_separator_line(self, line_character = '-' ):
        parts = []
        for x in self.column_widths:
            parts.append('+')
            parts.append(line_character * x)
        parts.append('+')
        return ''.join(parts)
    
    
    def grid_table_row_line(self, row):
        parts = []
        for width, x in zip(self.column_widths, row):
            parts.append('|')
            parts.append(' ')
            parts.append(x.rjust(width-1))
        parts.append('|')
        return ''.join(parts)
    
    

    def speed_iterate_over_particles(self):
        particles = Particles(self.total_number_of_points)
        particles.radius = 1.0 | nbody_system.length
        self.start_measurement()
        for x in particles:
            x.radius
        self.end_measurement()
    
    

    def speed_iterate_over_array(self):
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        self.start_measurement()
        for x in particles:
            x.radius
        self.end_measurement()
    
    

    def speed_iterate_over_particles2(self):
        particles = Particles(self.total_number_of_points)
        particles.radius = 1.0 | nbody_system.length
        self.start_measurement()
        y = particles.radius
        for x in range(self.total_number_of_points):
            Particle(x)
        self.end_measurement()
    
    def speed_copy_attributes_from_code(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        code.particles.add_particles(particles)
        channel = code.particles.new_channel_to(particles)
        self.start_measurement()
        channel.copy()
        self.end_measurement()
        code.stop()
        
    def speed_copy_attributes_from_code_to_empty(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        empty_particles = particles.empty_copy()
        
        code.particles.add_particles(particles)
        channel = code.particles.new_channel_to(empty_particles)
        self.start_measurement()
        channel.copy()
        self.end_measurement()
        code.stop()
        
    def speed_copy_mass_attribute_from_code_to_empty(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        empty_particles = particles.empty_copy()
        
        code.particles.add_particles(particles)
        channel = code.particles.new_channel_to(empty_particles)
        self.start_measurement()
        channel.copy_attribute("mass","zmass")
        self.end_measurement()
        code.stop()
    
    def speed_copy_position_and_velocity_attributes_from_code_to_empty(self):
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        empty_particles = particles.empty_copy()
        
        code.particles.add_particles(particles)
        channel1 = code.particles.new_channel_to(empty_particles)
        channel2 = empty_particles.new_channel_to(code.particles)
        self.start_measurement()
        channel1.copy_attributes(["x","y","z", "vx","vy","vz"])
        channel2.copy_attributes(["x","y","z"])
        self.end_measurement()
        code.stop()
        
    
    def speed_copy_to_superset(self):
        particles1 = new_plummer_model(self.total_number_of_points)
        particles2 = new_plummer_model(self.total_number_of_points)
        particles_all = ParticlesSuperset([particles1, particles2])
        
        empty_particles = particles_all.empty_copy()
        
        channel1 = particles_all.new_channel_to(empty_particles)
        channel2 = empty_particles.new_channel_to(particles_all)
        self.start_measurement()
        channel1.copy_attributes(["x","y","z"])
        channel2.copy_attributes(["x","y","z"])
        self.end_measurement()
    
    def speed_copy_to_set(self):
        particles_all = new_plummer_model(self.total_number_of_points * 2)
        empty_particles = particles_all.empty_copy()
        
        channel1 = particles_all.new_channel_to(empty_particles)
        channel2 = empty_particles.new_channel_to(particles_all)
        self.start_measurement()
        channel1.copy_attributes(["x","y","z"])
        channel2.copy_attributes(["x","y","z"])
        self.end_measurement()
        
        
    
    
if __name__ == '__main__':
    #channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    else:
        n = 100
    
    if len(sys.argv) > 2:
        name_of_the_method = sys.argv[2]
    else:
        name_of_the_method = None
        
    x = RunSpeedTests(n, name_of_the_method)
    x.run()








