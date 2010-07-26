




from amuse.lab import *


import subprocess
import os
import numpy
import time
import sys
import signal

from mpi4py import MPI

class TimeoutException(Exception):
    pass
    
class SkipException(Exception):
    pass
    
class RunSpeedTests(object):
    
    @late
    def report_lines(self):
        return []
        
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
        
    def __init__(self, total_number_of_points):
        self.total_number_of_points = total_number_of_points
    
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
            method = getattr(self, x)
            print >> sys.stderr , self.row_formatters[0](method), '...',
            try:        
                method()
            except TimeoutException, ex:
                print >> sys.stderr , "timed out,", ex
                continue
            except SkipException, ex:
                print >> sys.stderr , "skipped,", ex
                continue
            except Exception, ex:
                print ex
                self.t1=-1
                self.t0=0
                pass
            delta_time = self.t1-self.t0
            self.total_time = self.total_time + delta_time
            print >> sys.stderr , self.row_formatters[1](delta_time)
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
        new_plummer_sphere(self.total_number_of_points)
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
        input = new_plummer_sphere(self.total_number_of_points)
        if self.total_number_of_points > 20000:
            raise SkipException("too many points")
        self.start_measurement()
        input.scale_to_standard()
        self.end_measurement()
        
        
    def speed_calculate_potential_energy(self):
        
        if self.total_number_of_points > 20000:
            raise SkipException("too many points")
    
        input = new_plummer_sphere(self.total_number_of_points)
        self.start_measurement()
        input.potential_energy()
        self.end_measurement()
        
        
    def speed_calculate_kinetic_energy(self):
        input = new_plummer_sphere(self.total_number_of_points)
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
        particles = new_plummer_sphere(self.total_number_of_points)
        particles.radius = 0| nbody.length
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
        
    
    def speed_add_particles_to_code_SI(self):
        """plummer sphere"""
        converter = nbody.nbody_to_si(1 | units.parsec, self.total_number_of_points | units.MSun)
        code = BHTree(converter)
        particles = new_plummer_sphere(self.total_number_of_points, converter)
        particles.radius = 0| units.RSun
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
    
    
    def speed_copy_particles_from_code(self):
        code = BHTree()
        particles = new_plummer_sphere(self.total_number_of_points)
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
        particles = new_plummer_sphere(self.total_number_of_points)
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
        particles = new_plummer_sphere(self.total_number_of_points)
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
    
    
if __name__ == '__main__':
    #channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    else:
        n = 100
        
    x = RunSpeedTests(n)
    x.run()







