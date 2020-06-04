"""
Runs several tests to determine how fast common operations
are in AMUSE.

to profile (in the amuse root directory):

./amuse.sh -m cProfile -s cumulative test/reports/speed_report.py --n_order==4 --test==speed_copy_to_set > profile.txt

"""


#import numpypy


from amuse.lab import *
#from amuse.datamodel import *
#from amuse.units import nbody_system
from amuse.support.thirdparty import texttable

import traceback

import subprocess
import os
import numpy
import time
import sys
import signal

from optparse import OptionParser

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
        return ('action', 'duration\n(seconds)')
        
    @late
    def row_formatters(self):
        return (self.method_to_action, lambda x : x)
    
    @late
    def maximum_number_of_seconds_to_allow_per_test(self):
        return 1200
        
    def method_to_action(self, x):
        name = x.__name__
        name = name[len(self.method_prefix):]
        name = name.replace('_d_', '.')
        name = name.replace('_', ' ')
        return name
        
    def __init__(self, 
            total_number_of_points, 
            subset_number_of_points = -1, 
            name_of_the_method = None,
            include_code_tests = False,
            include_slow_tests = False,
            csv_output = False):
        self.total_number_of_points = total_number_of_points
        if subset_number_of_points < 0:
            subset_number_of_points = self.total_number_of_points
            
        self.subset_number_of_points = subset_number_of_points
        self.name_of_the_method = name_of_the_method
        self.include_code_tests = include_code_tests
        self.include_slow_tests = include_slow_tests
        self.csv_output = csv_output
    
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
            print(self.row_formatters[0](method), '...', end=' ', file=sys.stderr)
            try:        
                method()
            except TimeoutException as ex:
                print("timed out,", ex, file=sys.stderr)
                continue
            except SkipException as ex:
                print("skipped,", ex, file=sys.stderr)
                continue
            except Exception as ex:
                print(ex)
                traceback.print_exc()
                self.t1=-1
                self.t0=0
                pass
            delta_time = self.t1-self.t0
            self.total_time = self.total_time + delta_time
            print(self.row_formatters[1](delta_time), file=sys.stderr)
            self.report_lines.append((method,delta_time))
        
        if self.csv_output:
            self.output_csv_line()
        else:
            self.output_table()
    
    def output_csv_line(self):
        line = []
        line.append('N')
        line.append('M')
        for x in self.report_lines:
            line.append(self.method_to_action(x[0]))
        print('#' + ','.join(line))
        
        line = []
        line.append(str(self.total_number_of_points))
        line.append(str(self.subset_number_of_points))
        for x in self.report_lines:
            line.append(str(x[1]))
        print(','.join(line))
        
    def output_table(self):
        lines = []
        
        lines.append(':run date:')
        lines.append('    {0}'.format(time.asctime()))
        lines.append('')
        lines.append(':number of points:')
        lines.append('    {0}'.format(self.total_number_of_points))
        lines.append('')
        
        for x in lines:
            print(x)
            
        table = texttable.Texttable()
        #table.set_deco(texttable.Texttable.HEADER)
        table.set_cols_dtype([
            't',  # text 
            'f',  # float (decimal)
        ])
        table.set_cols_align(["l", "r"])
        rows = []
        rows.append(self.header_line)
        rows.extend(self.report_lines_as_strings)
        table.add_rows(rows)
        print(table.draw())
        
            
    def names_of_testing_methods(self):
        for x in sorted(dir(type(self))):
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
            yield list(map(lambda x, formatter : formatter(x), line, self.row_formatters))
        
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
        
        self.is_slow_test()
        
        self.start_measurement()
        input.scale_to_standard()
        self.end_measurement()
        
    def is_slow_test(self):
        if not self.include_slow_tests:
            raise SkipException("slow tests disabled")
            
    def is_single_particle_test(self):
        if not self.include_slow_tests:
            raise SkipException("single particle tests disabled")
            
    def speed_calculate_potential_energy(self):
        
        self.is_slow_test()
        input = new_plummer_model(self.total_number_of_points)
        self.start_measurement()
        input.potential_energy(G=nbody_system.G)
        self.end_measurement()
        
        
    def speed_calculate_kinetic_energy(self):
        input = new_plummer_model(self.total_number_of_points)
        self.start_measurement()
        input.kinetic_energy()
        self.end_measurement()
        
    
    def speed_start_and_stop_BHTree_code(self):
        self.is_code_test()
        self.start_measurement()
        code = BHTree()
        code.stop()
        self.end_measurement()
    
    def speed_add_particles_to_code(self):
        self.is_code_test()
        
        code = BHTree()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
        
    
    def speed_add_particles_to_code_SI(self):
        """plummer sphere"""
        self.is_code_test()
        
        converter = nbody.nbody_to_si(1 | units.parsec, self.total_number_of_points | units.MSun)
        code = BHTree(converter)
        particles = new_plummer_model(self.total_number_of_points, converter)
        particles.radius = 0| units.RSun
        self.start_measurement()
        code.particles.add_particles(particles)
        self.end_measurement()
        code.stop()
    
    
    def speed_copy_particles_from_code(self):
        self.is_code_test()
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
    
        self.is_code_test()
        self.is_slow_test()
    
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
    
        self.is_code_test()
        self.is_slow_test()
    
        code = Hermite()
        particles = new_plummer_model(self.total_number_of_points)
        particles.radius = 0| nbody.length
        code.particles.add_particles(particles)
        self.start_measurement()
        code.evolve_model(0.001 | nbody.time)
        self.end_measurement()
        code.stop()
       
    def is_code_test(self):
        if not self.include_code_tests:
            raise SkipException("code tests disabled") 

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
        
        self.is_single_particle_test()
        
        particles = Particles(self.total_number_of_points)
        particles.radius = 1.0 | nbody_system.length
        array = numpy.zeros(self.total_number_of_points, dtype=numpy.object)
        self.start_measurement()
        i = 0
        for x in particles:
            array[i] = x
            i += 1
        self.end_measurement()
    
    
    def speed_iterate_over_particles2(self):
        
        self.is_single_particle_test()
        
        class A(object):
            __slots__ = ('i', 'j', 'k' , 'l')
            __array_interface__ = {'shape':()}
            
            def __len__(self):
                raise AttributeError()
            def __iter__(self):
                raise AttributeError()
    
            
            def __init__(self, i, j = 10, k = 20, l = 24):
                self.i = i
                if i > 10:
                    self.j = j + 10
                else:
                    self.j = j
                
                self.k = k
                self.l = l
        
            def __getattr__(self, name_of_the_attribute):
                raise AttributeError("You tried to access attribute '{0}' but this attribute is not defined for this set.".format(name_of_the_attribute, ex))
    
                
        array = numpy.zeros(self.total_number_of_points, dtype=numpy.object)
        
        self.start_measurement()
        for x in range(self.total_number_of_points):
            array[x] = A(x,x,x,x)
        self.end_measurement()
    

    def speed_iterate_over_array(self):
        self.is_single_particle_test()
        
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        self.start_measurement()
        for x in particles:
            x.radius
        self.end_measurement()
    
    

    def speed_create_N_particles(self):
        self.is_single_particle_test()
        
        particles = Particles(self.total_number_of_points)
        particles.radius = 1.0 | nbody_system.length
        self.start_measurement()
        for x in range(self.total_number_of_points):
            Particle(x)
        self.end_measurement()
    
    def speed_copy_attributes_from_code(self):
        self.is_code_test()
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
        self.is_code_test()
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
        self.is_code_test()
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
        self.is_code_test()
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
        
    def speed_copy(self):
        particles_all = new_plummer_model(self.total_number_of_points)
        
        self.start_measurement()
        particles_all.copy_to_memory()
        self.end_measurement()
        
        
    def speed_copy_subset(self):
        particles_all = new_plummer_model(self.total_number_of_points)
        subset = particles_all[:self.subset_number_of_points]
        self.start_measurement()
        subset.copy_to_memory()
        self.end_measurement()
        
    
    def speed_select_array(self):
        particles_all = new_plummer_model(self.total_number_of_points)
        
        self.start_measurement()
        particles_selected = particles_all.select_array(lambda position: position.lengths() > 0.5 | nbody_system.length, ["position"])
        particles_selected.x
        self.end_measurement()
        
    def speed_select_with_get_item(self):
        particles_all = new_plummer_model(self.total_number_of_points)
        
        self.start_measurement()
        particles_selected = particles_all[particles_all.position.lengths() > 0.5 | nbody_system.length]
        self.end_measurement()

    def speed_iterate_over_quantity(self):
        
        lengths = numpy.arange(self.total_number_of_points) | nbody_system.length
        self.start_measurement()
        for x in range(self.total_number_of_points):
            lengths[x]
        self.end_measurement()
        
    def speed_add_particles(self):
        particles_to_add = new_plummer_model(self.total_number_of_points)
        step = self.total_number_of_points / 10
        particle_sets = []
        i = 0
        while i < self.total_number_of_points:
            j = i + step 
            if j > self.total_number_of_points:
                j = self.total_number_of_points
            particle_sets.append(particles_to_add[i:j].copy())
            i = j
            
        particles = Particles()
        self.start_measurement()
        for x in particle_sets:
            particles.add_particles(x)
        self.end_measurement()
        
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n", 
        default = 2,
        dest="total_number_of_points",
        help="lenght of particle set",
        type="int"
    )
    result.add_option(
        "-m", 
        default = -1,
        dest="subset_number_of_points",
        help="lenght of particles in subset functions",
        type="int"
    )
    result.add_option(
        "--test", 
        default = None,
        dest="name_of_the_method",
        help="name of the test method to run",
        type="string"
    )
    result.add_option(
        "--code", 
        action="store_true",
        default=False,
        dest="include_code_tests",
        help="also run tests with codes"

    )
    result.add_option(
        "--slow", 
        action="store_true",
        default=False,
        dest="include_slow_tests",
        help="also run slow tests, scale N**2 or worse, for example calculating the potential energy"
    )
    result.add_option(
        "--csv", 
        action="store_true",
        default=False,
        help="display run as one comma separated line",
        dest="csv_output"

    )
    
    return result
    
    
    
if __name__ == '__main__':
    options, arguments = new_option_parser().parse_args()
    x = RunSpeedTests(**options.__dict__)
    x.run()








