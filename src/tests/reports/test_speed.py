from amuse.support.testing.amusetest import get_path_to_results
from reports.interface import TestCode
import os
from pathlib import Path
import numpy
import time

# from amuse.rfi import channel
from amuse.rfi.core import *


class RunSpeedTests(object):
    
    def __init__(self):
        self.exefile = str(Path(__file__).parent / 'c_worker')
        self.number_of_gridpoints = [8]            
    
    def start(self):
        
        for number_of_points_in_one_dimension in self.number_of_gridpoints:
            result = self.run(number_of_points_in_one_dimension)
    
            print(', '.join([str(x) for x in result]))
                
    def run(self, number_of_points_in_one_dimension):
    
        instance = TestCode(self.exefile)

        total_number_of_points = number_of_points_in_one_dimension ** 3
        number_of_bytes = 4 + 8 + 8 + 8
        total_number_of_bytes = total_number_of_points * (number_of_bytes + 4)
        indices = numpy.array(range(total_number_of_points), dtype='int32')

        data_x = numpy.array(range(total_number_of_points), dtype='float64')
        data_y = numpy.array(range(total_number_of_points), dtype='float64')
        data_z = numpy.array(range(total_number_of_points), dtype='float64')

        errorcode = instance.set_number_of_points_in_one_dimension(number_of_points_in_one_dimension)
        if errorcode < 0:
            raise Exception("Could not allocate memory")

        t0 = time.time()
        instance.set_data(indices, data_x, data_y, data_z)
        t1 = time.time()
        dt = t1 - t0
        mbytes_per_second = total_number_of_bytes / dt / (1000.0 * 1000.0)

        t2 = time.time()
        instance.set_data_to_same(total_number_of_points, 0.0, 1.0, 2.0)
        t3 = time.time()
        
        instance.reset()
        instance.stop()
        
        return dt, total_number_of_points, mbytes_per_second, t3-t2, (dt - (t3-t2)) / (t3-t2)     
        
        
def test_speed():
    x = RunSpeedTests()
    x.number_of_gridpoints = [8]
    x.start()


if __name__ == '__main__':
    #channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD
    x = RunSpeedTests()
    x.number_of_gridpoints = [64, 128, 192]
    x.start()
