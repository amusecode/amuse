import numpy
import sys
import time

from amuse.community.hermite0.interface import HermiteInterface

try:
    from matplotlib import pyplot
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False

def calculate_speed(range_of_number_of_particles):
    result = []
    for n in range_of_number_of_particles: #range(8000,20000, 1000):
        hermite1 = HermiteInterface()
        hermite1.initialize_code()
        
        hermite2 = HermiteInterface()
        hermite2.initialize_code()
        
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
    
        t0 = time.time()        
        hermite1.new_particle(values
                              , values
                              , values
                              , values
                              , values
                              , values
                              , values
                              , values)
        t1 = time.time()
        d1 = t1 - t0
        #print d1, t1, t0
        
        t0 = time.time() 
        for i in range(n-1):
            hermite2.new_particle(values[i]
                                  , values[i]
                                  , values[i]
                                  , values[i]
                                  , values[i]
                                  , values[i]
                                  , values[i]
                                  , values[i])
        t1 = time.time()
        d2 = t1 - t0
        result.append((n, d1, d2, d2/d1))
                    
        hermite1.cleanup_code()
        hermite2.cleanup_code()
        
        del hermite1
        del hermite2
    return result
    
if __name__ == '__main__':
    measurements = calculate_speed([
        2500, 5000, 7500, 
        10000, 15000, 20000, 
        25000, 30000, 40000])
    
    figure = pyplot.figure(figsize = (10,10))
    plot = figure.add_subplot(1,1,1)
    
    for number, dt0, dt1, ratio in measurements:
        color = 'b'
        plot.plot([number],[ratio], color + 'o')
        
    plot.set_xlim(0.0, 40000)
    plot.set_ylim(0.0, 350)
    
    figure.savefig("speed.svg")    
    
    print(measurements)
