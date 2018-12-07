import Queue
import threading
import multiprocessing
from amuse.lab import *

###BOOKLISTSTART1###
code_queue = Queue.Queue()

def remote_worker_code():
    code = code_queue.get()
    evolve_single_star(code)
    code_queue.task_done()

def evolve_with_different_stellar_model(codes):
    for ci in codes:
        code_queue.put(ci)
    n_cpu = multiprocessing.cpu_count() 
    for i in range(n_cpu):
        th = threading.Thread(target=remote_worker_code)
        th.daemon = True
        th.start()
    code_queue.join() 	    # block until all tasks are done
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def evolve_single_star(code):
    stars = Particles(mass=10|units.MSun)
    stellar = code()
    stellar.particles.add_particles(stars)
    channel = stellar.particles.new_channel_to(stars)

    stellar.evolve_model(1|units.Myr)
    channel.copy()
    print "Star evolved to time=", stellar.model_time, \
          " M=", stars.mass, "R=", stars.radius
    stellar.stop()
###BOOKLISTSTOP2###
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", action="store_true",
                      dest="threaded", help="run threaded [%default]")
    return result

if __name__ in ('__main__', '__plot__'):

    o, arguments  = new_option_parser().parse_args()
    set_printing_strategy("custom",\
        preferred_units = [units.MSun, units.RSun, units.Myr],\
        precision = 6, prefix = "", separator = "[", suffix = "]")

    codes = [SeBa, MESA, SSE, EVtwin]
    
    if o.threaded:
        print "Run threaded"
        evolve_with_different_stellar_model(codes)
    else:
        print "Run sequentially"
        for ci in codes:
            evolve_single_star(ci)
