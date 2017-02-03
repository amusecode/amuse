import Queue
import threading
import multiprocessing
from amuse.lab import *

m = Queue.Queue()
z = Queue.Queue()
t = Queue.Queue()
c = Queue.Queue()

def remote_worker_code():
    while True:
        mass = m.get()
        metalicity = z.get()
        time = t.get()
        code = c.get()
        evolve_single_star(mass, metalicity, time, code)
        m.task_done()
        z.task_done()
        t.task_done()
        c.task_done()

def evolve_with_different_stellar_model(mass, metalicity, time, codes):
    for ci, zi in zip(codes, metalicity):
        m.put(mass)
        z.put(zi)
        t.put(time)
        c.put(ci)
    n_cpu = multiprocessing.cpu_count() #detect number of cores
    print "Number of CPUs:", n_cpu
    for i in range(n_cpu):
        th = threading.Thread(target=remote_worker_code)
        th.daemon = True
        th.start()
    m.join() #block until all tasks are done
    z.join() #block until all tasks are done
    t.join() #block until all tasks are done
    c.join() #block until all tasks are done

def evolve_single_star(mass, z, time, code):
    stars = Particles(mass=mass)
    stellar = code()
    stellar.parameters.metallicity = z
    stellar.particles.add_particles(stars)
    channel = stellar.particles.new_channel_to(stars)

    stellar.evolve_model(time)
    channel.copy()
    print "Final star at time=", stellar.model_time, "(z=", z, ") M=", stars.mass, "R=", stars.radius, "type=", stellar.particles.stellar_type
    stellar.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", action="store_true",
                      dest="threadded", help="run threadded [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    set_printing_strategy("custom",\
        preferred_units = [units.MSun, units.RSun, units.Myr],\
        precision = 6, prefix = "", separator = "[", suffix = "]")
    mass = 10|units.MSun
    metalicities = [0.01, 0.02, 0.03]
    time = 1|units.Myr
    codes = [SeBa, MESA, SSE]
    if o.threadded:
        evolve_with_different_stellar_model(mass, metalicities, time, codes)
    else:
        print "Run sequentially"
        for ci, zi in zip(codes, metalicities):
            evolve_single_star(mass, zi, time, ci)
