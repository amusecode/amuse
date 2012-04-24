import numpy

from amuse.lab import *

from .interface import MikkolaInterface
from .interface import Mikkola

def new_system_of_Hulse_Taylor_pulsar():
    stars = core.Stars(2)
    Hulse = stars[0]
    Hulse.mass = 1.441 | units.MSun
    Hulse.radius = 1.4e-5 | units.RSun
    Hulse.position = [-1546691.3719943422, 0, 0] | units.km
    Hulse.velocity = [0.0, -110.0, 0.0] | units.km/units.s
    
    Taylor = stars[1]
    Taylor.mass = 1.387 | units.MSun
    Taylor.radius = 1.4e-5 | units.RSun
    Taylor.position = [1606908.6280056578, 0, 0] | units.km
    Taylor.velocity = [0.0, 114.28262436914201, 0.0] | units.km/units.s
    return stars

def HulseTaylor3():
    instance = Mikkola()
    stars = self.new_system_of_Hulse_Taylor_pulsar()
    instance.particles.add_particles(stars)
    Hulse = stars[0]
    Taylor = stars[1]
    
    postion_at_start = Taylor.position.value_in(units.AU)[0]
    
    #orbital period
    #period_HTpulsar = 7.751939106 | units.hour
    #period_HTpulsar = 77.51939106 | units.hour
    # period for abseidal motion
#        period_HTpulsar = 85.0 | units.yr #4.2degrees/year
    period_HTpulsar = 1.0 | units.yr 
    instance.evolve_model(period_HTpulsar)
    instance.particles.copy_values_of_state_attributes_to(stars)

    postion_after_full_rotation = Taylor.position.value_in(units.AU)[0]
   
    self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 4)
    
    instance.evolve_model(1.5 * period_HTpulsar)
    
    instance.particles.copy_values_of_state_attributes_to(stars)
    
    postion_after_half_a_rotation = Taylor.position.value_in(units.AU)[0]
    self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
    
    instance.evolve_model(1.75 * period_HTpulsar)
     
    instance.particles.copy_values_of_state_attributes_to(stars)
    
    postion_after_half_a_rotation = Taylor.position.value_in(units.AU)[1]
    
    self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
    instance.cleanup_code()
    del instance


def new_option_parser():
    result = OptionParser()
    result.add_option("-M", dest="M",type="float",default=12.)
    result.add_option("-m", dest="m",type="float",default=10.)
    result.add_option("-a", dest="a",type="float",default=205.)
    result.add_option("-e", dest="e",type="float",default=0.6)
    result.add_option("-t", dest="tend",type="float",default=10.)
    result.add_option("-d", dest="dtdiag",type="float",default=c.very_long_time)
    return result

def main(a=205, e=0.6, M=13., m=11., tend=10., dtdiag=c.very_long_time) :

    bs = double_star(a, e, M, m)
    M = bs.p.get_mass()
    m = bs.s.get_mass()
    Porb = bs.OrbitalPeriod()
    bs.evolve(tend, dtdiag)
#    print "Th binary Orbial Velocity: ", bs.OrbitalVelocity()
    print "BINARY:", bs
    bs.terminate()

if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)
