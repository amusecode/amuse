import numpy

from amuse.test import amusetest
from amuse.units import units, nbody_system, constants

from amuse.ext.rotating_bridge import Rotating_Bridge, inertial_to_rotating, rotating_to_inertial

from amuse.ext.composition_methods import *

from amuse.datamodel import Particles

class drift_without_gravity(object):
    """
    This class evolves the motion of the test particles (no gravity) 
    """
    def __init__(self, particles, time= 0 |units.Myr):
        self.particles = particles
        self.model_time = time
    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = t_end
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()

class TestRotatingBridge(amusetest.TestCase):
    def test1(self):
        p0=Particles(1)
        p0.position = [0.0, 1.0, 0.0]
        p0.velocity = [1.0, 0.0, 0.0]
        
        pos1=p0.position.copy()
        vel1=p0.velocity.copy()
        
        omega=1.
        N=0.25
        dt=0.01
        method = SPLIT_6TH_SS_M13

        tend=N*2*numpy.pi/omega
        dt=dt*2*numpy.pi/omega

        pr=inertial_to_rotating(0,omega, p0)
        
        drift=drift_without_gravity(pr, time=0)
        sys=Rotating_Bridge(omega, timestep=dt, method=method)
        sys.add_system(drift)
        
        sys.evolve_model(tend)
        
        pi=rotating_to_inertial(tend, omega,pr)
        
        pos2=pi.position.copy()
        
        
        self.assertAlmostEqual(sys.model_time,tend)
        self.assertAlmostEqual(pos1+vel1*tend,pos2,12)
        
        pi.velocity=-pi.velocity
        
        pr=inertial_to_rotating(0,-omega,pi)
        drift=drift_without_gravity(pr,time=0)
        sys=Rotating_Bridge(-omega, timestep=dt,method=method)
        sys.add_system(drift)
        
        sys.evolve_model(tend)        

        pi=rotating_to_inertial(tend,-omega,pr)
        
        pos3=pi.position.copy()
        
        self.assertAlmostEqual(pos1,pos3,12)

    def test2(self):
        p0=Particles(1)
        p0.position = [0.0, 1.0, 0.0] | units.m
        p0.velocity = [1.0, 0.0, 0.0] | units.m/units.s
        
        pos1=p0.position.copy()
        vel1=p0.velocity.copy()
        
        omega=1. | units.s**-1
        N=0.25
        dt=0.01
        method = SPLIT_6TH_SS_M13

        tend=N*2*numpy.pi/omega
        dt=dt*2*numpy.pi/omega

        pr=inertial_to_rotating(0*dt,omega, p0)
        
        drift=drift_without_gravity(pr, time=0 | units.s)
        sys=Rotating_Bridge(omega, timestep=dt, method=method)
        sys.add_system(drift)
        
        sys.evolve_model(tend)
        
        pi=rotating_to_inertial(tend, omega,pr)
        
        pos2=pi.position.copy()
        
        
        self.assertAlmostEqual(sys.model_time,tend)
        self.assertAlmostEqual(pos1+vel1*tend,pos2,12)
        
        pi.velocity=-pi.velocity
        
        pr=inertial_to_rotating(0*dt,-omega,pi)
        drift=drift_without_gravity(pr,time=0| units.s)
        sys=Rotating_Bridge(-omega, timestep=dt,method=method)
        sys.add_system(drift)
        
        sys.evolve_model(tend)        

        pi=rotating_to_inertial(tend,-omega,pr)
        
        pos3=pi.position.copy()
        
        self.assertAlmostEqual(pos1,pos3,12)
