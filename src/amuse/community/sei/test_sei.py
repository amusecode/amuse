from amuse.test.amusetest import TestWithMPI
from amuse.units import nbody_system
from amuse.units import units
import os
import sys
import numpy
import math

from interface import SeiInterface,Sei

from amuse.support import data
class TestSeinterface(TestWithMPI):

    def test0(self):
        instance = SeiInterface()
        instance.initialization()
        instance.set_state(0,1,0,0,0,0,0)
        for i in range(0,10):
            instance.evolve(i)
            print instance.get_state(0)
        instance.stop()
    
class TestSei(TestWithMPI):

    def test0(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1 | units.AU)
        
        particle = data.Particles(1)
        particle.position = [1.0, 0.0, 0.0,]|units.AU
        particle.velocity = [0.0, 2.0*3.1415926535*1.0/365, 0.0] | units.AUd
        sei = Sei(convert_nbody)
        sei.initialization()
        sei.particles.add_particles(particle)
        print sei.particles.position.x.value_in(units.AU)
        for i in range(365):
            sei.evolve_model(i|units.day)
            print sei.particles.position.x.value_in(units.AU)
