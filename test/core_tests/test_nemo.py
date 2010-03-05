from amuse.support.data import core
from amuse.support.io import nemotsf
from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import values


#import unittest
from amuse.test import amusetest
import os.path
import numpy 

class Test(amusetest.TestCase):
    
    def Ntest1(self):
        directory = os.path.dirname(__file__)
        
        I = nemotsf.Tsf2Particles()
        I.convert_to_particles(os.path.join(directory, 'p10.txt'))
                               
        self.assertEquals(I.number_of_particles, 10)
        self.assertEquals(I.Particles.mass[0], 0.1|units.MSun)

    def Ntest2(self):
        directory = os.path.dirname(__file__)
        
        I = nemotsf.Tsf2Particles()
        I.convert_to_particles(os.path.join(directory, 'h2048.txt'))
                               
        self.assertEquals(I.number_of_particles, 2048)
        self.assertEquals(I.Particles.mass[0], 0.000488281|units.MSun)

    def test3(self):
        directory = os.path.dirname(__file__)
        convert_nbody = nbody_system.nbody_to_si(1|units.g, 1|units.m)
        I = nemotsf.Tsf2Particles()
        I.convert_to_particles(os.path.join(directory, 'p10.txt'), convert_nbody)
        self.assertAlmostEqual(I.Particles.mass, values.new_quantity(0.1*numpy.ones(10), units.g), constants.precision)

    def test4(self):

        directory = os.path.dirname(__file__)
        I = nemotsf.Tsf2Particles()
        I.convert_to_particles(os.path.join(directory, 'p10.txt'))
        self.assertAlmostEqual(I.Particles.mass[0].value_in(nbody_system.mass), 0.1,  constants.precision)

    def test5(self):
        pass
        #J = nemotsf.Particles2Tsf()
        #J.convert_to_tsf(I.Particles)
        #J.init_string()
