from amuse.support.data import core
from amuse.support.io import nemotsf
from amuse.support.units import units

import unittest
import os.path

class Test(unittest.TestCase):
    
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
        I = nemotsf.Tsf2Particles()
        I.convert_to_particles(os.path.join(directory, 'p10.txt'))
        J = nemotsf.Particles2Tsf()
        J.convert_to_tsf(I.Particles)
        J.init_string()
