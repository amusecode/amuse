from amuse.support.data import core






from amuse.test import amusetest
import os.path
import numpy
from amuse.io import nemotsf
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
class Test(amusetest.TestCase):
    
    def setUp(self):
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory, 'p10.txt'), "r") as f:
            self.p10_string = f.read()
            
    def test1(self):
        instance = nemotsf.Tsf2Particles()
        
        particles = instance.convert_to_particles(self.p10_string)
                               
        self.assertEquals(instance.number_of_particles, 10)
        self.assertEquals(particles.mass[0], 0.1|nbody_system.mass)

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1|units.g, 1|units.m)
        instance = nemotsf.Tsf2Particles()
        particles = instance.convert_to_particles(self.p10_string, convert_nbody)
        self.assertAlmostEqual(particles.mass, units.g.new_quantity(0.1*numpy.ones(10)), constants.precision)

    def test3(self):
        instance = nemotsf.Tsf2Particles()
        particles = instance.convert_to_particles(self.p10_string)
        self.assertAlmostEqual(particles.mass[0], 0.1 | nbody_system.mass,  constants.precision)

    def test4(self):
        reader =  nemotsf.Tsf2Particles()
        particles = reader.convert_to_particles(self.p10_string)
        writer = nemotsf.Particles2Tsf()
        string = writer.convert_to_string(particles)
        print string
        print self.p10_string
        self.assertTrue("double PhaseSpace[10][2][3]" in string)
        self.assertTrue("double Mass[10]" in string)
        
