import numpy

from amuse.test import amusetest
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.plummer import new_plummer_model, MakePlummerModel

class TestData(amusetest.TestCase):
    pass
   
class TestPlummer(TestData):
    def test1(self):
        numpy.random.seed(0)
        #print numpy.random.get_state()
        m = MakePlummerModel(2)
        m1, p, v = m.new_model()
        self.assertEqual(m1[0,0], 0.5)
        self.assertEqual(m1[1,0], 0.5)
        self.assertAlmostEqual(p[0,0], -0.729636617171, 5)
        self.assertAlmostEqual(p[1,0], -0.713272921751 , 5)
        self.assertAlmostEqual(p[0,1],  0.379570256435, 5)
        self.assertAlmostEqual(p[1,1],  -0.930290757081, 5)
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(6|units.kg, 7 | units.m) 
        stars =  new_plummer_model(2, convert_nbody)
        self.assertEqual(stars[0].mass.value_in(units.kg), 3.0)
        self.assertEqual(stars[1].mass.value_in(units.kg), 3.0)
        
    def test3(self):
        stars =  new_plummer_model(2, None)
        self.assertEqual(stars[0].mass.value_in(nbody_system.mass), 0.5)
        self.assertEqual(stars[1].mass.value_in(nbody_system.mass), 0.5)
        
    def test4(self):
        stars = new_plummer_model(2, do_scale = True)
        self.assertAlmostEqual(stars.kinetic_energy(),             0.25 | nbody_system.energy)
        self.assertAlmostEqual(stars.potential_energy(G=nbody_system.G), -0.50 | nbody_system.energy)
        self.assertAlmostEqual(stars.center_of_mass(),          [0,0,0] | nbody_system.length)
        self.assertAlmostEqual(stars.center_of_mass_velocity(), [0,0,0] | nbody_system.speed)
        self.assertAlmostEqual(stars.mass.sum(),                   1.00 | nbody_system.mass)
        self.assertAlmostEqual(stars.virial_radius(),              1.00 | nbody_system.length)
        
