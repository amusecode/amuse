import unittest

from stellar.single.EFT89.muse_stellar import EFT89
from stellar.single.SSE.muse_stellar import SSE
#import sys
#import os

class TestEvolution(unittest.TestCase):
   pass
        

class TestEFT89(TestEvolution):
    def test1(self):
        mass = 30.0
        interface = EFT89()
        interface.setup_module()
        interface.add_zams_star(1,mass)
        self.assertAlmostEqual(interface.get_mass(1), mass)
        self.assertEqual(interface.get_age(1), 0.0)
        dt = interface.get_time_step(1,dm=10)
        expected_masses = [30.0,30.0,30.0,30.0,30.0,30.0,0.8,0.8,0.8,0.8]  
        expected_radii = [6.97475696739,7.97478562361,8.6932129409,9.8109171988,12.4481130671,19.282293503,0.01,0.01,0.01,0.01]
        for x in range(0,10):
            interface.evolve(1,x)
            self.assertEqual(interface.get_age(1), x)
            print interface.get_luminosity(1)
            self.assertAlmostEqual(interface.get_mass(1), expected_masses[x])
            self.assertAlmostEqual(interface.get_radius(1), expected_radii[x])
        self.assertEqual(interface.get_age(1), 9)
        
        
class TestSSE(TestEvolution):
    def test1(self):
        mass = 30.0
        interface = SSE()
        interface.setup_module()
        interface.add_zams_star(1,mass)
        self.assertAlmostEqual(interface.get_mass(1), mass,1)
        self.assertEqual(interface.get_age(1), 0.0)
        dt = interface.get_time_step(1,dm=10)
        array = []
        expected_masses = [29.99, 29.82, 29.55, 29.32, 28.88, 28.46, 26.52, 8.57, 8.57, 8.57] 
        expected_radii = [7.74, 8.51, 9.66, 10.78, 13.30,16.77, 1274.92, 3.634e-05, 3.634e-05, 3.634e-05]
        for x in range(0,10):
            interface.evolve(1,x)
            self.assertEqual(interface.get_age(1), x)
            array.append(interface.get_radius(1))
            self.assertAlmostEqual(interface.get_mass(1), expected_masses[x],1)
            self.assertAlmostEqual(interface.get_radius(1), expected_radii[x],1)
        self.assertEqual(interface.get_age(1), 9)
        