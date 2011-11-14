import numpy

from amuse.test import amusetest
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.gasplummer import new_plummer_gas_model, MakePlummerGasModel

class TestPlummerGasModel(amusetest.TestCase):
    
    def test1(self):
        print "Test 1: testing low-level interface (no units or datamodel)"
        numpy.random.seed(345672)
        mpgm = MakePlummerGasModel(2)
        mass, x,y,z, vx,vy,vz, u = mpgm.new_model()
        self.assertEquals(mass[0], 0.5)
        self.assertEquals(mass[1], 0.5)
        self.assertAlmostEqual(x, [-0.02295788, 0.12829775])
        self.assertAlmostEqual(y, [-0.41054985, 0.14190860])
        self.assertAlmostEqual(z, [-0.50723639, 0.08937734])
        self.assertAlmostEqual(vx, [0.0, 0.0])
        self.assertAlmostEqual(vy, [0.0, 0.0])
        self.assertAlmostEqual(vz, [0.0, 0.0])
        self.assertAlmostEqual(u, [0.28413716, 0.39898137])
    
    def test2(self):
        print "Test 2: testing user interface, with convert_nbody -> SI units"
        convert_nbody = nbody_system.nbody_to_si(6|units.kg, 7 | units.m) 
        gas =  new_plummer_gas_model(2, convert_nbody)
        self.assertEquals(gas[0].mass.value_in(units.kg), 3.0)
        self.assertEquals(gas[1].mass.value_in(units.kg), 3.0)
    
    def test3(self):
        print "Test 3: testing user interface, without convert_nbody -> nbody units"
        gas =  new_plummer_gas_model(2, None)
        self.assertEquals(gas[0].mass.value_in(nbody_system.mass), 0.5)
        self.assertEquals(gas[1].mass.value_in(nbody_system.mass), 0.5)
    
    def test4(self):
        print "Test 4: test new_plummer_gas_model, model properties"
        numpy.random.seed(345672)
        gas = new_plummer_gas_model(100)
        
        self.assertEqual(len(gas), 100)
        self.assertAlmostEqual(gas.kinetic_energy(), 0.00 | nbody_system.energy)
        self.assertIsOfOrder(  gas.thermal_energy(), 0.25 | nbody_system.energy)
        self.assertAlmostEqual(gas.thermal_energy(), 0.238075609078 | nbody_system.energy)
        self.assertIsOfOrder(  gas.potential_energy(G=nbody_system.G), -0.50 | nbody_system.energy)
        self.assertAlmostEqual(gas.potential_energy(G=nbody_system.G), -0.447052244411 | nbody_system.energy)
        
        self.assertAlmostEqual(gas.center_of_mass(),          [0,0,0] | nbody_system.length)
        self.assertAlmostEqual(gas.center_of_mass_velocity(), [0,0,0] | nbody_system.speed)
        self.assertAlmostEqual(gas.total_mass(),                 1.00 | nbody_system.mass)
        self.assertIsOfOrder(gas.virial_radius(),                1.00 | nbody_system.length)
        self.assertAlmostEqual(gas.virial_radius(),     1.11843751206 | nbody_system.length)
    
    def test5(self):
        print "Test 5: test new_plummer_gas_model with do_scale"
        gas = new_plummer_gas_model(100, do_scale = True)
        
        self.assertEqual(len(gas), 100)
        self.assertAlmostEqual(gas.kinetic_energy(), 0.00 | nbody_system.energy)
        self.assertAlmostEqual(gas.thermal_energy(), 0.25 | nbody_system.energy)
        self.assertAlmostEqual(gas.potential_energy(G=nbody_system.G), -0.50 | nbody_system.energy)
        self.assertAlmostEqual(gas.center_of_mass(),          [0,0,0] | nbody_system.length)
        self.assertAlmostEqual(gas.center_of_mass_velocity(), [0,0,0] | nbody_system.speed)
        self.assertAlmostEqual(gas.total_mass(),                 1.00 | nbody_system.mass)
        self.assertAlmostEqual(gas.virial_radius(),              1.00 | nbody_system.length)
    
