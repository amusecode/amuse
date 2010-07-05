from amuse.test import amusetest

from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.support.units import nbody_system as nbody
from amuse.support.units import units

class TestEvrardModel(amusetest.TestCase):
    def test1(self):
        print "First test: making an Evrard gas sphere model."
        target_number_of_particles = 1000
        gas_parts = new_evrard_gas_sphere(target_number_of_particles, seed=1234)
        self.assertEquals(len(gas_parts), 978)
    
    def test2(self):
        print "Testing properties of an Evrard model."
        target_number_of_particles = 1000
        gas_parts = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed=1234)
        self.assertEquals(len(gas_parts), 978)
        self.assertAlmostEqual(gas_parts.kinetic_energy(),             0.00 | nbody.energy)
        self.assertAlmostEqual(gas_parts.potential_energy(G=nbody.G), -0.50 | nbody.energy)
        self.assertAlmostEqual(gas_parts.center_of_mass(),          [0,0,0] | nbody.length)
        self.assertAlmostEqual(gas_parts.center_of_mass_velocity(), [0,0,0] | nbody.speed)
        self.assertAlmostEqual(gas_parts.mass.sum(),                   1.00 | nbody.mass)
    
    def test3(self):
        print "Testing virial radius of an Evrard model."
        target_number_of_particles = 100
        gas_parts = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed=1234)
        self.assertAlmostEqual(gas_parts.virial_radius(),             1.00 | nbody.length)
    
