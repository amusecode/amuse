from amuse.test.amusetest import TestWithMPI
from amuse.lab import units

from interface import CraterInterface
from interface import Crater

class CraterInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = Crater()
        instance.set_target_density(3.0|units.g/units.cm**3)
        rho = instance.get_target_density()
        self.assertAlmostRelativeEqual(rho, 3.0 | units.g/units.cm**3, 8)

    def test2(self):
        instance = Crater(redirection="none")
        instance.initialize_code()
        instance.set_projectile_diameter(10.0|units.km)
        d = instance.get_projectile_diameter()
        self.assertAlmostRelativeEqual(d, 10000.0 | units.m, 8)
        instance.stop()

    def test3(self):
        instance = Crater(redirection="none")
        instance.initialize_code()
        instance.set_target_density(3.34|units.g/units.cm**3)
        instance.set_target_gravity(1.62|units.m/units.s**2)
        instance.set_target_type(2)
        instance.set_projectile_diameter(10.0|units.km)
        instance.set_projectile_density(2.6|units.g/units.cm**3)
        instance.set_projectile_type(3)
        instance.set_impact_angle(45|units.deg)
        instance.set_impact_velocity(30|units.kms)
        instance.evolve_model()
        d = instance.get_crater_diameter()
        self.assertAlmostRelativeEqual(d, 87.3167847257 | units.km, 8)
        instance.stop()

    def test4(self):
        instance = Crater(redirection="none")
        instance.initialize_code()
        instance.set_target_density(3.0|units.g/units.cm**3)
        instance.set_target_gravity(9.8|units.m/units.s**2)
        instance.set_target_type(3)
        instance.set_projectile_diameter(10.0|units.km)
        instance.set_projectile_density(3.0|units.g/units.cm**3)
        instance.set_projectile_type(3)
        instance.set_impact_angle(90|units.deg)
        instance.set_impact_velocity(30|units.kms)
        instance.evolve_model()
        d = instance.get_crater_diameter()
        self.assertAlmostRelativeEqual(d, 174.09185123 | units.km, 8)
        instance.stop()
        
        
        
