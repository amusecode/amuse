from amuse.test.amusetest import TestWithMPI
from amuse.lab import units, constants
from interface import Crater

instance = Crater(redirection="none")
instance.initialize_code()
print "crater initialized."
instance.set_target_density(3.0|units.g/units.cm**3)
instance.set_target_gravity(9.8|units.m/units.s**2)
instance.set_target_type(3)
instance.set_projectile_diameter(10.0|units.km)
instance.set_projectile_density(3.0|units.g/units.cm**3)
instance.set_projectile_type(3)
instance.set_impact_angle(90|units.deg)
instance.set_impact_velocity(30|units.kms)

instance.evolve_model()

print "crater D=", instance.get_crater_diameter().in_(units.km)
print "crater t=", instance.get_formation_time().in_(units.s)
print "crater type=", instance.get_crater_type()

instance.stop()

