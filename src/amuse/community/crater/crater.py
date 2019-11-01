from amuse.test.amusetest import TestWithMPI
from amuse.lab import units, constants

from interface import CraterInterface
from interface import Crater

instance = CraterInterface(redirect="none")

instance.set_target_density(3.0)
instance.set_target_gravity(9.8)
instance.set_projectile_density(2.6)
instance.set_projectile_diameter(100)
instance.set_impact_angle(90)
instance.set_impact_velocity(30)

#instance.evolve_model(1.0)
Dt = 333 # target diameter in km
rhoproj = 3.0 # density of the projectile
L = 1 # projectile size in km
v = 10 # units.kms
theta = 90 # impact angle in degrees (head on=180)
cratertype = 0 # output: 0:simple, 1: complex, 2: simple/complex
rhotarget = 3 # density of the target object 
g = 1 # no idea
targtype = 3 # target type 1:, 2:, 3:
Dyield = 0 # output
Dgault = 0# output
Tform = 0# output
Dfinal = 0 # final crater diameter in km
Lyield = 0# output
Lgault = 0# output

instance.crater_radius(rhoproj, L, v, theta,
                       rhotarget, g, targtype)

print "D=", instance.get_crater_diameter()
print "c=", instance.get_crater_type()

instance.stop()

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

