from amuse.community.bhtree.interface import BHTree
from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
instance = BHTree(convert_nbody,channel_type='ibis')
instance.parameters.epsilon_squared = 0.001 | units.AU**2
stars = core.Particles(2)
sun = stars[0]
sun.mass = 1.0 | units.MSun
sun.position = [0.0,0.0,0.0] | units.m
sun.velocity = [0.0,0.0,0.0] | units.m / units.s
sun.radius = 1.0 | units.RSun
earth = stars[1]
earth.mass = 5.9736e24 | units.kg
earth.radius = 6371.0 | units.km
earth.position = [1.0, 0.0, 0.0] | units.AU
earth.velocity = [0.0, 29783, 0.0] | units.ms
instance.particles.add_particles(stars)
channel = instance.particles.new_channel_to(stars)
print earth.position[0]
print earth.position.as_quantity_in(units.AU)[0]
instance.evolve_model(1.0 | units.yr)
print earth.position.as_quantity_in(units.AU)[0] # This is the outdated value! (should update_particles first)
channel.copy()
print earth.position.as_quantity_in(units.AU)[0]
instance.evolve_model(1.5 | units.yr)
channel.copy()
print earth.position.as_quantity_in(units.AU)[0]
instance.stop()
