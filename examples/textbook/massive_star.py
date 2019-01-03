# from amuse.community.mosse.interface import MOSSE
from amuse.community.seba.interface import SeBa

from amuse.units import units
from amuse.datamodel import Particle

star = Particle()

star.mass = 16 | units.MSun

# SE = MOSSE()
SE = SeBa()

SE.particles.add_particle(star)

print SE.particles[0]

SE.stop()
