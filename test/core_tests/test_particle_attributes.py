from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system

from amuse.support.data import core
from amuse.support.data import particle_attributes

from amuse.test import amusetest

class TestParticlesAttributes(amusetest.TestCase):
    
    def test1(self):
        particles = core.Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | nbody_system.length
        particles.velocity = [[-1, 0, 0], [1,0,0]] | nbody_system.length/nbody_system.time
        particles.mass = 0.4 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass.sum(), 0.8 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.4 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.08 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | nbody_system.length)
        particles.scale_to_standard()
        self.assertAlmostRelativeEquals(particles.mass.sum(), 1.0 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.25 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.5 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 1.0 | nbody_system.length)
        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        particles = core.Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | units.parsec
        particles.velocity = [[-1, 0, 0], [1,0,0]] | units.parsec / units.Myr
        particles.mass = 0.5 | units.MSun 
        
        self.assertAlmostRelativeEquals(particles.mass.sum(), 1.0 | units.MSun)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 1.0 * (0.5 | units.MSun) * (1 |units.parsec / units.Myr) **2 )
        self.assertAlmostRelativeEquals(particles.potential_energy(), -constants.G *  (0.5 | units.MSun) ** 2  / ([2,0,0] | units.parsec).length() )
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | units.parsec)
        
        particles.scale_to_standard(convert_nbody)
        self.assertAlmostRelativeEquals(particles.mass.sum(), convert_nbody.to_si(1.0 | nbody_system.mass))
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), convert_nbody.to_si(0.25 | nbody_system.energy))
        self.assertAlmostRelativeEquals(particles.potential_energy().as_quantity_in(units.J), convert_nbody.to_si(-0.5 | nbody_system.energy).as_quantity_in(units.J), 12)
        self.assertAlmostRelativeEquals(particles.virial_radius(), convert_nbody.to_si(1.0 | nbody_system.length))
