import numpy.random

from amuse.test import amusetest
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system

from amuse.ic.plummer import new_plummer_sphere
from amuse.datamodel import Particle, Particles
from amuse.datamodel import particle_attributes

class TestParticlesAttributes(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(2)
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
        particles = Particles(2)
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
        
    def test3(self):
        numpy.random.seed(123)
        plummer = new_plummer_sphere(10000)
        result = plummer.new_particle_from_cluster_core()
        self.assertTrue(isinstance(result, Particle))
        self.assertAlmostEqual(result.density, 3.55015420914 | nbody_system.mass * nbody_system.length**-3)
        self.assertAlmostEqual(result.radius, 0.136202583288 | nbody_system.length)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | nbody_system.speed, 1)
        
        plummer.vx += 42 | nbody_system.speed
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|nbody_system.speed)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | nbody_system.length
        result = plummer.new_particle_from_cluster_core()
        self.assertTrue(isinstance(result, Particle))
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e3 | nbody_system.mass * nbody_system.length**-3, 4)
        self.assertAlmostEqual(result.radius, 0.0136203656761 | nbody_system.length)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | nbody_system.speed, 1)
    
