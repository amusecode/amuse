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
        print "Test basic particle attributes and scale_to_standard - nbody units"
        particles = Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | nbody_system.length
        particles.velocity = [[-1, 0, 0], [1,0,0]] | nbody_system.length/nbody_system.time
        particles.mass = 0.4 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.total_mass(), 0.8 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.4 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.08 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | nbody_system.length)
        particles.scale_to_standard()
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 0.25 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), -0.5 | nbody_system.energy)
        self.assertAlmostRelativeEquals(particles.virial_radius(), 1.0 | nbody_system.length)
        
    def test2(self):
        print "Test basic particle attributes and scale_to_standard - SI units"
        convert_nbody = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        particles = Particles(2)
        particles.position = [[-1, 0, 0], [1,0,0]] | units.parsec
        particles.velocity = [[-1, 0, 0], [1,0,0]] | units.parsec / units.Myr
        particles.mass = 0.5 | units.MSun 
        
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0 | units.MSun)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 1.0 * (0.5 | units.MSun) * (1 |units.parsec / units.Myr) **2 )
        self.assertAlmostRelativeEquals(particles.potential_energy(), -constants.G *  (0.5 | units.MSun) ** 2  / ([2,0,0] | units.parsec).length() )
        self.assertAlmostRelativeEquals(particles.virial_radius(), 4.0 | units.parsec)
        
        particles.scale_to_standard(convert_nbody)
        self.assertAlmostRelativeEquals(particles.total_mass(), convert_nbody.to_si(1.0 | nbody_system.mass))
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), convert_nbody.to_si(0.25 | nbody_system.energy))
        self.assertAlmostRelativeEquals(particles.potential_energy().as_quantity_in(units.J), convert_nbody.to_si(-0.5 | nbody_system.energy).as_quantity_in(units.J), 12)
        self.assertAlmostRelativeEquals(particles.virial_radius(), convert_nbody.to_si(1.0 | nbody_system.length))
        
    def test3(self):
        print "Test new_particle_from_cluster_core - nbody units"
        numpy.random.seed(123)
        plummer = new_plummer_sphere(10000)
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | nbody_system.length
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostEqual(result.density, 3.55015420914 | nbody_system.density)
        
        plummer.vx += 42 | nbody_system.speed
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|nbody_system.speed)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | nbody_system.length
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e3 | nbody_system.density, 4)
    
    def test4(self):
        print "Test new_particle_from_cluster_core - SI units"
        numpy.random.seed(123)
        converter = nbody_system.nbody_to_si(1000.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(10000, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | units.parsec
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostEqual(result.density, 3.55015420914e3 | units.MSun * units.parsec**-3)
        
        plummer.vx += 42 | units.km / units.s
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|units.km / units.s)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | units.parsec
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e6 | units.MSun * units.parsec**-3, 4)
    
