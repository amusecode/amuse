import numpy.random

from amuse.test import amusetest
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support.interface import ConvertArgumentsException

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
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, stop_hop=False)
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
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, stop_hop=True)
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
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, stop_hop=False)
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
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, stop_hop=True)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e6 | units.MSun * units.parsec**-3, 4)
    
    def test5(self):
        print "Test new_particle_from_cluster_core - stop_hop or not"
        converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(100, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, stop_hop=False)
        
        # Hop wasn't stopped, will use same Hop instance:
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, stop_hop=False)
        
        nbody_plummer = new_plummer_sphere(100)
        # Hop wasn't stopped, unit_converters don't match:
        self.assertRaises(AttributeError, nbody_plummer.new_particle_from_cluster_core, 
            expected_message="Cannot combine units from different systems: m and length")
        
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, stop_hop=True)
        
        # Hop was stopped, new instance will be made with supplied unit_converter (None in this case):
        result = nbody_plummer.new_particle_from_cluster_core(stop_hop=False)
        
        self.assertRaises(ConvertArgumentsException, plummer.new_particle_from_cluster_core, unit_converter=converter,#,
            expected_message="error while converting parameter 'mass', error: Cannot express kg in mass, the units do not have the same bases")
        
        result = nbody_plummer.new_particle_from_cluster_core(stop_hop=True)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, stop_hop=True)
    


class TestParticlesDomainAttributes(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        
        self.assertAlmostRelativeEqual(particles.a.foo,  1 | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  2 | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  3 | units.s)
        
    def test2(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.bar = 2 | units.kg
        particles.foo = 3 | units.s
        self.assertEquals(
            sorted(particles.a.get_attribute_names_defined_in_store()), 
            ['foo']
        )
        self.assertEquals(
            sorted(particles.b.get_attribute_names_defined_in_store()), 
            ['bar']
        )
        
    def test3(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        particles.a.add_particle(Particle(foo = 2 | units.m))
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [1,1,2] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2,0] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3,0] | units.s)
        
        particles.add_particle(Particle(foo = 2 | units.s))
        self.assertAlmostRelativeEqual(particles.a.foo,  [1,1,2,0] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2,0,0] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3,0,2] | units.s)
        
        
    def test4(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        def set_a():
            particles.a = 1 | units.kg
        self.assertRaises(AttributeError, set_a)

    def test5(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        particles.a.foo = 1 | units.m
        particles.b.foo = 2 | units.kg
        particles.foo = 3 | units.s
        
        particles.a[0].foo = 3 | units.m
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [3,1] | units.m)
        
        particles[0].a.foo = 4 | units.m
        
        self.assertAlmostRelativeEqual(particles.a.foo,  [4,1] | units.m)
        self.assertAlmostRelativeEqual(particles.b.foo,  [2,2] | units.kg)
        self.assertAlmostRelativeEqual(particles.foo,  [3,3] | units.s)

    def test6(self):
        particles = Particles(2)
        particles.add_attribute_domain('a')
        particles.add_attribute_domain('b')
        def set_a():
            particles[0].a = 1 | units.kg
        self.assertRaises(AttributeError, set_a)
