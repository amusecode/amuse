import numpy
import math

from amuse.units import nbody_system, units, constants, quantities
from amuse.datamodel import Particles, Particle
from amuse.ic.plummer import new_plummer_model
from amuse.couple.bridge import Bridge, CalculateFieldForCodesUsingReinitialize

from amuse.test.amusetest import TestWithMPI
from amuse.community.fastkick.interface import FastKickInterface, FastKick


class TestFastKickInterface(TestWithMPI):

    number_of_workers = 2
    
    def test1(self):
        instance = FastKickInterface(number_of_workers = self.number_of_workers)
        instance.initialize_code()
        id1, error1 = instance.new_particle(mass = 11.0, x = 0.0, y = 0.0, z = 0.0)
        id2, error2 = instance.new_particle(mass = 21.0, x = 10.0, y = 0.0, z = 0.0)
        self.assertEquals(0, id1)
        self.assertEquals(0, error1)
        self.assertEquals(1, id2)
        self.assertEquals(0, error2)
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        fastkick = FastKickInterface(number_of_workers = self.number_of_workers)
        self.assertEquals(0, fastkick.set_eps2(0.101))
        self.assertEquals([0.101, 0], fastkick.get_eps2().values())
        self.assertEquals(0, fastkick.set_eps2(0.2))
        self.assertEquals([0.2, 0], fastkick.get_eps2().values())
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test3(self):
        fastkick = FastKickInterface(number_of_workers = self.number_of_workers)
        fastkick.initialize_code()
        fastkick.new_particle([10,10],[-1,1],[0,0], [0,0])
        self.assertEqual([-20.0, 0], fastkick.get_potential_at_point(0, 0,0,0).values())
        self.assertAlmostEqual(-10.0*math.sqrt(2.0), fastkick.get_potential_at_point(1.0, 0,0,0).values()[0])
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test4(self):
        numpy.random.seed(12345)
        plummer = new_plummer_model(100)
        fastkick = FastKickInterface(number_of_workers = self.number_of_workers)
        fastkick.initialize_code()
        fastkick.new_particle([1]*100, plummer.x.number, plummer.y.number, plummer.z.number)
        points = new_plummer_model(73)
        potential1, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax1, ay1, az1, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        fastkick.cleanup_code()
        
        fastkick.initialize_code()
        potential0, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax0, ay0, az0, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        self.assertAlmostEqual(potential0, 0)
        self.assertAlmostEqual(ax0, 0)
        self.assertAlmostEqual(ay0, 0)
        self.assertAlmostEqual(az0, 0)
        
        for p in plummer:
            fastkick.new_particle(1, p.x.number, p.y.number, p.z.number)
        potential2, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax2, ay2, az2, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        self.assertAlmostEqual(potential1, potential2)
        self.assertAlmostEqual(ax1, ax2)
        self.assertAlmostEqual(ay1, ay2)
        self.assertAlmostEqual(az1, az2)
        fastkick.cleanup_code()
        fastkick.stop()
    


class TestFastKick(TestWithMPI):
    
    number_of_workers = 2
    
    def test1(self):
        print "Testing FastKick (SI)"
        sun = Particle()
        sun.mass = 1.0|units.MSun
        sun.position = [0, 0, 0] | units.m
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        fastkick = FastKick(convert_nbody, number_of_workers = self.number_of_workers)
        fastkick.particles.add_particle(sun)
        ax, ay, az = fastkick.get_gravity_at_point(0|units.AU, 1|units.AU, 0|units.AU, 0|units.AU)
        fastkick.stop()
        self.assertAlmostRelativeEqual(ax, -4*math.pi**2 | units.AU / units.yr**2, 3)
        self.assertAlmostRelativeEqual(ay, 0 | units.m * units.s**-2, 6)
        self.assertAlmostRelativeEqual(az, 0 | units.m * units.s**-2, 6)
        
    def test2(self):
        print "Testing FastKick reset"
        particles = new_plummer_model(50)
        instance = FastKick(number_of_workers = self.number_of_workers)
        instance.parameters.epsilon_squared = 0.12345 | nbody_system.length**2
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 50)
        
        instance.reset()
        
        self.assertAlmostRelativeEquals(instance.parameters.epsilon_squared , 0.12345 | nbody_system.length**2)
        self.assertEquals(len(instance.particles), 0)
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 50)
    
    def test3(self):
        print "Testing FastKick parameters"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = FastKick(convert_nbody, number_of_workers = self.number_of_workers)
        
        self.assertAlmostEquals(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEquals(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        
        instance.stop()
    
    def test4(self):
        numpy.random.seed(12345)
        plummer = new_plummer_model(100)
        fastkick = FastKick(number_of_workers = self.number_of_workers)
        fastkick.initialize_code()
        fastkick.particles.add_particles(plummer)
        points = new_plummer_model(73)
        potential1 = fastkick.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax1, ay1, az1 = fastkick.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        fastkick.reset()
        
        potential0 = fastkick.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax0, ay0, az0 = fastkick.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        self.assertAlmostEqual(potential0, 0 | nbody_system.potential)
        self.assertAlmostEqual(ax0, 0 | nbody_system.acceleration)
        self.assertAlmostEqual(ay0, 0 | nbody_system.acceleration)
        self.assertAlmostEqual(az0, 0 | nbody_system.acceleration)
        fastkick.reset()
        
        for p in plummer:
            fastkick.particles.add_particle(p)
        potential2 = fastkick.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax2, ay2, az2 = fastkick.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        self.assertAlmostEqual(potential1, potential2)
        self.assertAlmostEqual(ax1, ax2)
        self.assertAlmostEqual(ay1, ay2)
        self.assertAlmostEqual(az1, az2)
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test5(self):
        print "Testing FastKick states"
        plummer = new_plummer_model(100)
        
        print "First do everything manually:"
        instance = FastKick(number_of_workers = self.number_of_workers)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print "commit_parameters(), (re)commit_particles(), and cleanup_code() should be called " \
            "automatically before new_xx_particle(), get_xx(), and stop():"
        instance = FastKick(number_of_workers = self.number_of_workers)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        eps = instance.parameters.epsilon_squared
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(plummer[:-1])
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        acc = instance.get_gravity_at_point(0*plummer[-1].x, plummer[-1].x, plummer[-1].y, plummer[-1].z)
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.particles.add_particle(plummer[-1])
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        acc = instance.get_gravity_at_point(0*plummer[-1].x, plummer[-1].x, plummer[-1].y, plummer[-1].z)
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
    
    def test6(self):
        plummer = new_plummer_model(100)
        points = new_plummer_model(73)
        instance = FastKick(number_of_workers=1)
        instance.particles.add_particles(plummer)
        potential1 = instance.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax1, ay1, az1 = instance.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        instance.stop()
        
        for n in [2, 3, 4, 5]:
            instance = FastKick(number_of_workers=n)
            instance.particles.add_particles(plummer)
            potential = instance.get_potential_at_point(0*points.x, points.x, points.y, points.z)
            ax, ay, az = instance.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
            instance.stop()
            
            self.assertAlmostEqual(potential, potential1, 13)
            self.assertAlmostEqual(ax, ax1, 13)
            self.assertAlmostEqual(ay, ay1, 13)
            self.assertAlmostEqual(az, az1, 13)
    
    def test7(self):
        instance = FastKick(number_of_workers = self.number_of_workers)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | nbody_system.length**2
        
        particles = Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        instance.particles.add_particles(particles)
        
        zero = 0.0 | nbody_system.length
        ax, ay, az = instance.get_gravity_at_point(zero, 1.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(ax, 0.0 | nbody_system.acceleration, 6)
        self.assertAlmostEqual(ay, 0.0 | nbody_system.acceleration, 6)
        self.assertAlmostEqual(az, 0.0 | nbody_system.acceleration, 6)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            ax0, ay0, az0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            ax1, ay1, az1 = instance.get_gravity_at_point(zero, x1, zero, zero)
            
            self.assertAlmostEqual(ay0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(az0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(ay1, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(az1, 0.0 | nbody_system.acceleration, 6)
            
            self.assertAlmostEqual(ax0, -ax1, 5)
            ax = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length**3 / nbody_system.time**2)
            self.assertAlmostEqual(ax, ax0, 2)
            self.assertAlmostEqual(potential0, potential1, 5)
        instance.stop()
    
    def test8(self):
        print "Testing FastKick for Bridge: potential energy"
        numpy.random.seed(12345)
        stars = new_plummer_model(100)
        cluster = new_gravity_code(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | nbody_system.length), ['x'] )
        second_half = stars - first_half
        cluster1 = new_gravity_code(first_half)
        cluster2 = new_gravity_code(second_half)
        kick_from_cluster1 = CalculateFieldForCodesUsingReinitialize(FastKick(), (cluster1,))
        kick_from_cluster2 = CalculateFieldForCodesUsingReinitialize(FastKick(), (cluster2,))
        bridgesys = Bridge()
        bridgesys.add_system(cluster1, (kick_from_cluster2,))
        bridgesys.add_system(cluster2, (kick_from_cluster1,))
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
        kick_from_cluster1.code.stop()
        kick_from_cluster2.code.stop()
    
    def test9(self):
        print "Testing FastKick for Bridge: evolving a binary"
        particles = Particles(2)
        particles.mass = [3.0, 1.0] | units.MSun
        particles.position = [0, 0, 0] | units.AU
        particles.velocity = [0, 0, 0] | units.km / units.s
        particles[1].x = 2.0 | units.AU
        particles[1].vy = (constants.G * (4.0 | units.MSun) / (2.0 | units.AU)).sqrt()
        particles.move_to_center()
        
        primary_sys = new_gravity_code(particles[:1])
        secondary_sys = new_gravity_code(particles[1:])
        
        primary = primary_sys.particles[0]
        P = 2 * math.pi * primary.x / primary.vy
        
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        kick_from_primary = CalculateFieldForCodesUsingReinitialize(FastKick(converter), (primary_sys,))
        kick_from_secondary = CalculateFieldForCodesUsingReinitialize(FastKick(converter), (secondary_sys,))
        
        bridgesys = Bridge(timestep = P / 64.0)
        bridgesys.add_system(primary_sys, (kick_from_secondary,))
        bridgesys.add_system(secondary_sys, (kick_from_primary,))
        
        position_at_start = primary.position.x
        bridgesys.evolve_model(P / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 2)
        
        bridgesys.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 2)
        
        bridgesys.evolve_model(P)
        kick_from_primary.code.stop()
        kick_from_secondary.code.stop()
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 2)
    

class GravityCodeForTesting(object):
    
    def __init__(self):
        self.particles = Particles()
        self.model_time = quantities.zero
    
    def evolve_model(self, t_end):
        self.particles.position += self.particles.velocity * (t_end - self.model_time)
        self.model_time = t_end
    
    @property
    def potential_energy(self):
        G = nbody_system.G if self.particles.x.unit == nbody_system.length else constants.G
        return self.particles.potential_energy(G=G)
    
    @property
    def kinetic_energy(self):
        return self.particles.kinetic_energy()
    

def new_gravity_code(particles):
    instance = GravityCodeForTesting()
    instance.particles.add_particles(particles)
    return instance