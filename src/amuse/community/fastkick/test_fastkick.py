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
    mode = "cpu"
    
    def test1(self):
        instance = FastKickInterface(mode=self.mode, number_of_workers=self.number_of_workers)
        instance.initialize_code()
        id1, error1 = instance.new_particle(mass = 11.0, x = 0.0, y = 0.0, z = 0.0)
        id2, error2 = instance.new_particle(mass = 21.0, x = 10.0, y = 0.0, z = 0.0)
        self.assertEqual(0, id1)
        self.assertEqual(0, error1)
        self.assertEqual(1, id2)
        self.assertEqual(0, error2)
        self.assertEqual(0, instance.commit_particles())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def test2(self):
        fastkick = FastKickInterface(mode=self.mode, number_of_workers=self.number_of_workers)
        self.assertEqual(0, fastkick.set_eps2(0.101))
        self.assertEqual([0.101, 0], list(fastkick.get_eps2().values()))
        self.assertEqual(0, fastkick.set_eps2(0.2))
        self.assertEqual([0.2, 0], list(fastkick.get_eps2().values()))
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test3(self):
        fastkick = FastKickInterface(mode=self.mode, number_of_workers=self.number_of_workers)
        fastkick.initialize_code()
        fastkick.new_particle([10,10],[-1,1],[0,0], [0,0])
        self.assertEqual(0, fastkick.commit_particles())
        self.assertEqual([-20.0, 0], list(fastkick.get_potential_at_point(0, 0,0,0).values()))
        self.assertAlmostEqual(-10.0*math.sqrt(2.0), list(fastkick.get_potential_at_point(1.0, 0,0,0).values())[0], 4)
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test4(self):
        numpy.random.seed(12345)
        plummer = new_plummer_model(100)
        fastkick = FastKickInterface(mode=self.mode, number_of_workers=self.number_of_workers)
        fastkick.initialize_code()
        fastkick.new_particle([1]*100, plummer.x.number, plummer.y.number, plummer.z.number)
        self.assertEqual(0, fastkick.commit_particles())
        points = new_plummer_model(73)
        potential1, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax1, ay1, az1, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        fastkick.cleanup_code()
        
        fastkick.initialize_code()
        self.assertEqual(0, fastkick.commit_particles())
        potential0, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax0, ay0, az0, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        self.assertAlmostEqual(potential0, 0)
        self.assertAlmostEqual(ax0, 0)
        self.assertAlmostEqual(ay0, 0)
        self.assertAlmostEqual(az0, 0)
        fastkick.cleanup_code()
        
        fastkick.initialize_code()
        for p in plummer:
            fastkick.new_particle(1, p.x.number, p.y.number, p.z.number)
        self.assertEqual(0, fastkick.commit_particles())
        potential2, error = fastkick.get_potential_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        ax2, ay2, az2, error = fastkick.get_gravity_at_point([0]*73, 
            points.x.number, points.y.number, points.z.number)
        self.assertAlmostEqual(potential1, potential2, 4)
        self.assertAlmostEqual(ax1, ax2, 4)
        self.assertAlmostEqual(ay1, ay2, 4)
        self.assertAlmostEqual(az1, az2, 4)
        fastkick.cleanup_code()
        fastkick.stop()
    


class TestFastKick(TestWithMPI):
    
    mode = "cpu"
    def new_fastkick_instance(self, convert_nbody=None, number_of_workers=2, **kwargs):
        return FastKick(convert_nbody, mode=self.mode, number_of_workers=number_of_workers, **kwargs)
    
    def test1(self):
        print("Testing FastKick (SI)")
        sun = Particle()
        sun.mass = 1.0|units.MSun
        sun.position = [0, 0, 0] | units.m
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        fastkick = self.new_fastkick_instance(convert_nbody)
        fastkick.particles.add_particle(sun)
        ax, ay, az = fastkick.get_gravity_at_point(0|units.AU, 1|units.AU, 0|units.AU, 0|units.AU)
        fastkick.stop()
        self.assertAlmostRelativeEqual(ax, -4*math.pi**2 | units.AU / units.yr**2, 3)
        self.assertAlmostRelativeEqual(ay, 0 | units.m * units.s**-2, 6)
        self.assertAlmostRelativeEqual(az, 0 | units.m * units.s**-2, 6)
        
    def test2(self):
        print("Testing FastKick reset")
        particles = new_plummer_model(50)
        instance = self.new_fastkick_instance()
        instance.parameters.epsilon_squared = 0.12345 | nbody_system.length**2
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)
        
        instance.reset()
        
        self.assertAlmostRelativeEquals(instance.parameters.epsilon_squared , 0.12345 | nbody_system.length**2)
        self.assertEqual(len(instance.particles), 0)
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)
    
    def test3(self):
        print("Testing FastKick parameters")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = self.new_fastkick_instance(convert_nbody)
        
        self.assertAlmostEqual(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEqual(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        
        instance.stop()
    
    def test4(self):
        numpy.random.seed(12345)
        plummer = new_plummer_model(100)
        fastkick = self.new_fastkick_instance()
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
        self.assertAlmostEqual(potential1, potential2, 5)
        self.assertAlmostEqual(ax1, ax2, 5)
        self.assertAlmostEqual(ay1, ay2, 5)
        self.assertAlmostEqual(az1, az2, 5)
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test5(self):
        print("Testing FastKick states")
        plummer = new_plummer_model(100)
        
        print("First do everything manually:")
        instance = self.new_fastkick_instance()
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEqual(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print("commit_parameters(), (re)commit_particles(), and cleanup_code() should be called " \
            "automatically before new_xx_particle(), get_xx(), and stop():")
        instance = self.new_fastkick_instance()
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        eps = instance.parameters.epsilon_squared
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(plummer[:-1])
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        acc = instance.get_gravity_at_point(0*plummer[-1].x, plummer[-1].x, plummer[-1].y, plummer[-1].z)
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.particles.add_particle(plummer[-1])
        self.assertEqual(instance.get_name_of_current_state(), 'UPDATE')
        acc = instance.get_gravity_at_point(0*plummer[-1].x, plummer[-1].x, plummer[-1].y, plummer[-1].z)
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
    
    def test6(self):
        plummer = new_plummer_model(100)
        points = new_plummer_model(73)
        instance = self.new_fastkick_instance(number_of_workers=1)
        instance.particles.add_particles(plummer)
        potential1 = instance.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax1, ay1, az1 = instance.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        instance.stop()
        
        expected_accuracy = 13 if self.mode == "cpu" else 5
        number_of_workers_range = [2, 3, 4, 5] if self.mode == "cpu" else [2]
        for n in number_of_workers_range:
            instance = self.new_fastkick_instance(number_of_workers=n)
            instance.particles.add_particles(plummer)
            potential = instance.get_potential_at_point(0*points.x, points.x, points.y, points.z)
            ax, ay, az = instance.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
            instance.stop()
            
            self.assertAlmostEqual(potential, potential1, expected_accuracy)
            self.assertAlmostEqual(ax, ax1, expected_accuracy)
            self.assertAlmostEqual(ay, ay1, expected_accuracy)
            self.assertAlmostEqual(az, az1, expected_accuracy)
    
    def test7(self):
        instance = self.new_fastkick_instance()
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
        print("Testing FastKick for Bridge: potential energy")
        numpy.random.seed(12345)
        stars = new_plummer_model(100)
        cluster = new_gravity_code(stars)
        
        first_half = stars.select_array(lambda x: (x > 0 | nbody_system.length), ['x'] )
        second_half = stars - first_half
        cluster1 = new_gravity_code(first_half)
        cluster2 = new_gravity_code(second_half)
        kick_from_cluster1 = CalculateFieldForCodesUsingReinitialize(self.new_fastkick_instance(), (cluster1,))
        kick_from_cluster2 = CalculateFieldForCodesUsingReinitialize(self.new_fastkick_instance(), (cluster2,))
        bridgesys = Bridge()
        bridgesys.add_system(cluster1, (kick_from_cluster2,))
        bridgesys.add_system(cluster2, (kick_from_cluster1,))
        
        self.assertAlmostRelativeEqual(cluster.potential_energy, bridgesys.potential_energy, 7)
        self.assertAlmostRelativeEqual(cluster.kinetic_energy, bridgesys.kinetic_energy)
        kick_from_cluster1.code.stop()
        kick_from_cluster2.code.stop()
    
    def test9(self):
        print("Testing FastKick for Bridge: evolving a binary")
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
        kick_from_primary = CalculateFieldForCodesUsingReinitialize(self.new_fastkick_instance(converter), (primary_sys,))
        kick_from_secondary = CalculateFieldForCodesUsingReinitialize(self.new_fastkick_instance(converter), (secondary_sys,))
        
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
    
    def test10(self):
        number_of_sources = 10
        number_of_points = 1000000
        mass, length, G = nbody_system.mass, nbody_system.length, nbody_system.G
        sources = Particles(mass=numpy.ones(number_of_sources)|mass, x=0|length, y=0|length, z=0|length)
        points = Particles(x=0|length, y=0|length, z=numpy.arange(number_of_points)+1.0|length)
        
        instance = self.new_fastkick_instance()
        instance.particles.add_particles(sources)
        potential = instance.get_potential_at_point(0*points.x, points.x, points.y, points.z)
        ax, ay, az = instance.get_gravity_at_point(0*points.x, points.x, points.y, points.z)
        self.assertAlmostEqual(ax, G * (0 | mass/length**2), 5)
        self.assertAlmostEqual(ay, G * (0 | mass/length**2), 5)
        self.assertAlmostRelativeEqual(az, -G*sources.total_mass()/points.z**2, 3)
        self.assertAlmostRelativeEqual(potential, -G*sources.total_mass()/points.z, 3)
        instance.stop()
    
    def test11(self):
        print("Test that a source is not included when calculating gravity on itself.")
        number_of_sources = 100
        mass, length, G = nbody_system.mass, nbody_system.length, nbody_system.G
        sources = Particles(mass=numpy.ones(number_of_sources)|mass, x=0|length, y=0|length, z=0|length)
        point = Particle(x=0|length, y=0|length, z=1.0|length)
        
        instance = self.new_fastkick_instance()
        instance.particles.add_particles(sources)
        potential = instance.get_potential_at_point(0|length, point.x, point.y, point.z)
        ax, ay, az = instance.get_gravity_at_point(0|length, point.x, point.y, point.z)
        self.assertAlmostEqual(ax, G * (0 | mass/length**2), 5)
        self.assertAlmostEqual(ay, G * (0 | mass/length**2), 5)
        self.assertAlmostRelativeEqual(az, -G*sources.total_mass()/point.z**2, 3)
        self.assertAlmostRelativeEqual(potential, -G*sources.total_mass()/point.z, 3)
        
        point.mass = 1e6 | mass
        instance.particles.add_particle(point)
        potential = instance.get_potential_at_point(0|length, point.x, point.y, point.z)
        ax, ay, az = instance.get_gravity_at_point(0|length, point.x, point.y, point.z)
        self.assertAlmostEqual(ax, G * (0 | mass/length**2), 5)
        self.assertAlmostEqual(ay, G * (0 | mass/length**2), 5)
        self.assertAlmostRelativeEqual(az, -G*sources.total_mass()/point.z**2, 3)
        self.assertAlmostRelativeEqual(potential, -G*sources.total_mass()/point.z, 3)
        instance.stop()
    
    def test12(self):
        print("Test FastKick potential energy")
        number_of_sources = 1000
        numpy.random.seed(12345)
        plummer = new_plummer_model(number_of_sources)
        
        fastkick = self.new_fastkick_instance()
        fastkick.particles.add_particles(plummer)
        self.assertAlmostRelativeEqual(fastkick.potential_energy, plummer.potential_energy(G=nbody_system.G), 6)
        fastkick.cleanup_code()
        fastkick.stop()
    
    def test13(self):
        number_of_sources = 10
        mass, length, G = nbody_system.mass, nbody_system.length, nbody_system.G
        sources = Particles(
                mass=numpy.ones(number_of_sources)|mass,
                x=0|length, 
                y=0|length, 
                z=0|length)
        
        instance = self.new_fastkick_instance()
        instance.particles.add_particles(sources)
        self.assertAlmostRelativeEqual(instance.particles.mass, 1 | mass)
        instance.particles.mass = 2 | mass
        self.assertAlmostRelativeEqual(instance.particles.mass, 2 | mass)

class TestFastKickGPU(TestFastKick):

    mode = "gpu"
    def new_fastkick_instance(self, convert_nbody=None, number_of_workers=1, **kwargs):
        return self.new_instance_of_an_optional_code(FastKick, convert_nbody, 
            mode=self.mode, number_of_workers=number_of_workers, **kwargs)

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

