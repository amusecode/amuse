import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.units import nbody_system, units, constants
from amuse.datamodel import Particles
from amuse.ic.plummer import new_plummer_model

try:
    from amuse.community.tupan.interface import TupanInterface, Tupan, MODULES_MISSING
except ImportError:
    MODULES_MISSING = True


class TestTupanInterface(TestWithMPI):

    def test01(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Test TupanInterface initialization"
        instance = TupanInterface()
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test02(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Test TupanInterface new_particle / get_state"
        instance = TupanInterface()
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())

        id, error = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(0, id)
        id, error = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(1, id)
        self.assertEquals(0, instance.commit_particles())

        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
        self.assertEquals(0,  retrieved_state1['__result'])
        self.assertEquals(0,  retrieved_state2['__result'])
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals( 0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test03(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Test TupanInterface particle property getters/setters"
        instance = TupanInterface()
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([0, 0], instance.new_particle(0.01,  1, 0, 0,  0, 1, 0, 0.1).values())
        self.assertEquals([1, 0], instance.new_particle(0.02, -1, 0, 0,  0,-1, 0, 0.1).values())
####        self.assertEquals(-1, instance.get_mass(1)['__result']) # Have to commit first
        self.assertEquals(0, instance.commit_particles())

        # getters
        mass, result = instance.get_mass(0)
        self.assertAlmostEquals(0.01, mass)
        self.assertEquals(0,result)
        radius, result = instance.get_radius(1)
        self.assertAlmostEquals(0.1, radius)
        self.assertEquals(0,result)
        self.assertEquals(-1, instance.get_mass(2)['__result']) # Particle not found
        self.assertEquals([ 1, 0, 0,  0], instance.get_position(0).values())
        self.assertEquals([-1, 0, 0,  0], instance.get_position(1).values())
        self.assertEquals([ 0, 1, 0,  0], instance.get_velocity(0).values())
        self.assertEquals([ 0,-1, 0,  0], instance.get_velocity(1).values())

        # setters
        self.assertEquals(0, instance.set_state(0, 0.01, 1,2,3, 4,5,6, 0.1))
        self.assertEquals([0.01, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_mass(0, 0.02))
        self.assertEquals([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_radius(0, 0.2))
        self.assertEquals([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.2, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_position(0, 10,20,30))
        self.assertEquals([0.02, 10.0,20.0,30.0, 4.0,5.0,6.0, 0.2, 0], instance.get_state(0).values())
        self.assertEquals(0, instance.set_velocity(0, 40,50,60))
        self.assertEquals([0.02, 10.0,20.0,30.0, 40.0,50.0,60.0, 0.2, 0], instance.get_state(0).values())

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test04(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Test TupanInterface parameters"
        instance = TupanInterface()
        self.assertEquals(0, instance.initialize_code())

        self.assertEquals([0.03125, 0], instance.get_eta().values())
        self.assertEquals(0, instance.set_eta(0.001))
        self.assertEquals([0.001, 0], instance.get_eta().values())

        self.assertEquals([0.0, 0], instance.get_begin_time().values())
        self.assertEquals(0, instance.set_begin_time(1.0))
        self.assertEquals([1.0, 0], instance.get_begin_time().values())

        self.assertEquals(["sia21h.dkd", 0], instance.get_integrator_method().values())
        self.assertEquals(0, instance.set_integrator_method("bogus"))
        self.assertEquals(["bogus", 0], instance.get_integrator_method().values())
        self.assertEquals(-1, instance.commit_parameters())
        self.assertEquals(0, instance.set_integrator_method("sakura"))
        self.assertEquals(["sakura", 0], instance.get_integrator_method().values())

        self.assertEquals(0, instance.commit_parameters())

        self.assertEquals(0, instance.set_pn_order(7))
        self.assertEquals([7, 0], instance.get_pn_order().values())
        self.assertEquals(-1, instance.commit_parameters())
        self.assertEquals(0, instance.set_clight(1024))
        self.assertEquals([1024, 0], instance.get_clight().values())

        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test05(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Test TupanInterface evolve_model, binary"
        instance = TupanInterface()
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())

        self.assertEquals([0, 0], instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0, 0.001).values())
        self.assertEquals([1, 0], instance.new_particle(0.5, -0.5, 0, 0,  0,-0.5, 0, 0.001).values())
        self.assertEquals(0, instance.commit_particles())

        P = 2 * math.pi
        self.assertEquals(0, instance.evolve_model(P / 2)) # half an orbit
        for result, expected in zip(instance.get_position(0).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 2)

        self.assertEquals(0, instance.evolve_model(P)) # full orbit
        for result, expected in zip(instance.get_position(0).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 2)

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()



class TestTupan(TestWithMPI):

    default_converter = nbody_system.nbody_to_si(1.0e4 | units.MSun, 1.0 | units.AU)

    def new_sun_earth_system(self):
        particles = Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        return particles

    def test01(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan initialization"
        instance = Tupan(self.default_converter, )
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()

    def xtest02(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan parameters"
        instance = Tupan(self.default_converter, )
        instance.initialize_code()

        self.assertEquals(instance.parameters.epsilon_squared,
            instance.unit_converter.to_si(0.0 | nbody_system.length**2))
        self.assertEquals(instance.parameters.timestep_parameter, 0.125)

        for par, value in [('epsilon_squared_star_star', 0.0 | nbody_system.length**2),
                ('epsilon_squared_star_blackhole', 0.0 | nbody_system.length**2),
                ('epsilon_squared_blackhole_blackhole', 0.0 | nbody_system.length**2),
                ('initial_timestep_parameter', 1.0e-4),
                ('timestep_parameter_stars', 0.1),
                ('timestep_parameter_supermassive_black_holes', 0.4),
                ('timestep_parameter_intermediate_mass_black_holes', 0.4),
                ('max_relative_energy_error', 5.0e-5),
                ('maximum_timestep', 1.0/1024.0 | nbody_system.time),
                ('smbh_mass', 1.0 | nbody_system.mass)]:
            self.assertEquals(instance.unit_converter.to_si(value),
                getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEquals(instance.unit_converter.to_si(3.0 | value.unit),
                getattr(instance.parameters, par))

        # epsilon_squared is an alias for epsilon_squared_star_star, so epsilon_squared also has become 3:
        self.assertEquals(instance.parameters.epsilon_squared,
            instance.unit_converter.to_si(3.0 | nbody_system.length**2))
        instance.parameters.epsilon_squared = 0.1 | nbody_system.length**2
        self.assertEquals(instance.parameters.epsilon_squared,
            instance.unit_converter.to_si(0.1 | nbody_system.length**2))
        # timestep_parameter is an alias for timestep_parameter_stars, so timestep_parameter also has become 3:
        self.assertEquals(instance.parameters.timestep_parameter, 3.0)
        instance.parameters.timestep_parameter = 0.01
        self.assertEquals(instance.parameters.timestep_parameter, 0.01)

        self.assertEquals(instance.parameters.include_smbh, False)
        instance.parameters.include_smbh = True
        self.assertEquals(instance.parameters.include_smbh, True)
        self.assertEquals(instance.parameters.calculate_postnewtonian, True)
        instance.parameters.calculate_postnewtonian = False
        self.assertEquals(instance.parameters.calculate_postnewtonian, False)

        self.assertEquals(instance.parameters.drink, "Vodka martini. Shaken, not stirred.")

        instance.stop()

    def test03(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan particles"
        instance = Tupan(self.default_converter, )
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()

        self.assertAlmostEquals(instance.particles.mass, [1.0, 3.0037e-6] | units.MSun)
        self.assertAlmostEquals(instance.particles.radius, 1.0 | units.RSun)
        self.assertAlmostEquals(instance.particles.position,
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU)
        self.assertAlmostEquals(instance.particles.velocity,
            [[0.0, 0.0, 0.0], [0.0, 29.7885, 0.0]] | units.km / units.s, 3)

        instance.cleanup_code()
        instance.stop()

    def xtest04(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan evolve_model, 2 particles orbiting the SMBH"
        particles = Particles(2)
        particles.mass = 1.0 | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = ((constants.G * (10001.0 | units.MSun) / (1.0 | units.AU)).sqrt() +
                           (constants.G * (10000.0 | units.MSun) / (1.0 | units.AU)).sqrt())
        particles.move_to_center()
        print particles

        instance = Tupan(self.default_converter, )
        instance.initialize_code()
#        instance.parameters.include_smbh = True
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]

        P = 2 * math.pi * primary.x / primary.vy
        P_corrected = (P) * (2.0 / (1.0 + math.sqrt(1.0001)))

        position_at_start = primary.position.x
        instance.evolve_model(P_corrected / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 3)

        instance.evolve_model(P_corrected / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 3)

        instance.evolve_model(P_corrected)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)

        instance.cleanup_code()
        instance.stop()

    def test05(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan evolve_model, 2 particles, no SMBH"
        particles = Particles(2)
        particles.mass = 1.0 | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * (2.0 | units.MSun) / (2.0 | units.AU)).sqrt()
        particles.move_to_center()
        print particles

        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Tupan(converter, )
        instance.initialize_code()
        instance.parameters.integrator_method = "sakura"
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]

        P = 2 * math.pi * primary.x / primary.vy

        position_at_start = primary.position.x
        instance.evolve_model(P / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 3)

        instance.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 3)

        instance.evolve_model(P)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)

        instance.cleanup_code()
        instance.stop()

    def test06(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan evolve_model, earth-sun system, no SMBH"
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Tupan(converter, )
        instance.initialize_code()
#        instance.parameters.smbh_mass = 0.0 | units.MSun
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        earth = instance.particles[1]

        position_at_start = earth.position.x
        instance.evolve_model(0.25 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, earth.position.y, 3)

        instance.evolve_model(0.5 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, -earth.position.x, 3)

        instance.evolve_model(1.0 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, earth.position.x, 3)

        instance.cleanup_code()
        instance.stop()

    def test07(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing effect of Tupan parameter epsilon_squared"
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = self.new_sun_earth_system()
        particles.rotate(0.0, 0.0, -math.pi/4)
        particles.move_to_center()

        tan_initial_direction = particles[1].vy/particles[1].vx
        self.assertAlmostEquals(tan_initial_direction, math.tan(math.pi/4))
        tan_final_direction =  []
        for log_eps2 in range(-9,10,2):
            instance = Tupan(converter, )
            instance.initialize_code()
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
#            instance.parameters.smbh_mass = 0.0 | units.MSun
            instance.commit_parameters()
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            tan_final_direction.append(instance.particles[1].velocity[1]/
                instance.particles[1].velocity[0])
            instance.cleanup_code()
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEquals(tan_final_direction[0], math.tan(3 * math.pi / 4.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEquals(tan_final_direction[-1], tan_initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(tan_final_direction[i+1]-tan_final_direction[i]) for i in range(len(tan_final_direction)-1)]
        self.assertEquals(delta[len(tan_final_direction)/2 -1], max(delta))

    def xtest08(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan get_gravity_at_point and get_potential_at_point"
        instance = Tupan()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
#        instance.parameters.smbh_mass = 0.0 | nbody_system.mass

        particles = Particles(2)
        particles.mass = 1.0 | nbody_system.mass
        particles.radius =  0.0 | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)

        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, 1.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(zero, x1, zero, zero)

            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fy1, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz1, 0.0 | nbody_system.acceleration)

            self.assertAlmostEqual(fx0, -1.0 * fx1)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0)
            self.assertAlmostEqual(potential0, potential1)

        instance.stop()

    def test09(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan evolve_model and getters energy, plummer sphere, no SMBH"
        converter = nbody_system.nbody_to_si(1.0e2 | units.MSun, 1.0 | units.parsec)
        instance = Tupan(converter, )
        instance.parameters.timestep_parameter = 1.0/256
        instance.initialize_code()
#        instance.parameters.smbh_mass = 0.0 | units.MSun
        instance.commit_parameters()
        numpy.random.seed(987654321)
        instance.particles.add_particles(new_plummer_model(100, convert_nbody=converter))
        instance.commit_particles()

        kinetic_energy = instance.kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 2.12292810174e+37 | units.J, 10)
        self.assertAlmostRelativeEqual(potential_energy, -4.33511391248e+37 | units.J, 10)

        initial_total_energy = kinetic_energy + potential_energy
        instance.evolve_model(0.1 | nbody_system.time)
        kinetic_energy = instance.kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 2.1362368884e+37 | units.J, 4)
        self.assertAlmostRelativeEqual(potential_energy, -4.34842269914e+37 | units.J, 4)

        self.assertAlmostRelativeEqual(potential_energy + kinetic_energy, initial_total_energy, 4)

        instance.cleanup_code()
        instance.stop()

    def xtest10(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan collision_detection"
        particles = Particles(7)
        particles.mass = 0.00000001 | nbody_system.mass
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed

        instance = Tupan()
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)

        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
        self.assertEquals(len(collisions.particles(0)), 3)
        self.assertEquals(len(collisions.particles(1)), 3)
        self.assertEquals(len(particles - collisions.particles(0) - collisions.particles(1)), 1)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) <
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True, True, True])

        sticky_merged = Particles(len(collisions.particles(0)))
        sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
        sticky_merged.radius = collisions.particles(0).radius
        for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
            merged.position = (p1 + p2).center_of_mass()
            merged.velocity = (p1 + p2).center_of_mass_velocity()

        print instance.model_time
        print instance.particles
        instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
        instance.particles.add_particles(sticky_merged)

        instance.evolve_model(1.0 | nbody_system.time)
        print
        print instance.model_time
        print instance.particles
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEquals(len(collisions.particles(0)), 1)
        self.assertEquals(len(collisions.particles(1)), 1)
        self.assertEquals(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) <
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()

    def test11(self):
        if MODULES_MISSING:
            self.skip("Failed to import a module required for Tupan")
        print "Testing Tupan properties"
        numpy.random.seed(12345)
        particles = new_plummer_model(100, do_scale=True)
        particles.position += [1, 2, 3] | nbody_system.length
        cluster_velocity = [4, 5, 6] | nbody_system.speed
        particles.velocity += cluster_velocity
        external_kinetic_energy = (0.5 | nbody_system.mass) * cluster_velocity.length_squared()

        instance = Tupan()
        instance.particles.add_particles(particles)

        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 0.25 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(potential_energy, -0.5 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position,
            [1, 2, 3] | nbody_system.length, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity,
            [4, 5, 6] | nbody_system.speed, 10)
        initial_total_energy = kinetic_energy + potential_energy

        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEqual(instance.model_time, 0.1 | nbody_system.time, 3)
        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy+potential_energy, -0.25 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position,
            [1.4, 2.5, 3.6] | nbody_system.length, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity,
            [4, 5, 6] | nbody_system.speed, 3)

        instance.cleanup_code()
        instance.stop()

