import numpy

from amuse.test.amusetest import TestCase
from amuse.support.exceptions import AmuseWarning, AmuseException
from amuse.units import units, nbody_system, generic_unit_system
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.units.quantities import zero, VectorQuantity
from amuse.datamodel import Particle, Particles, ParticlesSuperset
from amuse.community.interface.gd import GravitationalDynamics
from amuse.support.codes.stopping_conditions import StoppingConditions
from amuse.ext.sink import SinkParticles, new_sink_particles

from amuse.community.gadget2.interface import Gadget2
from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.ic.plummer import new_plummer_model


class TestSinkParticles(TestCase):
    
    def test1(self):
        print "Testing SinkParticles initialization from new (blank) particle"
        sinks = SinkParticles(Particles(2), sink_radius=[42.0,43.0]|units.RSun)
        
        self.assertEqual(sinks.sink_radius, [42.0, 43.0] | units.RSun)
        self.assertEqual(sinks.mass, 0.0 | units.MSun)
        self.assertEqual(sinks.position, [0.0, 0.0, 0.0] | units.parsec)
        self.assertEqual(sinks.velocity, [0.0, 0.0, 0.0] | units.km/units.s)
        
        sinks = SinkParticles(Particles(2), sink_radius=24.0|units.RSun, mass=[1.0,2.0]|units.MSun, 
            position=[1,2,3]|nbody_system.length, velocity=[4,5,6]|nbody_system.length/nbody_system.time)
        self.assertEqual(sinks.sink_radius, 24.0 | units.RSun)
        self.assertEqual(sinks.mass, [1.0, 2.0] | units.MSun)
        self.assertEqual(sinks.position, [1.0, 2.0, 3.0] | nbody_system.length)
        self.assertEqual(sinks.velocity, [4.0, 5.0, 6.0] | nbody_system.length/nbody_system.time)
    
    def test2(self):
        print "Testing SinkParticles initialization from existing particles"
        original = Particles(3)
        self.assertRaises(AttributeError, SinkParticles, original, expected_message=
            "You tried to access attribute 'radius' but this attribute is not defined for this set.")
        original.radius = 42.0 | units.RSun
        original.mass = 10.0 | units.MSun
        original.position = [[i, -i, 2*i] for i in range(3)] | units.parsec
        sinks = SinkParticles(original)
        self.assertEqual(sinks.sink_radius, 42.0 | units.RSun)
        self.assertEqual(sinks.mass, 10.0 | units.MSun)
        self.assertEqual(sinks.position, [[0,0,0], [1,-1,2], [2,-2,4]] | units.parsec)
        self.assertRaises(AttributeError, getattr, sinks, "bogus", expected_message=
            "You tried to access attribute 'bogus' but this attribute is not defined for this set.")
    
    def test3(self):
        print "Testing SinkParticles initialization from existing particles in set"
        particles = Particles(10)
        self.assertRaises(AttributeError, SinkParticles, particles[[4, 7]], expected_message=
            "You tried to access attribute 'radius' but this attribute is not defined for this set.")
        particles.radius = 42.0 | units.RSun
        particles.mass = range(1,11) | units.MSun
        particles.position = [[i, 2*i, 3*i] for i in range(10)] | units.parsec
        
        sinks = SinkParticles(particles[[4]])
        self.assertEqual(sinks.mass, 5.0 | units.MSun)
        self.assertEqual(sinks.sink_radius, 42.0 | units.RSun)
        self.assertEqual(sinks.radius, 42.0 | units.RSun)
        self.assertEqual(sinks.position, [4.0, 8.0, 12.0] | units.parsec)
        
        sinks = SinkParticles(particles[[4, 7]], sink_radius=[1,2]|units.AU)
        self.assertEqual(sinks.sink_radius, [1.0, 2.0] | units.AU)
        self.assertEqual(sinks.radius, 42.0 | units.RSun)
        self.assertEqual(sinks.mass, [5.0, 8.0] | units.MSun)
        self.assertEqual(sinks.position, [[4, 8, 12], [7, 14, 21]] | units.parsec)
        
        self.assertEqual(set(['key', 'mass', 'radius', 'x', 'y', 'z', 'sink_radius', 'vx','vy','vz']), 
            set(str(sinks).split("\n")[0].split()))
        self.assertEqual(set(['key', 'mass', 'radius', 'x', 'y', 'z']), 
            set(str(particles).split("\n")[0].split()))
    
    def test4(self):
        print "Testing SinkParticles accrete"
        particles = Particles(10)
        particles.radius = 42.0 | units.RSun
        particles.mass = range(1,11) | units.MSun
        particles.position = [[i, 2*i, 3*i] for i in range(10)] | units.parsec
        particles.velocity = [[i, 0, -i] for i in range(10)] | units.km/units.s
        particles.age = range(10) | units.Myr
        copy = particles.copy()
        
        sinks = SinkParticles(particles[[3, 7]], sink_radius=[4,5]|units.parsec)
        self.assertEqual(sinks.sink_radius, [4.0, 5.0] | units.parsec)
        self.assertEqual(sinks.mass, [4.0, 8.0] | units.MSun)
        self.assertEqual(sinks.position, [[3, 6, 9], [7, 14, 21]] | units.parsec)
        
        sinks.accrete(particles)
        self.assertEqual(len(particles), 6) # 4 particles were accreted
        self.assertEqual(sinks.mass, [12.0, 24.0] | units.MSun) # mass of sinks increased
        self.assertEqual(sinks.get_intersecting_subset_in(particles).mass, 
            [12.0, 24.0] | units.MSun) # original particles' masses match
        self.assertEqual(particles.total_mass(), copy.total_mass()) # total mass is conserved
        self.assertEqual(particles.center_of_mass(), copy.center_of_mass()) # center of mass is conserved
        self.assertEqual(particles.center_of_mass_velocity(), copy.center_of_mass_velocity()) # center of mass velocity is conserved
        
        sinks.sink_radius = [4.0, 8.0] | units.parsec
        sinks.accrete(particles)
        self.assertEqual(len(particles), 4) # another 2 particles were accreted
        self.assertEqual(sinks.mass, [12.0, 40.0] | units.MSun) # mass of sinks increased
        self.assertEqual(sinks.get_intersecting_subset_in(particles).mass, 
            [12.0, 40.0] | units.MSun) # original particles' masses match
        self.assertEqual(particles.total_mass(), copy.total_mass()) # total mass is conserved
        self.assertEqual(particles.center_of_mass(), copy.center_of_mass()) # center of mass is conserved
        self.assertEqual(particles.center_of_mass_velocity(), copy.center_of_mass_velocity()) # center of mass velocity is conserved
    
    def test5(self):
        print "Testing SinkParticles accrete, one particle within two sinks' radii"
        particles = Particles(10)
        particles.radius = 42.0 | units.RSun
        particles.mass = range(1,11) | units.MSun
        particles.position = [[i, 2*i, 3*i] for i in range(10)] | units.parsec
        particles.velocity = [[i, 0, -i] for i in range(10)] | units.km/units.s
        particles.age = range(10) | units.Myr
        copy = particles.copy()
        
        sinks = SinkParticles(particles[[3, 7]], sink_radius=[4,12]|units.parsec)
        self.assertEqual(sinks.sink_radius, [4.0, 12.0] | units.parsec)
        self.assertEqual(sinks.mass, [4.0, 8.0] | units.MSun)
        self.assertEqual(sinks.position, [[3, 6, 9], [7, 14, 21]] | units.parsec)
        
        sinks.accrete(particles)
        self.assertEqual(len(particles), 4) # 6 particles were accreted
        self.assertEqual(sinks.mass, [12.0, 40.0] | units.MSun) # mass of sinks increased
        self.assertEqual(sinks.get_intersecting_subset_in(particles).mass, 
            [12.0, 40.0] | units.MSun) # original particles' masses match
        self.assertEqual(particles.total_mass(), copy.total_mass()) # total mass is conserved
        self.assertEqual(particles.center_of_mass(), copy.center_of_mass()) # center of mass is conserved
        self.assertEqual(particles.center_of_mass_velocity(), copy.center_of_mass_velocity()) # center of mass velocity is conserved
    

class TestNewSinkParticles(TestCase):
    
    def test1(self):
        print "Test the documentation for new_sink_particles"
        print new_sink_particles.__doc__
    
    def test2(self):
        print "Demonstrate new_sink_particles usage"
        cloud = Particles(100)
        cloud.mass = 1 | units.MSun
        cloud.position = [[0, 0, 0], [100, 100, 100], [200, 200, 200], [300, 300, 300]]*25 | units.parsec
        cloud.velocity = [[0, 0, 0], [1, 1, 1]]*50 | units.km / units.s
        unit_converter = ConvertBetweenGenericAndSiUnits(1|units.m, 1|units.kg, 1|units.s)
        sph_code = Stub(unit_converter)
        sph_code.parameters.stopping_condition_maximum_density = 1 | units.kg / units.m**3
        sph_code.gas_particles.add_particles(cloud)
        density_limit_detection = sph_code.stopping_conditions.density_limit_detection
        density_limit_detection.enable()
        
        sph_code.evolve_model(1 | units.Myr)
        self.assertTrue(density_limit_detection.is_set())
        self.assertEqual(len(density_limit_detection.particles()), 3)
        self.assertEqual(density_limit_detection.particles().position, 
            [[100, 100, 100], [200, 200, 200], [300, 300, 300]] | units.parsec)
        print density_limit_detection.particles()
        
        clumps = density_limit_detection.particles().copy()
        sph_code.gas_particles.remove_particles(clumps)
        
        sinks = new_sink_particles(clumps, sink_radius=1|units.parsec)
        self.assertEqual(sinks.sink_radius, 1.0 | units.parsec)
        self.assertEqual(sinks.mass, 1.0 | units.MSun)
        self.assertEqual(sinks.position, 
            [[100, 100, 100], [200, 200, 200], [300, 300, 300]] | units.parsec)
        self.assertEqual(len(sph_code.gas_particles), 97)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass() + clumps.total_mass(), 100 | units.MSun, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), 97 | units.MSun, 10)
        
        sinks.accrete(sph_code.gas_particles)
        self.assertAlmostRelativeEqual(sinks.mass, [25, 25, 25] | units.MSun, 10)
        self.assertEqual(len(sph_code.gas_particles), 25)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass() + clumps.total_mass(), 100 | units.MSun, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), 25 | units.MSun, 10)
    
    def test3(self):
        print "Demonstrate new_sink_particles usage (using Gadget2)"
        UnitLength = 1.0 | units.kpc
        UnitMass = 1.0e10 | units.MSun
        UnitVelocity = 1.0 | units.km / units.s
        convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
        converter = ConvertBetweenGenericAndSiUnits(UnitLength, UnitMass, UnitVelocity)
        number_gas_particles = 1000
        gas = new_evrard_gas_sphere(number_gas_particles, convert_nbody, do_scale=True, seed=12345)
        
        sph_code = Gadget2(converter)
        sph_code.initialize_code()
        sph_code.parameters.stopping_condition_maximum_density = 10 * UnitMass / UnitLength**3
        sph_code.gas_particles.add_particles(gas)
        self.assertIsOfOrder(max(sph_code.gas_particles.density), UnitMass / UnitLength**3)
        
        density_limit_detection = sph_code.stopping_conditions.density_limit_detection
        density_limit_detection.enable()
        
        sph_code.evolve_model(10.0 | units.Myr)
        self.assertTrue(density_limit_detection.is_set())
        self.assertTrue(sph_code.model_time < 10.0 | units.Myr)
        print "density_limit exceeded at t =", sph_code.model_time.as_quantity_in(units.Myr)
        self.assertEquals(len(density_limit_detection.particles()), 1)
        self.assertTrue(density_limit_detection.particles().density > 
                10 * UnitMass / UnitLength**3)
        
        clumps = density_limit_detection.particles().copy()
        sph_code.gas_particles.remove_particles(clumps)
        clumps_in_code = sph_code.dm_particles.add_particles(clumps)
        
        sinks = new_sink_particles(clumps_in_code)
        self.assertEqual(sinks.sink_radius, clumps.radius)
        self.assertAlmostRelativeEqual(sinks.mass, UnitMass / number_gas_particles, 10)
        self.assertAlmostRelativeEqual(sinks.position, clumps.position, 10)
        self.assertEqual(len(sph_code.gas_particles), number_gas_particles - 1)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass(), UnitMass, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), UnitMass - sinks.total_mass(), 10)
        self.assertEqual(set(sinks.get_attribute_names_defined_in_store()) - set(["sink_radius"]), 
            set(sph_code.particles.get_attribute_names_defined_in_store()))
        
        sinks.accrete(sph_code.gas_particles)
        self.assertAlmostRelativeEqual(sinks.mass, 3 * UnitMass / number_gas_particles, 10)
        self.assertEqual(len(sph_code.gas_particles), number_gas_particles - 3)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass(), UnitMass, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), UnitMass - sinks.total_mass(), 10)
        
        sinks.accrete(sph_code.particles) # Nothing happens: gas already gone, and cannot accrete itself
        self.assertAlmostRelativeEqual(sinks.mass, 3 * UnitMass / number_gas_particles, 10)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass(), UnitMass, 10)
        
        steps = 0
        while True:
            sph_code.evolve_model(sph_code.model_time + (0.1 | units.Myr))
            sinks.sink_radius = 4 * clumps_in_code.radius
            sinks.accrete(sph_code.gas_particles)
            steps += 1
            if density_limit_detection.is_set():
                break
        
        self.assertEqual(len(sph_code.gas_particles), number_gas_particles - 7)
        self.assertAlmostRelativeEqual(sinks.mass, 7 * UnitMass / number_gas_particles, 10)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass(), UnitMass, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), UnitMass - sinks.total_mass(), 10)
        
        self.assertTrue(density_limit_detection.is_set())
        self.assertEqual(steps, 5)
        self.assertTrue(sph_code.model_time < 10.0 | units.Myr)
        print "density_limit exceeded at t =", sph_code.model_time.as_quantity_in(units.Myr)
        self.assertEquals(len(density_limit_detection.particles()), 5)
        self.assertTrue((density_limit_detection.particles().density > 
                10 * UnitMass / UnitLength**3).all())
        
        clumps = density_limit_detection.particles().copy()
        sph_code.gas_particles.remove_particles(clumps)
        clumps_in_code = sph_code.dm_particles.add_particles(clumps)

        sinks.add_sinks(clumps_in_code, sink_radius=0.1|units.kpc)
        self.assertEqual(sinks[1:].sink_radius, 0.1 | units.kpc)
        self.assertEqual(len(sph_code.gas_particles), number_gas_particles - 12)
        self.assertAlmostRelativeEqual(sinks.mass[1:], UnitMass / number_gas_particles, 10)
        self.assertAlmostRelativeEqual(sinks.position[1:], clumps.position, 10)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass(), UnitMass, 10)
        self.assertAlmostRelativeEqual(sph_code.gas_particles.total_mass(), UnitMass - sinks.total_mass(), 10)
        
        sinks.accrete(sph_code.gas_particles)
        self.assertEqual(len(sph_code.gas_particles), number_gas_particles - 66)
        self.assertAlmostRelativeEqual(sph_code.particles.total_mass().as_quantity_in(units.MSun), UnitMass, 10)
        self.assertAlmostRelativeEqual(sinks.mass, [7.0, 13.0, 15.0, 9.0, 11.0, 11.0] * UnitMass / number_gas_particles, 10)
    


class StubInterface(object):
    
    def __init__(self, **options):
        self.maximum_density = 1 | units.kg / units.m**3
        self._gas_particles = Particles()
        self._dm_particles = Particles()
        self._all_particles = ParticlesSuperset([self._gas_particles, self._dm_particles])
    
    def before_get_parameter(self):
        pass
        
    def before_set_parameter(self):
        pass
        
    def initialize_code(self):
        return 0
    
    synchronize_model = commit_particles = recommit_particles = commit_parameters = initialize_code
    
    def new_particle(self, mass, x, y, z, vx, vy, vz, *args):
        next_id = len(self._dm_particles)
        temp = Particles(len(mass))
        temp.mass = mass
        temp.x = x
        temp.y = y
        temp.z = z
        temp.vx = vx
        temp.vy = vy
        temp.vz = vz
        temp.id = list(range(next_id, next_id + len(mass)))
        self._dm_particles.add_particles(temp)
        return [temp.id, temp.id]
    
    def new_gas_particle(self, mass, x, y, z, vx, vy, vz, *args):
        next_id = len(self._gas_particles) + 1000000
        temp = Particles(len(mass))
        temp.mass = mass
        temp.x = x
        temp.y = y
        temp.z = z
        temp.vx = vx
        temp.vy = vy
        temp.vz = vz
        temp.id = list(range(next_id, next_id + len(mass)))
        self._gas_particles.add_particles(temp)
        return [temp.id, temp.id]
    
    def delete_particle(self, indices):
        for index in indices:
            for id, particle in zip(self._all_particles.id, self._all_particles):
                if id == index:
                    self._all_particles.remove_particle(particle)
        return 0
    
    def get_mass(self, indices):
        return [[mass for index in indices for id, mass in zip(self._all_particles.id, 
            self._all_particles.mass) if index == id], [0]*len(indices)]
    
    def set_mass(self, indices, masses):
        for index, mass in zip(indices, masses):
            for id, particle in zip(self._all_particles.id, self._all_particles):
                if id == index:
                    particle.mass = mass
                    break
        return 0
    
    def get_position(self, indices):
        return [[x for index in indices for id, x in zip(self._all_particles.id, self._all_particles.x) if index == id], 
            [y for index in indices for id, y in zip(self._all_particles.id, self._all_particles.y) if index == id], 
            [z for index in indices for id, z in zip(self._all_particles.id, self._all_particles.z) if index == id], 
            [0]*len(indices)]
    
    def get_velocity(self, indices):
        return [[vx for index in indices for id, vx in zip(self._all_particles.id, self._all_particles.vx) if index == id], 
            [vy for index in indices for id, vy in zip(self._all_particles.id, self._all_particles.vy) if index == id], 
            [vz for index in indices for id, vz in zip(self._all_particles.id, self._all_particles.vz) if index == id], 
            [0]*len(indices)]
    
    def has_stopping_condition(self, type):
        return 1 if type == 6 else 0
    
    def get_stopping_condition_maximum_density_parameter(self):
        return self.maximum_density
    
    def set_stopping_condition_maximum_density_parameter(self, value):
        self.maximum_density = value
    
    is_stopping_condition_set = is_stopping_condition_enabled = has_stopping_condition
    
    def get_number_of_stopping_conditions_set(self):
        return 3
    
    def get_stopping_condition_info(self, sc_indices):
        return [6]*len(sc_indices), [1]*len(sc_indices)
    
    def get_stopping_condition_particle_index(self, sc_index, sc_sub_index):
        return range(len(self._gas_particles) + 1000000 - len(sc_index), len(self._gas_particles) + 1000000)
    
    def enable_stopping_condition(self, type):
        pass
    
    def evolve_model(self, time):
        return 0
    

class Stub(GravitationalDynamics):
    
    def __init__(self, unit_converter = None, **options):
        self.stopping_conditions = StoppingConditions(self)
        
        GravitationalDynamics.__init__(
            self,
            StubInterface(**options),
            unit_converter,
            **options
        )
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_stopping_condition_maximum_density_parameter",
            "set_stopping_condition_maximum_density_parameter", 
            "stopping_condition_maximum_density", 
            "maximum density of a gas particle", 
            default_value = -1.0 | generic_unit_system.density
        )
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method("new_gas_particle", 
            (nbody_system.mass, nbody_system.length, nbody_system.length, nbody_system.length,
                nbody_system.speed, nbody_system.speed, nbody_system.speed), 
            (object.INDEX, object.ERROR_CODE))
    
    def define_particle_sets(self, object):
        object.define_super_set('particles', ['dm_particles','gas_particles'], 
            index_to_default_set = 0)
        
        object.define_set('dm_particles', 'index_of_the_particle')
        object.set_new('dm_particles', 'new_particle')
        object.set_delete('dm_particles', 'delete_particle')
        object.add_getter('dm_particles', 'get_mass', names=("mass",))
        object.add_setter('dm_particles', 'set_mass', names=("mass",))
        object.add_getter('dm_particles', 'get_position', names=("x", "y", "z"))
        object.add_getter('dm_particles', 'get_velocity', names=("vx", "vy", "vz"))
        
        object.define_set('gas_particles', 'index_of_the_particle')
        object.set_new('gas_particles', 'new_gas_particle')
        object.set_delete('gas_particles', 'delete_particle')
        object.add_getter('gas_particles', 'get_mass', names=("mass",))
        object.add_getter('gas_particles', 'get_position', names=("x", "y", "z"))
        object.add_getter('gas_particles', 'get_velocity', names=("vx", "vy", "vz"))
        
        object.add_query('particles', 'get_stopping_condition_particle_index')
    

