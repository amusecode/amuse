from amuse.test import amusetest

from amuse.support.exceptions import AmuseException



from amuse.support.interface import InCodeComponentImplementation

import numpy
import time
import sys
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support import data
class TestParticles(amusetest.TestCase):
    
    def test1(self):
        
        particles = data.Particles(2)
        for i in range(3):
            particles.mass = (i * 1.0) | units.kg
            particles.savepoint((i + 1) * 1.0 | units.s)
        
        masses = particles[0].get_timeline_of_attribute("mass")
        self.assertEquals(len(masses), 3)
        self.assertEquals(masses[0][0], 1.0 | units.s)
        self.assertEquals(masses[1][0], 2.0 | units.s)
        self.assertEquals(masses[2][0], 3.0 | units.s)

    def test2(self):
        particles = data.Particles(2)
        particles.mass = [1,1]|units.kg
        heavy = particles.select(lambda m: m>1|units.kg, ["mass"])
        self.assertTrue(heavy.is_empty())
        
    def test3(self):
        particles = data.Particles(2)
        for i in range(3):
            particles.mass = (i * 2.0) | units.kg
            particles.savepoint((i + 1) * 1.0 | units.s)
        
        state0 = particles.get_state_at_timestamp(2.0 | units.s)
        self.assertEquals(state0[0].mass, 2.0 | units.kg)
        state0 = particles.get_state_at_timestamp(2.4 | units.s)
        self.assertEquals(state0[0].mass, 2.0 | units.kg)
        state0 = particles.get_state_at_timestamp(2.6 | units.s)
        self.assertEquals(state0[0].mass, 4.0 | units.kg)
        
    
    def test4(self):
        
        particles = data.Particles(2)
        particles.mass = 1.0 | units.kg
        particles.vy = 1.0 | units.m / units.s
        particles.vx = 0.0 | units.m / units.s
        particles.vz = 0.0 | units.m / units.s
        particles.x = [0.0, 1.0] | units.m
        particles.y = 0.0 | units.m
        particles.z = 0.0 | units.m
        self.assertEquals(particles.kinetic_energy(), 1.0 | units.J)
        self.assertEquals(particles.potential_energy(), -1.0 * constants.G * (1.0 | units.kg ** 2 / units.m))
        self.assertEquals(particles.center_of_mass().x, 0.5 | units.m)
        self.assertEquals(particles.center_of_mass().y, 0.0 | units.m)
        self.assertEquals(particles.center_of_mass().z, 0.0 | units.m)
        self.assertEquals(particles.center_of_mass_velocity().x, 0.0 | units.m / units.s)
        self.assertEquals(particles.center_of_mass_velocity().y, 1.0 | units.m / units.s)
        self.assertEquals(particles.center_of_mass_velocity().z, 0.0 | units.m / units.s)
        
    def test5(self):
        
        particles = data.Particles(4)
        particles.mass = 1.0 | units.kg
        particles[2].mass = 2.0 | units.kg
        subset = particles.select_array(lambda mass: mass > 1.0 | units.kg, ["mass"])
        self.assertEquals(len(subset), 1)
        copyof_subset = subset.copy()
        self.assertEquals(len(copyof_subset), 1)
        

    def test6(self):
        
        particles = data.Particles(2)
        particles.mass = 1.0 | units.kg
        particles.vy = 1.0 | units.m / units.s
        particles.vx = 0.0 | units.m / units.s
        particles.vz = 0.0 | units.m / units.s
        particles.x = [0.0, 1.0] | units.m
        particles.y = 0.0 | units.m
        particles.z = 0.0 | units.m
        
        Ek=0. | units.J
        Ep=0. | units.J
        for x in particles:
            Ek+=x.mass*x.specific_kinetic_energy()
            Ep+=x.mass*x.potential()/2
        self.assertEquals(particles.kinetic_energy(), Ek)
        self.assertEquals(particles.potential_energy(), Ep)
        
        self.assertEquals((particles.mass*particles.specific_kinetic_energy()).sum(), Ek)
        self.assertEquals(0.5*(particles.mass*particles.potential()).sum(), Ep)

    
    def test7(self):
        
        particles = data.Particles(10)
        for i in range(10):
            particles[i].mass = (i * 1.0) | units.kg

    def test8(self):
        particles = data.Particles(4)
        particles.mass = [1.,2.,3.,4.] | units.kg
        subset = particles[0:2]
        self.assertEquals(len(subset), 2)
        self.assertTrue(str(subset).find('kg') > 0)

    def test9(self):
        particles = data.Particles(4)
        particles.mass = [1.,2.,3.,4.] | units.kg
        particles.vy = [1.,2.,3.,4.] | units.m / units.s
        particles.vx = [0.,1.,2.,3.] | units.m / units.s
        particles.vz = [1.,2.,3.,4.] | units.m / units.s
        particles.x = [0., 1.,2.,3.] | units.m
        particles.y = [4., 3.,2.,1.] | units.m
        particles.z = [4., 3.,2.,1.] | units.m
        
        Ek=0. | units.J
        Ep=0. | units.J
        for x in particles:
            Ek+=x.mass*x.specific_kinetic_energy()
            Ep+=x.mass*x.potential()/2
        self.assertEquals(particles.kinetic_energy(), Ek)
        self.assertEquals(particles.potential_energy(), Ep)
        
        self.assertEquals((particles.mass*particles.specific_kinetic_energy()).sum(), Ek)
        self.assertEquals(0.5*(particles.mass*particles.potential()).sum(), Ep)

       
class TestStars(amusetest.TestCase):

    def test1(self):
        stars = data.Stars(2)
        stars[0].mass = 10 | units.g
        stars[0].position = units.m(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].position = units.m(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m, stars.center_of_mass().x)

    def test2(self):
        stars = data.Stars(2)
        stars[0].mass = 10 | units.g
        stars[0].velocity = (units.m / units.s)(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].velocity = (units.m / units.s)(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m / units.s, stars.center_of_mass_velocity().x)
        self.assertEquals(1.0 | units.m / units.s, stars.center_of_mass_velocity().y)
        
    
    def test3(self):
        stars = data.Stars(2)
        stars[0].mass = 10 | units.g
        stars[0].velocity = (units.m / units.s)(numpy.array([1.0,2.0,1.0])) 
        stars[0].position = units.m(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].velocity = (units.m / units.s)(numpy.array([0.0,0.0,0.0]))
        stars[1].position = units.m(numpy.array([0.0,0.0,0.0]))
        
        
        self.assertEquals(stars.mass[0], 10|units.g)
        
        self.assertEquals(stars.position[0], [1.0, 2.0, 1.0] | units.m)
        self.assertEquals(stars.velocity[0], [1.0, 2.0, 1.0] | units.m / units.s)


    
    def test4(self):
        stars = data.Stars(2)
        stars[0].x = 1.0  | units.km
        stars[0].y = 2000.0 | units.m
        stars[0].z = 3500.0 | units.m
        
        self.assertEquals(stars.position[0], [1000.0, 2000.0, 3500.0] | units.m)    



class TestParticlesWithBinding(amusetest.TestCase):
    class TestLegacyCode(object):
            
        def __init__(self):
            self.masses = {}
            
        def get_mass(self, id):
            masses = []
            errors = []
            for x in id:
                masses.append(self.masses[x])
                errors.append(0)
            return ( masses, errors, )
        
        def set_mass(self, id, mass):
            for i,m in zip(id,mass):
                self.masses[i] = m
                
            return ( [0] * len(id),)
            
        def new_particle(self, mass):
            ids = []
            errors = []
            
            for x in mass:
                id = len(self.masses)
                self.masses[len(self.masses)]  = x
                ids.append(id)
                errors.append(0)
                
            return (ids, errors)
        
        def delete_particle(self, ids):
            errors = []
            for x in ids:
                del self.masses[x]
                errors.append(0)
            return errors
            
        def get_number_of_particles(self):
            return (len(self.masses), 0)
            
        set_state = set_mass
        get_state = get_mass
        
    class TestInterface(InCodeComponentImplementation):
        
        def __init__(self):
            InCodeComponentImplementation.__init__(self, TestParticlesWithBinding.TestLegacyCode())
        
        def define_methods(self, handler):
            handler.add_method('get_mass',(handler.NO_UNIT,), (units.g, handler.ERROR_CODE))
            handler.add_method('set_mass',(handler.NO_UNIT, units.g,), (handler.ERROR_CODE,))
            handler.add_method('new_particle',(units.g,), (handler.INDEX, handler.ERROR_CODE))
            handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
            handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        def define_particle_sets(self, handler):
            handler.define_set('particles', 'id')
            handler.set_new('particles', 'new_particle')
            handler.set_delete('particles', 'delete_particle')
            handler.add_setter('particles', 'set_mass')
            handler.add_getter('particles', 'get_mass', names = ('mass',))
        
    
    def test1(self):
        interface = self.TestInterface()
        interface.particles.add_particles_to_store(
            [1,2],
            ["mass"],
            [[3.0, 4.0] | units.kg]
        )
        
        remote_particles = interface.particles
        local_particles = remote_particles.copy()
        
        channel = remote_particles.new_channel_to(local_particles)
        
        self.assertEquals(interface.masses[0], 3000)
        
        interface.masses[0] = 3500
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3.5)
        self.assertEquals(local_particles[0].mass.value_in(units.kg), 3.0)
        
        channel.copy_attributes(["mass"])
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3.5)
        self.assertEquals(local_particles[0].mass.value_in(units.kg), 3.5)
        
        
    
    
    def test2(self):
        interface = self.TestInterface()
        interface.particles.add_particles_to_store(
            [1, 2],
            ["mass"],
            [[3.0, 4.0] | units.kg]
        )
        
        remote_particles = interface.particles
        
        self.assertEquals(len(remote_particles), 2)
        
        interface.particles.add_particles_to_store(
            [3, 4],
            ["mass"],
            [[5.0, 6.0] | units.kg]
        )
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        interface.particles.remove_particles_from_store((1,3))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
    
    def test3(self):
        interface = self.TestInterface()
        
        local_particles = data.Stars(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        local_particles2 = data.Stars(2)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0])
       
        remote_particles.add_particles(local_particles2)
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        keys = remote_particles.get_all_keys_in_store()
        interface.particles.remove_particles_from_store((keys[0], keys[2]))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
        
    
    def test4(self):
        interface = self.TestInterface()
        
        local_particles = data.Stars(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        particle = data.Particle()
        particle.mass = units.g.new_quantity(10.0)
        local_particles.add_particle(particle)
        
        self.assertEquals(len(local_particles), 3)
        
    
    
    def test5(self):
        interface = self.TestInterface()
        
        local_particles1 = data.Stars(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particles2 = data.Stars(3)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0, 7.0])
        
        local_particles1.add_particles(local_particles2)
        
        self.assertEquals(len(local_particles1), 5)
        self.assertEquals(len(local_particles2), 3)
        
        local_particles1.synchronize_to(remote_particles)
        
        local_particles1.remove_particles_from_store([local_particles1.get_all_keys_in_store()[0]])
        local_particles1.synchronize_to(remote_particles)
        self.assertEquals(len(remote_particles), 4)
        
    
    def test6(self):
        interface = self.TestInterface()
        
        local_particles1 = data.Stars(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particle = data.Particle()
        local_particle.mass = units.kg.new_quantity(5.0)
        
        local_particles1.add_particle(local_particle)
        
        self.assertEquals(len(local_particles1), 3)
        self.assertEquals(len(remote_particles), 2)
        
        local_particles1.synchronize_to(remote_particles)
        
        self.assertEquals(len(remote_particles), 3)
        
        remote_particles.remove_particle(local_particle)
        self.assertEquals(len(remote_particles), 2)
        
        
        
        
        
        
        
class TestParticlesWithUnitsConverted(amusetest.TestCase):
    
    def test1(self):
        stars = data.Stars(2)
        stars[0].mass = 10 | units.g
        stars[1].mass = 20 | units.g
        
        class LengthMassConverter(object):
            "source == length, m"
            "target == mass, g"
            source_unit = units.m
            target_unit = units.g
            
            def from_source_to_target(self, quantity):
                value = quantity.value_in(self.source_unit)
                return self.target_unit.new_quantity(value) 
                
            def from_target_to_source(self, quantity):
                value = quantity.value_in(self.target_unit)
                return self.source_unit.new_quantity(value) 
        
        converted_stars = data.ParticlesWithUnitsConverted(stars, LengthMassConverter())
        
        self.assertEquals(stars[0].mass, 10 | units.g)
        self.assertEquals(converted_stars[0].mass, 10 | units.m)
        
        
        converted_stars[0].mass = 30 | units.m
        
        self.assertEquals(stars[0].mass, 30 | units.g)
        
    
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        stars = data.Stars(1)
        stars[0].mass = 10 | nbody_system.mass
        stars[0].x = 10.0 | nbody_system.length
        stars[0].y = 20.0 | nbody_system.length
        stars[0].z = 30.0 | nbody_system.length
        
        
        converted_stars = data.ParticlesWithUnitsConverted(
            stars, 
            convert_nbody.as_converter_from_si_to_nbody())
        
        self.assertEquals(stars[0].mass, 10 | nbody_system.mass)
        self.assertAlmostEqual(converted_stars[0].mass, 100.0 | units.kg, 5)
        
        converted_star = converted_stars[0]
        
        expected = units.m.new_quantity([50.0, 100.0, 150.0])
        for i in range(3):
            self.assertAlmostEqual(converted_star.position[i], expected[i], 6)
            
        converted_star.position = [100.0, 200.0, 300.0] | units.m
        star = stars[0]
        
        expected = nbody_system.length([20.0, 40.0, 60.0])
        for i in range(3):
            self.assertAlmostEqual(star.position[i], expected[i], 6)
            
        
        


class TestParticlesWithChildren(amusetest.TestCase):
    
    def test1(self):
        
        all = data.Particles(3)
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        
        self.assertEquals(child1.parent.key, parent.key)
        
        children = parent.children()
        
        self.assertEquals(len(children), 2)
        
    def test2(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = data.Particles(3)
        all.mass = [4.0, 3.0, 1.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        code1.particles.add_particles(parent.as_set())
        code2.particles.add_particles(parent.children())
        
        self.assertEquals(len(code1.particles), 1)
        self.assertEquals(len(code2.particles), 2)
        
        code1.legacy_interface.set_mass([0], [10000.0])
        code2.legacy_interface.set_mass([0], [9000.0])
        
        self.assertEquals(parent.mass, 4.0 | units.kg)
        self.assertEquals(child1.mass, 3.0 | units.kg)
        
        code1.particles.copy_values_of_all_attributes_to(all)
        
        self.assertEquals(parent.mass, 10.0 | units.kg)
        self.assertEquals(child1.mass, 3.0 | units.kg)
        
        
        code2.particles.copy_values_of_all_attributes_to(all)
        
        self.assertEquals(parent.mass, 10.0 | units.kg)
        self.assertEquals(child1.mass, 9.0 | units.kg)
        
        
    
    def test3(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = data.Particles(5)
        all.mass = [4.0, 3.0, 1.0, 6.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        self.assertEquals(parent.parent, None)
        self.assertEquals(child1.parent, parent)
        self.assertEquals(child2.parent, parent)
        
        
        all_except_children = all.difference(parent.children())
        code1.particles.add_particles(all_except_children)
        code2.particles.add_particles(parent.children())
        
        self.assertEquals(len(code1.particles), 3)
        self.assertEquals(len(code2.particles), 2)
        
    def test4(self):
        all = data.Particles(5)
        all.mass = [1.0, 2.0, 3.0, 4.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        child3 = all[3]
        child4 = all[4]
        
        parent.add_child(child1)
        child1.add_child(child2)
        child2.add_child(child3)
        child3.add_child(child4)
        
        self.assertEquals(len(parent.children()), 1)
        self.assertEquals(len(parent.descendents()), 4)
        self.assertEquals(len(child1.descendents()), 3)
        self.assertEquals(len(child2.descendents()), 2)
        self.assertEquals(len(child3.descendents()), 1)
        
    def test5(self):
        all = data.Particles(5)
        all.mass = [1.0, 2.0, 3.0, 4.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        child3 = all[3]
        child4 = all[4]
        
        parent.add_child(child1)
        child1.add_child(child2)
        child2.add_child(child3)
        child3.add_child(child4)
        
        copy = all.copy()
        
        self.assertEquals(copy[0].parent, None)
        self.assertEquals(copy[1].parent, copy[0])
        self.assertEquals(copy[2].parent, copy[1])
        self.assertEquals(copy[3].parent, copy[2])
        self.assertEquals(copy[4].parent, copy[3])
        self.assertEquals(copy[1].parent, parent)
        self.assertEquals(len(copy[0].descendents()), 4)
        
class TestParticlesSuperset(amusetest.TestCase):
    
    def test1(self):
        print "Test1: getting attributes of a particle superset."
        superset = data.Particles(2)
        superset.x = [1.0, 2.0] | units.m
        set2 = data.Particles(2)
        set2.x = [4.0, 5.0] | units.m
        particle1 = data.Particle()
        particle1.x = 3.0 | units.m
        particle2 = data.Particle()
        particle2.x = 6.0 | units.m
        for x in [particle1, set2, particle2]:
            superset = data.ParticlesSuperset([superset, x.as_set()])
        self.assertTrue(isinstance(superset, data.ParticlesSuperset))
        self.assertEqual(len(superset),6)
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]|units.m))
        # Check whether it returns the right value for the right key when the key order is 'random'
        dictionary = dict(zip(superset.key, superset.x))
        sorted_keys = sorted(superset.get_all_keys_in_store())
        sorted_values = superset.get_values_in_store(sorted_keys,['x'])[0]
        for key, value in zip(sorted_keys, sorted_values):
            self.assertEqual(dictionary[key],value)
    
    def test2(self):
        print "Test2: setting attributes of a particle superset."
        superset = data.Particles(2)
        set2 = data.Particles(2)
        particle1 = data.Particle()
        particle2 = data.Particle()
        for x in [particle1, set2, particle2]:
            superset = data.ParticlesSuperset([superset, x.as_set()])
        self.assertTrue(isinstance(superset, data.ParticlesSuperset))
        self.assertEqual(len(superset),6)
        superset.x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] | units.m
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]|units.m))
        superset.y = 9.0 | units.m
        self.assertEqual(superset.y, ([9.0, 9.0, 9.0, 9.0, 9.0, 9.0]|units.m))
        superset.z = [-1.0, 1.0] | units.m
        self.assertEqual(superset.z, ([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]|units.m))
        # Check whether it sets the value of the right particle when the key order is 'random'
        sorted_keys = sorted(superset.get_all_keys_in_store())
        sorted_values = superset.get_values_in_store(sorted_keys,['x'])[0]
        superset.set_values_in_store(sorted_keys,['zz'],[sorted_values])
        dictionary = dict(zip(superset.key, superset.zz))
        for key, value in zip(sorted_keys, sorted_values):
            print dictionary[key],value
            self.assertEqual(dictionary[key],value)
            
    def test3(self):
        set1 = data.Particles(2)
        set2 = data.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        set2.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0, 5.0, 12.0, 21.0] | units.m ** 2))
        
    def test4(self):
        set1 = data.Particles(2)
        set2 = data.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o)
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0, 5.0, 12.0, 21.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0, 7.0, 14.0, 23.0] | units.m ** 2))
        
    def test5(self):
        set1 = data.Particles(2)
        set2 = data.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o)
        superset = data.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0] | units.m ** 2))
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0] | units.m ** 2))
        
    def test6(self):
        set1 = data.Particles(2)
        set2 = data.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        set2.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0] | units.m ** 2))
        superset = data.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0] | units.m ** 2))
    
    def test7(self):
        def xtimesy1(x,y):
            raise Exception("error querying function on empty set")
            
        set1 = data.Particles(2)
        set2 = data.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_function_attribute("xtimesy", lambda p : p.x * p.y)
        set2.add_function_attribute("xtimesy", xtimesy1)
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy(), ([3.0, 8.0] | units.m ** 2))
        superset = data.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy(), ([3.0, 8.0] | units.m ** 2))
        
    def test8(self):
        set1 = data.Particles(2)
        set2 = data.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o, lambda all, p, o : p.x * p.x - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o, lambda all, p, o : p.x * p.x - o)
        superset = data.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset[1].xtimesypluso(0.0 | units.m ** 2), (4.0 | units.m ** 2))
        self.assertEquals(list(superset)[1].xtimesypluso(0.0 | units.m ** 2), (4.0 | units.m ** 2))
        
    
    def test9(self):
        
        particles1 = data.Particles(2)
        particles1.mass = 10 | units.kg
        particles2 = data.Particles(2)
        particles2.mass = 20 | units.kg
        superset = data.ParticlesSuperset([particles1, particles2])
            
        self.assertTrue(hasattr(superset, 'mass'))
        self.assertFalse(hasattr(superset, 'radius'))
        particles1.radius = 10 | units.m
        self.assertFalse(hasattr(superset, 'radius'))
        particles2.radius = 20 | units.m
        self.assertTrue(hasattr(superset, 'radius'))
    
class TestSliceParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test: slice a particle set."
        number_of_particles = 10
        original_set = data.Particles(number_of_particles)
        self.assertTrue(isinstance(original_set, data.Particles))
        self.assertEqual(len(original_set),number_of_particles)
        print "Defining all kind of slices of the particle set..."
        subset1 = original_set[:2]  # contains first two particles
        subset2 = original_set[2:]  # contains the rest of the particles
        odd = original_set[1::2]    # contains all particles with odd indices
        even = original_set[::2]    # contains all particles with even indices
        reverse = original_set[::-1]# contains all particles in reverse order
        all = original_set[:]       # contains all particles
        one = original_set[3]       # contains one particle (Particle)
        another = original_set[5:6] # contains one particle (ParticlesSubset)
        empty = original_set[2:2]   # contains no particle
        print "Checking their type (slicing returns a subset, indexing returns a particle)... ",
        self.assertTrue(isinstance(subset1, data.ParticlesSubset))
        self.assertTrue(isinstance(subset2, data.ParticlesSubset))
        self.assertTrue(isinstance(odd,     data.ParticlesSubset))
        self.assertTrue(isinstance(even,    data.ParticlesSubset))
        self.assertTrue(isinstance(reverse, data.ParticlesSubset))
        self.assertTrue(isinstance(all,     data.ParticlesSubset))
        self.assertTrue(isinstance(one,     data.Particle))
        self.assertTrue(isinstance(another, data.ParticlesSubset))
        self.assertTrue(isinstance(empty,   data.ParticlesSubset))
        print "ok!"
        print "Checking their length... ",
        self.assertEqual(len(subset1), 2)
        self.assertEqual(len(subset2), number_of_particles-2)
        self.assertEqual(len(odd),     int((number_of_particles/2.0)))
        self.assertEqual(len(even),    int((0.5+number_of_particles/2.0)))
        self.assertEqual(len(reverse), number_of_particles)
        self.assertEqual(len(all),     number_of_particles)
        self.assertEqual(len([one]),   1)
        self.assertEqual(len(another), 1)
        self.assertEqual(len(empty),   0)
        print "ok!"
    
class TestAddParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test1: create a particle subset by adding a particle to a set."
        original_set = data.Particles(4)
        original_set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        set = original_set[:2]
        particle = original_set[3]
        self.assertTrue(isinstance(set, data.ParticlesSubset))
        self.assertTrue(isinstance(particle, data.Particle))
        
        new_set = set + particle
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)+1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        
        set += particle
        self.assertTrue(isinstance(set, data.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print "Test2: create a particle subset by adding a set to a set."
        original_set = data.Particles(5)
        original_set.x = [1.0, 2.0, -789.0, 3.0, 4.0] | units.m
        set1 = original_set[:2]
        set2 = original_set[3:]
        self.assertTrue(isinstance(set1, data.ParticlesSubset))
        self.assertTrue(isinstance(set2, data.ParticlesSubset))
        
        new_set = set1 + set2
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)+len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        
        set1 += set2
        self.assertTrue(isinstance(set1, data.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print "Test3: create a particle superset by adding a particle to a set."
        set = data.Particles(2)
        set.x = [1.0, 2.0] | units.m
        particle = data.Particle()
        particle.x = 3.0 | units.m
        
        superset = data.ParticlesSuperset([set, particle.as_set()])
        self.assertTrue(isinstance(superset, data.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+1)
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0]|units.m))
        
        set2 = data.Particles(2)
        set2.x = [3.0, 4.0] | units.m
        superset = data.ParticlesSuperset([set, set2])
        self.assertTrue(isinstance(superset, data.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+len(set2))
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test4(self):
        print "Test4: check if the particle is already part of the set."
        set = data.Particles(2)
        particle = data.Particle()
        set = data.ParticlesSuperset([set, particle.as_set()])
        self.assertRaises(AmuseException, data.ParticlesSuperset, [set, particle.as_set()], 
            expected_message = "Unable to add a particle, because it was already part of this set.")
        self.assertEqual(len(set),3)
        other_set = data.Particles(2)
        other_set = data.ParticlesSuperset([other_set, particle.as_set()])
        self.assertRaises(AmuseException, data.ParticlesSuperset, [set, other_set], 
            expected_message = "Unable to add a particle, because it was already part of this set.")
    
    def test5(self):
        print "Test5: recursive addition, create a new superset from supersets."
        particle = data.Particle()
        set1 = data.Particles(2)
        set2 = data.Particles(2)
        set3 = data.Particles(2)
        set4 = data.Particles(2)
        superset1 = data.ParticlesSuperset([set1, set2])
        superset2 = data.ParticlesSuperset([set3, set4])
        for x in [particle, set3, superset2]:
            supersuperset = data.ParticlesSuperset([superset1, x.as_set()])
            self.assertTrue(isinstance(supersuperset, data.ParticlesSuperset))
            self.assertEqual(len(supersuperset),len(superset1)+len(x.as_set()))
            supersuperset.mass = 1 | units.kg
            self.assertEqual(x.mass, 1 | units.kg)
    
    def test6(self):
        print "Test6: check if the particle belongs to the same particle set as self."
        set1 = data.Particles(2)
        set2 = data.Particles(2)
        particle = set2[0]
        self.assertRaises(AmuseException, lambda: set1 + set2, 
            expected_message = "Can't create new subset from particles belonging to "
            "separate particle sets. Try creating a superset instead.")
        self.assertRaises(AmuseException, lambda: set1 + particle, 
            expected_message = "Can't create new subset from particles belonging to "
            "separate particle sets. Try creating a superset instead.")
    
    def test7(self):
        print "Test7: add a particle (set) to a particle."
        original_set = data.Particles(4)
        particle1 = original_set[0]
        particle2 = original_set[1]
        set = original_set[2:]
        self.assertTrue(isinstance(particle1, data.Particle))
        self.assertTrue(isinstance(particle2, data.Particle))
        self.assertTrue(isinstance(set, data.ParticlesSubset))
        new_set = particle1 + particle2
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),2)
        new_set = particle1 + set
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),3)
    
class TestSubtractParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test1: create a particle subset by removing a particle from a set."
        set = data.Particles(4)
        set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        particle = set[2]
        self.assertTrue(isinstance(particle, data.Particle))
        
        new_set = set - particle
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)-1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        
        set -= particle
        self.assertTrue(isinstance(set, data.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print "Test2: create a particle subset by removing a set from a set."
        set1 = data.Particles(5)
        set1.x = [1.0, 2.0, -789.0, 3.0, 4.0, -456.0] | units.m
        set2 = set1[2::3]
        self.assertTrue(isinstance(set1, data.Particles))
        self.assertTrue(isinstance(set2, data.ParticlesSubset))
        
        new_set = set1 - set2
        self.assertTrue(isinstance(new_set, data.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)-len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        
        set1 -= set2
        self.assertTrue(isinstance(set1, data.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print "Test3: check if the particle is actually part of the set."
        set = data.Particles(2)
        particle = data.Particle()
        self.assertRaises(AmuseException, lambda: set - particle, 
            expected_message = "Unable to subtract a particle, because it is not part of this set.")
    
    def test4(self):
        print "Test4: recursive subtraction, remove particles until the set is empty."
        set = data.Particles(10)
        self.assertEqual(len(set), 10)
        while len(set):
            set -= set[0]
        self.assertEqual(len(set), 0)
    
    def test5(self):
        print "Test5: check if it's possible to subtract particle(s) from a particle."
        particle = data.Particle()
        self.assertRaises(AmuseException, lambda: particle - particle, 
            expected_message = "Cannot subtract particle(s) from a particle.")
        particle2 = data.Particle()
        self.assertRaises(AmuseException, lambda: particle - particle2, 
            expected_message = "Cannot subtract particle(s) from a particle.")
    

class TestIterateOverParticles(amusetest.TestCase):
    
    def iterate_over_array(self, particles):
        for x in particles:
            x.radius
    
    def iterate_over_particles(self, particles):
        for x in particles:
            x.key


    def test1(self):
        self.total_number_of_points = 10000
    
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        particles = data.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles1(particles)
        t1 = time.time()
        dt0 = t1 - t0
        
        particles = data.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles2(particles)
        t1 = time.time()
        dt1 = t1 - t0
        
        print dt1, dt0, dt1 / dt0
    
        self.assertTrue((dt1 / dt0) < 20)

    def iterate_over_particles1(self, particles):
        for x in particles:
            x.key
    
    def iterate_over_particles2(self, particles):
        for x in particles:
            x.radius
    
    def test2(self):
        self.total_number_of_points = 10000
    
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        
        particles = data.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        #self.iterate_over_array(particles)
        for key in particles.get_all_keys_in_store():
            key
        t1 = time.time()
        dt0 = t1 - t0
        
        particles = data.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles2(particles)
        t1 = time.time()
        dt1 = t1 - t0
        
        print dt1, dt0, dt1 / dt0
    
        self.assertTrue((dt1 / dt0) < 400)
