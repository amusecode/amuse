import unittest

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.data import binding
from amuse.support.data import attributes

import numpy

class TestBase(unittest.TestCase):
    pass
    
class TestParticles(TestBase):
    
    def test1(self):
        
        particles = core.Particles(2)
        for i in range(3):
            particles.mass = (i * 1.0) | units.kg
            particles.savepoint((i + 1) * 1.0 | units.s)
        
        particles = particles.previous_state()
        
        masses = particles[0].get_timeline_of_attribute("mass")
        self.assertEquals(len(masses), 3)
        self.assertEquals(masses[0][0], 1.0 | units.s)
        self.assertEquals(masses[1][0], 2.0 | units.s)
        self.assertEquals(masses[2][0], 3.0 | units.s)
        
class TestStars(TestBase):

    def test1(self):
        stars = core.Stars(2)
        stars[0].mass = 10 | units.g
        stars[0].position = units.m(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].position = units.m(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m, stars.center_of_mass().x)

    def test2(self):
        stars = core.Stars(2)
        stars[0].mass = 10 | units.g
        stars[0].velocity = (units.m / units.s)(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].velocity = (units.m / units.s)(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m / units.s, stars.center_of_mass_velocity().x)
        self.assertEquals(1.0 | units.m / units.s, stars.center_of_mass_velocity().y)
        
    
    def test3(self):
        stars = core.Stars(2)
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
        stars = core.Stars(2)
        stars[0].x = 1.0  | units.km
        stars[0].y = 2000.0 | units.m
        stars[0].z = 3500.0 | units.m
        
        self.assertEquals(stars.position[0], [1000.0, 2000.0, 3500.0] | units.m)    



class TestParticlesWithBinding(TestBase):
    class TestInterface(object):
        attribute_definitions = [
            attributes.AttributeDefinition(
                name = "mass",
                setup_parameters = ["mass"],
                setter = ("set_mass", ["mass"]),
                getter = ("get_mass", ["mass"]),
                description = "mass of a star",
                unit = units.g,
                default = 1 | units.g         
            ),
        ]
        class InCodeAttributeStorage(binding.InCodeAttributeStorage):
            new_particle_method = binding.NewParticleMethod(
                "new_particle", 
                (
                    ("mass", "mass", units.g),
                )
            )
            
            getters = (
                binding.ParticleGetAttributesMethod(
                    "get_mass",
                    (
                        ("mass", "mass", units.g),
                    )
                ),
            )
            
            
            setters = (
                binding.ParticleGetAttributesMethod(
                    "set_mass",
                    (
                        ("mass", "mass", units.g),
                    )
                ),
            )
            
        def __init__(self):
            self.particles = core.Particles()
            self.particles._private.attribute_storage = self.InCodeAttributeStorage(self)
            self.masses = {}
            
        def get_mass(self, id):
            masses = []
            errors = []
            for x in id:
                masses.append(self.masses[x])
                errors.append(0)
            return {
                "mass": masses,
                "__result": errors,
            }
        
        def set_mass(self, id, mass):
            for i,m in zip(id,mass):
                self.masses[i] = m
                
            return {
                "__result": [0] * len(id),
            }
            
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
        
        
    
    def test1(self):
        interface = self.TestInterface()
        interface.particles._set_particles(
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
        interface.particles._set_particles(
            [1, 2],
            ["mass"],
            [[3.0, 4.0] | units.kg]
        )
        
        remote_particles = interface.particles
        
        self.assertEquals(len(remote_particles), 2)
        
        interface.particles._set_particles(
            [3, 4],
            ["mass"],
            [[5.0, 6.0] | units.kg]
        )
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        interface.particles._remove_particles((1,3))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
    
    def test3(self):
        interface = self.TestInterface()
        
        local_particles = core.Stars(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        local_particles2 = core.Stars(2)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0])
       
        remote_particles.add_particles(local_particles2)
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        keys = remote_particles._get_keys()
        interface.particles._remove_particles((keys[0], keys[2]))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
        
    
    def test4(self):
        interface = self.TestInterface()
        
        local_particles = core.Stars(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        particle = core.Particle()
        particle.mass = units.g.new_quantity(10.0)
        local_particles.add_particle(particle)
        
        self.assertEquals(len(local_particles), 3)
        
    
    
    def test5(self):
        interface = self.TestInterface()
        
        local_particles1 = core.Stars(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particles2 = core.Stars(3)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0, 7.0])
        
        local_particles1.add_particles(local_particles2)
        
        self.assertEquals(len(local_particles1), 5)
        self.assertEquals(len(local_particles2), 3)
        
        local_particles1.synchronize_to(remote_particles)
        
        local_particles1._remove_particles([local_particles1._get_keys()[0]])
        local_particles1.synchronize_to(remote_particles)
        self.assertEquals(len(remote_particles), 4)
        
    
    def test6(self):
        interface = self.TestInterface()
        
        local_particles1 = core.Stars(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particle = core.Particle()
        local_particle.mass = units.kg.new_quantity(5.0)
        
        local_particles1.add_particle(local_particle)
        
        self.assertEquals(len(local_particles1), 3)
        self.assertEquals(len(remote_particles), 2)
        
        local_particles1.synchronize_to(remote_particles)
        
        self.assertEquals(len(remote_particles), 3)
        
        remote_particles.remove_particle(local_particle)
        self.assertEquals(len(remote_particles), 2)
        
        
        
        
        
        
        
class TestParticlesWithUnitsConverted(TestBase):
    
    def test1(self):
        stars = core.Stars(2)
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
        
        converted_stars = core.ParticlesWithUnitsConverted(stars, LengthMassConverter())
        
        self.assertEquals(stars[0].mass, 10 | units.g)
        self.assertEquals(converted_stars[0].mass, 10 | units.m)
        
        
        converted_stars[0].mass = 30 | units.m
        
        self.assertEquals(stars[0].mass, 30 | units.g)
        
    
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10 | units.kg , 5 | units.m )
        
        stars = core.Stars(1)
        stars[0].mass = 10 | nbody_system.mass
        stars[0].x = 10.0 | nbody_system.length
        stars[0].y = 20.0 | nbody_system.length
        stars[0].z = 30.0 | nbody_system.length
        
        
        converted_stars = core.ParticlesWithUnitsConverted(
            stars, 
            convert_nbody.as_converter_from_si_to_nbody())
        
        self.assertEquals(stars[0].mass, 10 | nbody_system.mass)
        print converted_stars[0].mass.number 
        print (100.0 | units.kg).number 
        self.assertAlmostEqual(converted_stars[0].mass.value_in(units.kg), 100.0, 5)
        
        converted_star = converted_stars[0]
        
        expected = [50.0, 100.0, 150.0]
        for i in range(3):
            self.assertAlmostEqual(converted_star.position[i].value_in(units.m), expected[i], 6)
            
        converted_star.position = [100.0, 200.0, 300.0] | units.m
        star = stars[0]
        
        expected = [20.0, 40.0, 60.0]
        for i in range(3):
            self.assertAlmostEqual(star.position[i].value_in(nbody_system.length), expected[i], 6)
            
        
        


class TestParticlesWithChildren(TestBase):
    
    def test1(self):
        
        all = core.Particles(3)
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        
        self.assertEquals(parent.parent_id, 0 | units.none)
        
        children = parent.children()
        
        self.assertEquals(len(children), 2)
        
    def test2(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = core.Particles(3)
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
        
        code1.set_mass([0], [10000.0])
        code2.set_mass([0], [9000.0])
        
        self.assertEquals(parent.mass , 4.0 | units.kg)
        self.assertEquals(child1.mass , 3.0 | units.kg)
        
        code1.update_particles(all)
        
        self.assertEquals(parent.mass , 10.0 | units.kg)
        self.assertEquals(child1.mass , 3.0 | units.kg)
        
        
        code2.update_particles(all)
        
        self.assertEquals(parent.mass , 10.0 | units.kg)
        self.assertEquals(child1.mass , 9.0 | units.kg)
        
        
        code1.update_particles(all)
        code2.update_particles(all)
        
    
    def test3(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = core.Particles(5)
        all.mass = [4.0, 3.0, 1.0, 6.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        all_except_children = all.difference(parent.children())
        code1.particles.add_particles(all_except_children)
        code2.particles.add_particles(parent.children())
        
        self.assertEquals(len(code1.particles), 3)
        self.assertEquals(len(code2.particles), 2)
        
    def test4(self):
        all = core.Particles(5)
        all.mass = [1.0, 2.0, 3.0, 4.0 , 5.0] | units.kg
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
        
        
        
        
        
        
        
        
        
