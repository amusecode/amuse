import unittest

from amuse.support.units import units
from amuse.support.data import core
from amuse.support.data import binding
from amuse.support.data import attributes

import numpy

class TestBase(unittest.TestCase):
    pass
        
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
    class TestInterface(binding.InterfaceWithObjectsBinding):
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
        class InCodeAttributeStorage(binding.InCodeAttributeStorage2):
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
            binding.InterfaceWithObjectsBinding.__init__(self)
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
        local_particles2 = core.Stars(2)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0])
        
        local_particles1.add_particles(local_particles2)
        
        self.assertEquals(len(local_particles1), 4)
        self.assertEquals(len(local_particles2), 2)
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        remote_particles.add_particles(local_particles2)
        
        self.assertEquals(len(remote_particles), 6)
        
        
        
        
        
        
        
        
        
        
        
