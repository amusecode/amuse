import unittest

from amuse.support.units import units
from amuse.support.data import binding
from amuse.support.data import core
from amuse.support.data import attributes

import numpy

class TestBase(unittest.TestCase):
    
    class Test(object):
        
        test_property = binding.CodeProperty("get_property", units.m)
        test_property2 = binding.CodeProperty("property_value", units.m)
        
        def get_property(self):
            return (self.property_value, self.errorcode)
        
    
    def test1(self):
        o = self.Test()
        o.errorcode = 0
        o.property_value = 2.0
        
        self.assertEquals(o.test_property, 2.0 | units.m)
        try:
            o.errorcode = -1
            o.test_property
            self.fail("method returned an errorcode, so the property should raise an exception")
        except Exception, ex:
            self.assertEquals(str(ex), "calling 'get_property' to get the value for property 'test_property' resulted in an error (errorcode -1)")
            
        try:
            o.test_property = 4.0 | units.m
            self.fail("values of properties cannot be set")
        except Exception, ex:
            self.assertEquals(str(ex), "property 'test_property' of a 'Test' object is read-only, you cannot change it's value")
    
    def test2(self):
        o = self.Test()
        o.property_value = 3.0
        
        self.assertEquals(o.test_property2, 3.0 | units.m)
        
        

class TestParticlesWithBinding(TestBase):
    class TestInterface(binding.InterfaceWithObjectsBinding):
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
        
        def get_colliding_particles(self):
            result = {}
            result['c1'] = self.colliding_particle1
            result['c2'] = self.colliding_particle2  
            result['__result'] = 0
            return result
        
        colliding_particles_method =  binding.ParticleQueryMethod(
            "get_colliding_particles",
            ("c1","c2")
        )
        
        def colliding_particles(self):
            subset = self.colliding_particles_method._run(self, self.particles)
            return subset
            
            
    def test1(self):
        
        interface = self.TestInterface()
        
        local_particles = core.Stars(4)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0, 5.0, 6.0])
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        interface.colliding_particle1 = 0
        interface.colliding_particle2 = 2
        
        subset = interface.colliding_particles()
        self.assertEquals(len(subset), 2)
        self.assertEquals(subset[0].mass, units.kg.new_quantity(3.0))
        self.assertEquals(subset[1].mass, units.kg.new_quantity(5.0))
        
        
        
        
        
        
