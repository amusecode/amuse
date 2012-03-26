import numpy
import time
import sys
import pickle

from amuse.test import amusetest
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support.exceptions import AmuseException
from amuse.support.interface import InCodeComponentImplementation

from amuse import datamodel

class TestParticlesProperties(amusetest.TestCase):
    
    def test1(self):
        
        particles = datamodel.Particles(2)
        particles.mass = 10 | units.kg
        
        self.assertTrue(hasattr(particles, 'collection_attributes'))
        
        particles.collection_attributes.timestamp = 1 | units.yr
        self.assertEquals(particles.collection_attributes.timestamp,  1 | units.yr)
        
        particles.collection_attributes.a  = 2
        self.assertEquals(particles.collection_attributes.a,  2)
    
    
    def test2(self):
        
        particles = datamodel.Particles(2)
        particles.collection_attributes.timestamp = 1 | units.yr
        
        self.assertEquals(str(particles.collection_attributes), "timestamp: 1 yr")
        
        particles.collection_attributes.a = 2
        self.assertEquals(str(particles.collection_attributes), "timestamp: 1 yr\na: 2")
        
    def test3(self):
        
        particles1 = datamodel.Particles(2)
        particles1.collection_attributes.timestamp = 1 | units.yr
        particles1.collection_attributes.a = 2
        particles2 = particles1.copy()
        
        self.assertEquals(particles2.collection_attributes.timestamp,  1 | units.yr)
        self.assertEquals(particles2.collection_attributes.a,  2)
        self.assertEquals(str(particles2.collection_attributes), "timestamp: 1 yr\na: 2")
        
    def test4(self):
        
        particles1 = datamodel.Particles(2)
        particles1.collection_attributes.timestamp = 1 | units.yr
        particles1.collection_attributes.a = 2
        pickled_string = pickle.dumps(particles1)
        particles2 = pickle.loads(pickled_string)
        self.assertEquals(particles2.collection_attributes.timestamp,  1 | units.yr)
        self.assertEquals(particles2.collection_attributes.a,  2)
        
        
        
