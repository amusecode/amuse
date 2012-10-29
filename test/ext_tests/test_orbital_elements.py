import sys
import os
import numpy

from amuse.test import amusetest


from amuse.ext import orbital_elements

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel



class KeplerTests(amusetest.TestCase):

    def test1(self):
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length
        )
        
        self.assertEquals(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [1,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,numpy.sqrt(2),0] | nbody_system.speed)
    
    
    def test2(self):
        #test going around in a circular orbit
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 90,
        )
        
        self.assertEquals(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,1,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [-numpy.sqrt(2),0,0] | nbody_system.speed)
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 180,
        )
        
        self.assertEquals(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [-1,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,-numpy.sqrt(2),0] | nbody_system.speed)
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 270,
        )
        
        self.assertEquals(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,-1,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [numpy.sqrt(2),0,0] | nbody_system.speed)
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 45,
        )
        
        self.assertEquals(len(binary), 2)
        print binary[1]
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0.5 * numpy.sqrt(2),0.5 * numpy.sqrt(2),0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [-1,1,0] | nbody_system.speed)
        
    def xtest3(self):
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = orbital_elements.new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0.5
        )
        
        self.assertEquals(len(binary), 2)
        print binary[1], numpy.sqrt(2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,numpy.sqrt(2),0] | nbody_system.speed)
    
