import sys
import os
import numpy

from amuse.test import amusetest


from amuse.ext.orbital_elements import new_binary_from_orbital_elements,orbital_elements_from_binary

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel

from numpy import random

class KeplerTests(amusetest.TestCase):

    def test1(self):
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length
        )
        
        self.assertEquals(len(binary), 2)
        
        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity
        
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [1,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,numpy.sqrt(2),0] | nbody_system.speed)
    
    
    def test2(self):
        #test going around in a circular orbit
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 90,
        )
        
        self.assertEquals(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,1,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [-numpy.sqrt(2),0,0] | nbody_system.speed)
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 180,
        )
        
        self.assertEquals(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [-1,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,-numpy.sqrt(2),0] | nbody_system.speed)
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 270,
        )
        
        self.assertEquals(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,-1,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [numpy.sqrt(2),0,0] | nbody_system.speed)
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 45,
        )
        
        self.assertEquals(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0.5 * numpy.sqrt(2),0.5 * numpy.sqrt(2),0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [-1,1,0] | nbody_system.speed)


    def xtest3(self):
        mass1 = 1 | nbody_system.mass 
        mass2 = 1 | nbody_system.mass
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0.5
        )
        
        self.assertEquals(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,numpy.sqrt(2),0] | nbody_system.speed)
    
    def test4(self):
        numpy.random.seed(3456789)
        N=100
        
        mass1=random.random(N) | nbody_system.mass 
        mass2=random.random(N) | nbody_system.mass
        semi_major_axis=(-numpy.log(random.random(N))) | nbody_system.length 
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180       

        for arg in zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                  longitude_of_the_ascending_node,argument_of_periapsis):
          arg_=orbital_elements_from_binary(new_binary_from_orbital_elements(*arg))
          for i,(copy,org) in enumerate(zip(arg_,arg)):
            self.assertAlmostEquals(copy,org)

    def test5(self):
        numpy.random.seed(4567893)
        N=100
        
        mass1=random.random(N) | units.MSun 
        mass2=random.random(N) | units.MSun
        semi_major_axis=(-numpy.log(random.random(N))) | units.AU 
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180       

        for arg in zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                  longitude_of_the_ascending_node,argument_of_periapsis):
          arg_=orbital_elements_from_binary(new_binary_from_orbital_elements(*arg,G=constants.G),G=constants.G)
          for i,(copy,org) in enumerate(zip(arg_,arg)):
            self.assertAlmostEquals(copy,org)
