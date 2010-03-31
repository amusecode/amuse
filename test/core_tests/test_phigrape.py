from amuse.support.data import core
from amuse.support.io import phigrape
from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import values


#import unittest
from amuse.test import amusetest
import os.path
import numpy

class Test(amusetest.TestCase):

    def xtest1(self):
        directory = os.path.dirname(__file__)
        convert_nbody = nbody_system.nbody_to_si(1|units.kg, 1|units.m)
        instance = phigrape.Inp2Particles()
        instance.convert_to_particles(os.path.join(directory, 'plummer_100.ini'), convert_nbody)
        #16k = 2**14
        self.assertAlmostEqual(instance.Particles.mass,
                               values.new_quantity(1.0/2**14*numpy.ones(101), units.kg),
                               8)

    def test2(self):
        directory = os.path.dirname(__file__)
        instance = phigrape.Inp2Particles()
        instance.convert_to_particles(os.path.join(directory, 'plummer_100.ini'))
        rev_instance = phigrape.Particles2Inp()
        rev_instance.convert_to_inp(I.Particles, 'plummer_back_100.ini')

        control_instance = phigrape.Inp2Particles()
        control_instance_.convert_to_particles(os.path.join(directory, 'plummer_back_100.ini'))

        self.assertAlmostEqual(control_instance.Particles.mass, instance.Particles.mass, 16)
        self.assertAlmostEqual(control_instance.Particles[1].position, instance.Particles[1].position, 16)
        self.assertAlmostEqual(control_instance.Particles[1].velocity, instance.Particles[1].velocity, 16)
