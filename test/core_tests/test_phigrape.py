from amuse.test import amusetest
import os.path
import numpy
from amuse.io import phigrape
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
from amuse import datamodel
class Test(amusetest.TestCase):


    def test(self):
        directory = os.path.dirname(__file__)
        instance = phigrape.Inp2Particles()
        instance.convert_to_particles(os.path.join(directory, 'plummer_100.ini'))
        rev_instance = phigrape.Particles2Inp()
        rev_instance.convert_to_inp(instance.Particles,os.path.join(directory, 'plummer_back_100.ini'))

        control_instance = phigrape.Inp2Particles()
        control_instance.convert_to_particles(os.path.join(directory, 'plummer_back_100.ini'))

        self.assertAlmostEqual(control_instance.Particles.mass, instance.Particles.mass, 16)
        self.assertAlmostEqual(control_instance.Particles[1].position, instance.Particles[1].position, 16)
        self.assertAlmostEqual(control_instance.Particles[1].velocity, instance.Particles[1].velocity, 16)
