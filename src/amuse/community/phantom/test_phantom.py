from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.phantom.interface import PhantomInterface, Phantom

from amuse.datamodel import Particles
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic import plummer
from amuse.ic.plummer import new_plummer_model

from amuse.test.suite.codes_tests.gd_tests import (
    _TestGravitationalDynamicsInterface,
)
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestPhantomInterface(TestWithMPI):
    def gravity_code_interface(self):
        return PhantomInterface

    def reference_includes(self):
        return "Price"

    def starting_particle_index(self):
        return 1

    def test_initialise(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        instance.stop()

    def test_literature_reference(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        reference_string = self.reference_includes()
        self.assertTrue(
            reference_string in instance.all_literature_references_string()
        )
        instance.stop()

    def test_add_and_retrieve_particles(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        # Phantom won't work with fewer than 7 particles!
        n = 7
        values = [1.0 * i for i in range(1, n)]
        instance.new_particle(
            values,
            values,
            values,
            values,
            values,
            values,
            values,
        )
        error = instance.commit_particles()
        self.assertEqual(error, 0)
        retrieved_state = instance.get_state(self.starting_particle_index())
        self.assertEqual(1.0, retrieved_state['mass'])
        retrieved_state = instance.get_state(
            self.starting_particle_index()+n-2
        )
        instance.cleanup_code()
        # For any particle other than a sink, Phantom has one fixed mass!
        self.assertEqual(1.0,  retrieved_state['mass'])
        instance.stop()

    def test_add_and_retrieve_sph_particles(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()

        # Phantom won't work with fewer than 7 particles!
        n = 100
        values = [1.0 * i for i in range(1, n)]
        instance.new_sph_particle(
            values,
            values,
            values,
            values,
            values,
            values,
            values,
            values,
        )
        error = instance.commit_particles()
        self.assertEqual(error, 0)
        retrieved_state = instance.get_state(self.starting_particle_index())
        self.assertEqual(1.0, retrieved_state['mass'])
        retrieved_state = instance.get_state(
            self.starting_particle_index()+n-2
        )
        instance.cleanup_code()
        # For any particle other than a sink, Phantom has one fixed mass!
        self.assertEqual(1.0,  retrieved_state['mass'])
        instance.stop()

    def test_parameters(self):
        interface = self.gravity_code_interface()
        # instance = self.new_instance_of_an_optional_code(interface)
        instance = interface(redirection="none")
        instance.initialize_code()
        gamma, error = instance.get_gamma()
        self.assertEqual(0, error)
        self.assertEqual(1., gamma)
        ieos, error = instance.get_ieos()
        self.assertEqual(0, error)
        self.assertEqual(1, ieos)



class TestPhantom(TestWithMPI):
    def test_initialise(self):
        instance = Phantom()
        instance.stop()

    def test_add_gasparticles(self):
        n_particles = 10
        instance = Phantom()
        gas = Particles(n_particles)
        gas.mass = 1 | nbody_system.mass
        gas.x = numpy.arange(n_particles) | nbody_system.length
        gas.y = numpy.arange(n_particles) | nbody_system.length
        gas.z = 0 | nbody_system.length
        gas.velocity = [0, 0, 0] | nbody_system.speed
        gas.u = 0 | nbody_system.speed**2
        instance.gas_particles.add_particles(gas)
        self.assertEqual(10, len(instance.gas_particles))
        instance.stop()
