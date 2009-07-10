import unittest
import sys

from amuse.legacy.hermite0 import muse_dynamics_mpi as mpi_interface
from amuse.legacy.hermite0 import muse_dynamics_interface as swig_interface

        

class TestMPIInterface(unittest.TestCase):
    def test1(self):
        hermite = mpi_interface.Hermite()
        hermite.setup_module()
        new_state = hermite.dynamics_state()
        new_state.id = 1
        new_state.mass = 11.0
        hermite.add_particle(new_state)
        retrieved_state = hermite.get_state(1)
        self.assertEquals(new_state.mass,  retrieved_state.mass)
        self.assertEquals(hermite.get_number(), 1)
        hermite.cleanup_module()
        del hermite
    def test2(self):
        hermite = mpi_interface.Hermite()
        hermite.eps2 = 0.101
        self.assertEquals(0.101, hermite.eps2)
    def test3(self):
        hermite = mpi_interface.Hermite()
        hermite.flag_collision = 1
        self.assertEquals(1, hermite.flag_collision )

class TestSwigInterface(unittest.TestCase):
    def test1(self):
        hermite = swig_interface
        hermite.setup_module()
        new_state = hermite.dynamics_state()
        new_state.id = 1
        new_state.mass = 11.0
        hermite.add_particle(new_state)
        retrieved_state = hermite.get_state(1)
        self.assertEquals(new_state.mass,  retrieved_state.mass)
        self.assertEquals(hermite.get_number(), 1)
        hermite.cleanup_module()
    def test2(self):
        hermite = mpi_interface.Hermite()
        hermite.eps2 = 0.101
        self.assertEquals(0.101, hermite.eps2)
    def test3(self):
        hermite = mpi_interface.Hermite()
        hermite.flag_collision = 1
        self.assertEquals(1, hermite.flag_collision )
