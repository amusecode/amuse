
import sys

from amuse.legacy.bhtree import muse_dynamics_mpi as mpi_interface
from amuse.legacy.support import core as legacy_core

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

import numpy

from legacy_support import TestWithMPI

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
    

class TestMPIInterface(TestWithMPI):
    
    def setUp(self):
        super(TestMPIInterface, self).setUp()
        nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        
    def test1(self):
        instance = mpi_interface.BHTreeInterface()
        instance.setup_module()
        instance.add_particle(1, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(instance.get_number(), 1)
        instance.cleanup_module()
        del instance
        
    def test2(self):
        instance = mpi_interface.BHTreeInterface(debug_with_gdb=False)
        instance.eps2 = 0.101
        self.assertEquals(0.101, instance.eps2)
        instance.eps2 = 0.110
        self.assertEquals(0.110, instance.eps2)
        del instance
        
    def test3(self):
        instance = mpi_interface.BHTreeInterface()
        instance.flag_collision = 1
        self.assertEquals(1, instance.flag_collision)
        instance.flag_collision = 0
        self.assertEquals(0, instance.flag_collision)
        del instance
        
    def test4(self):
        class BHTree2(mpi_interface.BHTreeInterface):
            channel_factory = legacy_core.MultiprocessingMPIChannel
            pass
        
        instance = BHTree2()
        instance.setup_module()
        instance.add_particle(1, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(instance.get_number(), 1)
        instance.cleanup_module()
        del instance

    def test5(self):
        import socket
        hostname = socket.gethostname()
        
        instance = mpi_interface.BHTreeInterface(hostname = hostname)
        instance.setup_module()
        instance.cleanup_module()
        
        del instance
        
