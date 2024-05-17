from amuse import datamodel
from amuse.rfi.core import *
from amuse.support import exceptions
from amuse.support.interface import InCodeComponentImplementation
from amusetest import TestWithMPI
from compile_tests.grid_implementation.interface import ForTestingInterface
from amuse.units import units

import os
from pathlib import Path


class ForTesting(InCodeComponentImplementation):

    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)

    def define_grids(self, object):
        object.define_grid('grid0', grid_class=datamodel.RectilinearGrid, state_guard="before_new_set_instance")
        object.set_grid_range('grid0', 'get_grid0_range')
        object.add_getter('grid0', 'get0',  names=['a'])
        object.add_setter('grid0', 'set0',  names=['a'])


class TestCImplementationInterface(TestWithMPI):

    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'grid_worker')

    def test1(self):
        instance = ForTestingInterface(self.exefile)
        error = instance.set0(1)
        a_out, error = instance.get0()
        instance.stop()
        self.assertEqual(a_out, 1)
        self.assertEqual(error, 0)

    def test2(self):
        instance = ForTesting(self.exefile)
        print(instance.grid0)
        instance.grid0.a = 12. | units.m
        self.assertEqual(instance.grid0.a, 12. | units.m)
        instance.stop()

    def test3(self):
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance1.grid0.a = 12. | units.m
        instance1.grid0.new_channel_to(instance2.grid0).copy_all_attributes()
        self.assertEqual(instance2.grid0.a, 12. | units.m)
        instance1.stop()
        instance2.stop()
