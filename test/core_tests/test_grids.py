from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.interface import CodeInterface

import numpy
import inspect
import collections

def combine_indices(index0, index1):
    if isinstance(index0, collections.Sequence):
        result = []
        result.extend(index0)
        result.append(index)
        return result
    elif isinstance(index0, int) or isinstance(index0, long):
        return (index0, index1)
    elif isinstance(index0, slice):
        start = index0.start
        stop = index0.stop
        step = index0.step
        if step is None:
            step = 1
        if isinstance(index1, int) or isinstance(index1, long):
            return start + (index1 * step)
        else:
            raise Exception("not handled yet")
    else:
        return (index0, index1)
    
class TestGrids(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        print grid.mass
        self.assertEquals(grid.mass[0][1][2], 2.0 | units.kg)
        
    def test2(self):
        a = numpy.arange(12).reshape(3,4)
        print a, a[0][0:2]
        self.assertEquals(a[combine_indices(0,1)], a[0][1])
        self.assertEquals(a[combine_indices(1,0)], a[1][0])
        self.assertTrue(numpy.all(a[combine_indices(1,numpy.s_[0:2])] == a[1][0:2]))
        
    def test3(self):
        a = numpy.arange(12).reshape(3,4)
        print a, a[0:2][0], a[combine_indices(numpy.s_[0:2],0)]
        self.assertTrue(numpy.all(a[combine_indices(numpy.s_[0:2],0)] ==  a[0:2][0]))
