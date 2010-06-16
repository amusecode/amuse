from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.interface import CodeInterface

import numpy
import inspect
import collections

def unpack_slice(slice):
    start = 0 if slice.start is None else slice.start
    stop = None if slice.stop is None else slice.stop
    step = 1 if slice.step is None else slice.step
    return start, stop, step

def combine_slices(slice0, slice1):
    start0,stop0,step0 = unpack_slice(slice0)
    start1,stop1,step1 = unpack_slice(slice1)
            
    newstart = start0 + (start1 * step0)
    
    if stop0 is None:
        if stop1 is None:
            newstop = None
        else:
            newstop = (stop1 * step0) + start0
    elif stop1 is None:
        newstop = stop0
    else:
        newstop = min((stop1 * step0) + start0, stop0)
    
    newstep = step0 * step1
    return newstart, newstop, newstep

def combine_indices(index0, index1):
    if isinstance(index0, tuple):
        if len(index0) == 1:
           index0 = index0[0]
        else:
            result = []
            result.extend(index0[:-1])
            continuation = combine_indices(index0[-1], index1)
            if isinstance(continuation, collections.Sequence):
                result.extend(continuation)
            else: 
                result.append(continuation)
            return tuple(result)
    
    if isinstance(index0, int) or isinstance(index0, long):
        return (index0, index1)
    elif isinstance(index0, slice):
        if isinstance(index1, int) or isinstance(index1, long):
            start,stop,step = unpack_slice(index0)
            return start + (index1 * step)
        else:
            start,stop,step = combine_slices(index0, index1)
            return numpy.s_[start:stop:step]
    else:
        raise Exception("index must be a integer, slice or sequence")
    
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
        self.assertTrue(a[combine_indices(numpy.s_[0:2],0)].shape, a[0:2][0].shape)
        self.assertTrue(numpy.all(a[combine_indices(numpy.s_[0:2],0)] ==  a[0:2][0]))

    def test4(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1][:]
        indirect = a[combine_indices(1,numpy.s_[:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test5(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[0:2][:]
        indirect = a[combine_indices(numpy.s_[0:2],numpy.s_[:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test6(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1:3][1:]
        indirect = a[combine_indices(numpy.s_[1:3],numpy.s_[1:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test7(self):
        a = numpy.arange(30).reshape(6,5)
        direct =  a[1:5:2][1:]
        indirect = a[combine_indices(numpy.s_[1:5:2],numpy.s_[1:])]
        print a
        print "direct", direct
        print "indirect", indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test8(self):
        a = numpy.arange(30)
        direct =  a[2:14:3][1:5:2]
        indirect = a[combine_indices(numpy.s_[2:14:3],numpy.s_[1:5:2])]
        print a
        print "direct", direct
        print "indirect", indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    
    def test9(self):
        a = numpy.arange(100)
        for s in range(0,40):
            for e in range(40,101):
                for step in range(1,5):
                    direct =  a[s:e:step][1:5:2]
                    indirect = a[combine_indices(numpy.s_[s:e:step],numpy.s_[1:5:2])]
                    print direct, indirect
                    self.assertEquals(indirect.shape, direct.shape)
                    self.assertTrue(numpy.all(indirect ==  direct))
        #self.assertTrue(False)
    
    
    def test10(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3][2][1]
        indirect = a[combine_indices(combine_indices(3,2),1)]
        print direct, indirect, combine_indices(combine_indices(3,2),1), a[(3,2,1)]
        self.assertEquals(indirect, direct)