from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.interface import CodeInterface

from amuse.support.data.indexing import *

import numpy
import inspect
import collections

class TestGrids(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        print grid.mass
        self.assertEquals(grid.mass[0][1][2], 2.0 | units.kg)
    
    def test2(self):
        grid = core.Grid(5,4,3)
        grid.mass = units.kg.new_quantity(numpy.arange(5*4*3).reshape(5,4,3))
        self.assertEquals(grid.number_of_dimensions(), 3)
        subgrid = grid[1]
        self.assertEquals(subgrid.number_of_dimensions(), 2)
        print subgrid.mass
        self.assertEquals(subgrid.mass.number.shape, (4,3))
    
    def test3(self):
        grid = core.Grid(5,4,3)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        self.assertEquals(grid.number_of_dimensions(), 3)
       
        subgrid = grid[1][2]
        print subgrid.mass, type(subgrid._private.indices), subgrid._private.indices
        self.assertEquals(subgrid.number_of_dimensions(), 1)
        self.assertEquals(subgrid.mass.number.shape, (3,))
        self.assertTrue(numpy.all(values[1][2] == subgrid.mass.value_in(units.kg)))
        
    def test4(self):
        grid = core.Grid(5,4,3)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        self.assertEquals(grid.number_of_dimensions(), 3)
       
        gridpoint = grid[1][2][1]
        print values[(1,2,1)]
        print gridpoint.mass, gridpoint.index
        self.assertEquals(gridpoint.mass, 19.0 | units.kg)
        gridpoint = grid[1][2][2]
        self.assertEquals(gridpoint.mass, 20.0 | units.kg)
        
    def test5(self):
        grid = core.Grid(5,4,3)
        grid.add_calculated_attribute("squared_mass", lambda m : m * m, ["mass",])
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        gridpoint = grid[1][2][1]
        self.assertEquals(gridpoint.mass, 19.0 | units.kg)
        self.assertEquals(gridpoint.squared_mass, (19.0 | units.kg) ** 2)
        subgrid = grid[1][2]
        self.assertTrue(numpy.all(subgrid.squared_mass == ([18.0, 19.0, 20.0] | units.kg) ** 2))
        
    def test6(self):
        grid = core.Grid(5,4,3)
        grid.add_function_attribute("sum_mass", lambda grid, x : grid.mass.sum() + x, lambda grid, gridpoint, x : gridpoint.mass + x)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        gridpoint = grid[1][2][1]
        self.assertEquals(gridpoint.mass, 19.0 | units.kg)
        self.assertEquals(gridpoint.sum_mass(2.0 | units.kg), (21.0 | units.kg) )
        subgrid = grid[1][2]
        self.assertTrue(numpy.all(subgrid.sum_mass(2 | units.kg) == (18.0 + 19.0 + 20.0 + 2.0 | units.kg)))
        
    def test7(self):
        grid = core.Grid(5,4,3)
        grid.add_vector_attribute("position", ["x","y","z"])
        x = numpy.arange(5*4*3).reshape(5,4,3)
        y = x + 100.0
        z = x + 200.0
        grid.x = units.m.new_quantity(x)
        grid.y = units.m.new_quantity(y)
        grid.z = units.m.new_quantity(z)
        gridpoint = grid[1][2][1]
        print gridpoint.position
        self.assertEquals(gridpoint.position[0], 19 | units.m)
        self.assertEquals(gridpoint.position[1], 119 | units.m)
        self.assertEquals(gridpoint.position[2], 219 | units.m)
        subgrid = grid[1][2]
        print subgrid.position
        self.assertEquals(subgrid.position[1][0], 19 | units.m)
        self.assertEquals(subgrid.position[1][1], 119 | units.m)
        self.assertEquals(subgrid.position[1][2], 219 | units.m)
        
class TestIndexing(amusetest.TestCase):
    def test1(self):
        self.assertEquals(2, number_of_dimensions_after_index(3, 1))
        self.assertEquals(3, number_of_dimensions_after_index(3, numpy.s_[0:3]))
        self.assertEquals(1, number_of_dimensions_after_index(3, combine_indices(3,2)))
        self.assertEquals(0, number_of_dimensions_after_index(3, combine_indices(combine_indices(3,2),1)))

    def test2(self):
        a = numpy.arange(12).reshape(3,4)
        print a, a[0][0:2]
        self.assertEquals(a[combine_indices(0,1)], a[0][1])
        self.assertEquals(a[combine_indices(1,0)], a[1][0])
        self.assertTrue(numpy.all(a[combine_indices(1,numpy.s_[0:2])] == a[1][0:2]))
        indirect = combine_indices(0,1)
        self.assertEquals(number_of_dimensions(a, indirect), 0)
        
        
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
        
    def test11(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3]
        indirect = a[combine_indices(3,None)]
        print direct, indirect, combine_indices(combine_indices(3,2),1), a[(3,2,1)]
        
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
