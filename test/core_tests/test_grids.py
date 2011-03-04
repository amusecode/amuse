from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.interface import InCodeComponentImplementation

from amuse.support.data.indexing import *

from amuse.support.data.grids import *

import numpy
import inspect
import collections

class TestGrids(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        self.assertEquals(grid.mass[0][1][2], 2.0 | units.kg)
        self.assertEquals(grid[0][1][2].mass, 2.0 | units.kg)
        self.assertEquals(len(grid.mass), 5)
        
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
        
    
    def test8(self):
        grid0 = core.Grid(5,4,3)
        x = numpy.arange(5*4*3).reshape(5,4,3)
        y = x + 100.0
        grid0.x = units.m.new_quantity(x)
        grid0.y = units.m.new_quantity(y)
        
        grid1 = core.Grid(5,4,3)
        x = numpy.arange(5*4*3).reshape(5,4,3)
        x = x + 200.0
        y = x + 200.0
        grid1.x = units.m.new_quantity(x)
        grid1.y = units.m.new_quantity(y)
    
        self.assertTrue(numpy.all(grid0[1][2].x != grid1[1][2].x))
        
        channel = grid0.new_channel_to(grid1)
        channel.copy_attributes(["x",])
        
        self.assertTrue(numpy.all(grid0[1][2].x == grid1[1][2].x))
        self.assertTrue(numpy.all(grid0[1][2].y != grid1[1][2].y))
        

        

    def test9(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        print grid.x
        self.assertEquals(grid[0][0][0].x, 0.1 | units.m)
        self.assertEquals(grid[0][0][0].y, 0.125 | units.m)
        self.assertEquals(grid[0][0][0].z, 0.25 | units.m)
        self.assertEquals(grid[...,0,0].x, [0.1,0.3,0.5,0.7,0.9] | units.m)
        self.assertEquals(grid[0,0,...].z, [0.25, 0.75] | units.m)
        
        cellsize = grid.cellsize()
        self.assertAlmostRelativeEquals(cellsize[0], 0.2 | units.m)
        self.assertAlmostRelativeEquals(cellsize[1], 0.25 | units.m)
        self.assertAlmostRelativeEquals(cellsize[2], 0.5 | units.m)
    
    

    def test11(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        iarray,jarray,karray = grid.indices()
        for i in range(5):
            for j in range(4):
                for k in range(2):
                    self.assertEquals(iarray[i][j][k], i)
                    self.assertEquals(jarray[i][j][k], j)
                    self.assertEquals(karray[i][j][k], k)
        print grid[0].indices()
        iarray,jarray,karray = grid.indices()
        i = 0
        for j in range(4):
            for k in range(2):
                self.assertEquals(iarray[i][j][k], i)
                self.assertEquals(jarray[i][j][k], j)
                self.assertEquals(karray[i][j][k], k)
       
        iarray,jarray,karray = grid[...,0,0].indices()
        j = 0
        k = 0
        print iarray
        for i in range(5):
            self.assertEquals(iarray[i], i)
            self.assertEquals(jarray[i], j)
            self.assertEquals(karray[i], k)
        print grid[0,0,...].indices()
        iarray,jarray,karray = grid[3,2,...].indices()
        i = 3
        j = 2
        for k in range(2):
            self.assertEquals(iarray[k], i)
            self.assertEquals(jarray[k], j)
            self.assertEquals(karray[k], k)
        iarray,jarray,karray = grid[2,...,1].indices()
        i = 2
        k = 1
        for j in range(4):
            self.assertEquals(iarray[j], i)
            self.assertEquals(jarray[j], j)
            self.assertEquals(karray[j], k)
    
    
    def test12(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        print grid.indices()[0]
        print grid[grid.x>0.4| units.m].indices()
        print grid.x[grid.x>0.4| units.m]
        self.assertEquals(grid[grid.x > 0.4| units.m].x.shape, (24,))
        self.assertEquals(grid[grid.x>0.4| units.m].x, grid.x[grid.x>0.4|units.m])
        iarray,jarray,karray = grid.indices()
        self.assertEquals(grid[grid.x>0.4| units.m].indices()[0],  iarray[grid.x>0.4| units.m])
        self.assertEquals(grid[grid.x>0.4| units.m].indices()[1],  jarray[grid.x>0.4| units.m])
        self.assertEquals(grid[grid.x>0.4| units.m].indices()[2],  karray[grid.x>0.4| units.m])  
        
    def test13(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEquals(grid[0].shape, (4,2))
        self.assertEquals(grid[0][0].shape, (2,))
        self.assertEquals(grid[...,2,1].shape, (5,))
        
    def test14(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEquals(grid[0].x, grid.x[0])
        self.assertEquals(grid[0][1].x, grid.x[0][1])
        self.assertEquals(grid[1][2][1].x, grid.x[1][2][1])
        self.assertEquals(grid[...,2,1].x, grid.x[...,2,1])
        self.assertEquals(grid[1,...,1].x, grid.x[1,...,1])
        self.assertEquals(grid[1,2,...].x, grid.x[1,2,...])
        self.assertEquals(grid[...,...,1].x, grid.x[...,...,1])
        self.assertEquals(grid[2,...,...].x, grid.x[2,...,...])
        self.assertEquals(grid[...,3,...].x, grid.x[...,3,...])
        self.assertEquals(grid[...,3,...].y, grid.y[...,3,...])
        self.assertEquals(grid[...,3,...].z, grid.z[...,3,...])
    
    def test15(self):
        
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        nk = nj = ni =0
        for plane1 in grid:
            nk += 1
            for plane2 in plane1:
                nj += 1
                for plane3 in plane2:
                    ni += 1
        self.assertEquals(nk, 5)
        self.assertEquals(nj, 4 * 5)
        self.assertEquals(ni, 2 * 4 * 5)
                    
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
        
        
class TestGridAttributes(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        print grid.get_minimum_position()
        self.assertAlmostRelativeEquals(grid.get_minimum_position()+ (1.0  | units.m),  ([0.0, 0.0, 0.0] | units.m) + (1.0  | units.m))
        self.assertAlmostRelativeEquals(grid.get_maximum_position(),  [1.0, 1.0, 1.0] | units.m)
        self.assertAlmostRelativeEquals(grid.get_volume(),  1.0 | units.m ** 3)
        self.assertTrue(grid.contains([0.5,0.5,0.5] | units.m))
        
    def test2(self):
        grid = core.Grid.create((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertTrue(numpy.all(grid.contains([[0.5,0.5,0.5] , [0.1,0.1,0.1]]| units.m)))
        self.assertFalse(numpy.all(grid.contains([[1.1,0.5,0.5] , [0.1,1.1,0.1]]| units.m)))
        
class TestGridSampling(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        self.assertEquals(sample.index , [1,1,1])
        sample = grid.samplePoint([2.5,2.5,2.5]| units.m)
        self.assertEquals(sample.index , [1,1,1])
        sample = grid.samplePoint([3.5,3.5,3.5]| units.m)
        self.assertEquals(sample.index , [1,1,1])
        
        for x in range(0,200):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m)
            self.assertEquals(sample.index , [0,2,3])
            
        for x in range(200,400):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m)
            print sample.index
            self.assertEquals(sample.index , [1,3,4])

    def test2(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        self.assertEquals(sample.index_for_000_cell , [1,1,1])
        sample = grid.samplePoint([2.5,2.5,2.5]| units.m)
        self.assertEquals(sample.index_for_000_cell , [0,0,0])
        sample = grid.samplePoint([3.5,3.5,3.5]| units.m)
        self.assertEquals(sample.index_for_000_cell , [1,1,1])
        sample = grid.samplePoint([4.5,4.5,4.5]| units.m)
        self.assertEquals(sample.index_for_000_cell , [1,1,1])
        self.assertEquals(sample.index , [2,2,2])
        
        for x in range(0,100):
            
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m)
            self.assertEquals(sample.index_for_000_cell , [-1,1,2])
        for x in range(100,300):
            
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m)
            self.assertEquals(sample.index_for_000_cell , [0,2,3])
            
        for x in range(300,400):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m)
            print sample.index
            self.assertEquals(sample.index_for_000_cell , [1,3,4])
            
    
    def test3(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        self.assertEquals(sample.index_for_000_cell , [1,1,1])
        self.assertEquals(sample.surrounding_cell_indices , [
            [1,1,1],
            [2,1,1],
            [1,2,1],
            [1,1,2],
            [2,1,2],
            [1,2,2],
            [2,2,1],
            [2,2,2],
        ])
    
    def test4(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        print sample.surrounding_cells
        self.assertEquals(sample.surrounding_cells[0].position , [3.0,3.0,3.0] | units.m )    
        self.assertEquals(sample.surrounding_cells[1].position , [5.0,3.0,3.0] | units.m )   
        self.assertEquals(sample.surrounding_cells[-1].position , [5.0,5.0,5.0] | units.m )        

    
    def test5(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        masses = sample.get_values_of_attribute("mass")
        print masses
        self.assertEquals(masses[0] , 3.0 | units.kg ) 
        self.assertEquals(masses[1] , 5.0 | units.kg ) 
        self.assertEquals(masses[2] , 3.0 | units.kg )  
        self.assertEquals(masses[3] , 3.0 | units.kg )  
        self.assertEquals(masses[4] , 5.0 | units.kg )  
        self.assertEquals(masses[5] , 3.0 | units.kg )  
        self.assertEquals(masses[6] , 5.0 | units.kg )  
        self.assertEquals(masses[7] , 5.0 | units.kg ) 
        factors = sample.weighing_factors
        print factors
        self.assertEquals(factors[0] , 1.0 | units.none ) 
        self.assertEquals(factors[1] , 0.0 | units.none ) 
        self.assertEquals(factors[2] , 0.0 | units.none )  
        self.assertEquals(factors[3] , 0.0 | units.none )  
        self.assertEquals(factors[4] , 0.0 | units.none )  
        self.assertEquals(factors[5] , 0.0 | units.none )  
        self.assertEquals(factors[6] , 0.0 | units.none )  
        self.assertEquals(factors[7] , 0.0 | units.none ) 
        
        self.assertAlmostRelativeEquals(sample.mass , 3.0 | units.kg ) 
            
    def test6(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        for xpos in numpy.arange(3.0,5.0,0.1):
            sample = grid.samplePoint([xpos,3.0,3.0]| units.m)
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,3.0,3.0]| units.m)
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,5.0,3.0]| units.m)
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,3.0,5.0]| units.m)
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
        

            sample = grid.samplePoint([4.0,4.0,4.0]| units.m)
            self.assertAlmostRelativeEquals(sample.mass , (4.0 | units.kg)) 

    def test7(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m)
        self.assertTrue(sample.isvalid)
        sample = grid.samplePoint([11.0,3.0,3.0]| units.m)
        self.assertFalse(sample.isvalid)
        sample = grid.samplePoint([3.0,-1.0,3.0]| units.m)
        self.assertFalse(sample.isvalid)
        
    
    def test8(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m, must_return_values_on_cell_center = True)    
        self.assertEquals(sample.position, [3.0,3.0,3.0]| units.m)    
        self.assertEquals(sample.mass, 3.0 | units.kg)
        sample = grid.samplePoint([3.5,3.0,3.0]| units.m, must_return_values_on_cell_center = True)    
        self.assertEquals(sample.position, [3.0,3.0,3.0]| units.m)  
        self.assertEquals(sample.mass, 3.0 | units.kg)
        
    def test9(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m, must_return_values_on_cell_center = False)    
        self.assertEquals(sample.position, [3.0,3.0,3.0]| units.m)    
        self.assertEquals(sample.mass, 3.0 | units.kg)
        sample = grid.samplePoint([3.5,3.0,3.0]| units.m, must_return_values_on_cell_center = False)    
        self.assertEquals(sample.position, [3.5,3.0,3.0]| units.m)  
        self.assertEquals(sample.mass, 3.5 | units.kg)
        
    
class TestGridSamplingMultiplePoints(amusetest.TestCase):
    
    def test1(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        samples = grid.samplePoints([[3.0,3.0,3.0], [4.0,3.0,3.0]]| units.m)
        self.assertEquals(len(samples), 2)
        self.assertEquals(samples.position[0] , [3.0,3.0,3.0]| units.m)
        self.assertEquals(samples.position[0] , samples[0].position)
        self.assertEquals(samples.position[1] , samples[1].position)
        self.assertEquals(samples.mass , [3.0, 4.0] | units.kg)
    
    def test2(self):
        grid = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        samples = grid.samplePoints([[3.5,3.0,3.0], [4.5,3.0,3.0]]| units.m)
        self.assertEquals(len(samples), 2)
        self.assertEquals(samples.mass , [3.5, 4.5] | units.kg)
        
    def test3(self):
        grid1 = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid1.mass = grid1.x.value_in(units.m) | units.kg
        grid2 = core.Grid.create((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid2.position += (10.0,0,0) | units.m
        grid2.mass = grid2.x.value_in(units.m) | units.kg
        samples = SamplePointsOnMultipleGrids((grid1, grid2), [[3.0,3.0,3.0], [4.0,3.0,3.0], [13,3,3]]| units.m)
        self.assertEquals(len(samples), 3)
        self.assertEquals(samples.mass , [3.0, 4.0, 13.0] | units.kg)
