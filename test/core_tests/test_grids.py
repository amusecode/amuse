from amuse.test import amusetest

from amuse.support.interface import InCodeComponentImplementation

from amuse.datamodel.indexing import *

from amuse.datamodel.grids import *

import numpy
import inspect
import collections
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel

class TestGrids(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        self.assertEqual(grid.mass[0][1][2], 2.0 | units.kg)
        self.assertEqual(grid[0][1][2].mass, 2.0 | units.kg)
        self.assertEqual(len(grid.mass), 5)
        
    def test2(self):
        grid = datamodel.Grid(5,4,3)
        grid.mass = units.kg.new_quantity(numpy.arange(5*4*3).reshape(5,4,3))
        self.assertEqual(grid.number_of_dimensions(), 3)
        subgrid = grid[1]
        self.assertEqual(subgrid.number_of_dimensions(), 2)
        self.assertEqual(subgrid.mass.number.shape, (4,3))
    
    def test3(self):
        grid = datamodel.Grid(5,4,3)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        self.assertEqual(grid.number_of_dimensions(), 3)
       
        subgrid = grid[1][2]
        self.assertEqual(subgrid.number_of_dimensions(), 1)
        self.assertEqual(subgrid.mass.number.shape, (3,))
        self.assertTrue(numpy.all(values[1][2] == subgrid.mass.value_in(units.kg)))
        
    def test4(self):
        grid = datamodel.Grid(5,4,3)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        self.assertEqual(grid.number_of_dimensions(), 3)
       
        gridpoint = grid[1][2][1]
        self.assertEqual(gridpoint.mass, 19.0 | units.kg)
        gridpoint = grid[1][2][2]
        self.assertEqual(gridpoint.mass, 20.0 | units.kg)
        
    def test5(self):
        grid = datamodel.Grid(5,4,3)
        grid.add_calculated_attribute("squared_mass", lambda m : m * m, ["mass",])
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        gridpoint = grid[1][2][1]
        self.assertEqual(gridpoint.mass, 19.0 | units.kg)
        self.assertEqual(gridpoint.squared_mass, (19.0 | units.kg) ** 2)
        subgrid = grid[1][2]
        self.assertTrue(numpy.all(subgrid.squared_mass == ([18.0, 19.0, 20.0] | units.kg) ** 2))
        
    def test6(self):
        grid = datamodel.Grid(5,4,3)
        grid.add_function_attribute("sum_mass", lambda grid, x : grid.mass.sum() + x, lambda grid, gridpoint, x : gridpoint.mass + x)
        values = numpy.arange(5*4*3).reshape(5,4,3)
        grid.mass = units.kg.new_quantity(values)
        gridpoint = grid[1][2][1]
        self.assertEqual(gridpoint.mass, 19.0 | units.kg)
        self.assertEqual(gridpoint.sum_mass(2.0 | units.kg), (21.0 | units.kg) )
        subgrid = grid[1][2]
        self.assertTrue(numpy.all(subgrid.sum_mass(2 | units.kg) == (18.0 + 19.0 + 20.0 + 2.0 | units.kg)))
        
    def test7(self):
        grid = datamodel.Grid(5,4,3)
        grid.add_vector_attribute("position", ["x","y","z"])
        x = numpy.arange(5*4*3).reshape(5,4,3)
        y = x + 100.0
        z = x + 200.0
        grid.x = units.m.new_quantity(x)
        grid.y = units.m.new_quantity(y)
        grid.z = units.m.new_quantity(z)
        gridpoint = grid[1][2][1]
        self.assertEqual(gridpoint.position[0], 19 | units.m)
        self.assertEqual(gridpoint.position[1], 119 | units.m)
        self.assertEqual(gridpoint.position[2], 219 | units.m)
        subgrid = grid[1][2]
        self.assertEqual(subgrid.position[1][0], 19 | units.m)
        self.assertEqual(subgrid.position[1][1], 119 | units.m)
        self.assertEqual(subgrid.position[1][2], 219 | units.m)
        
    
    def test8(self):
        grid0 = datamodel.Grid(5,4,3)
        x = numpy.arange(5*4*3).reshape(5,4,3)
        y = x + 100.0
        grid0.x = units.m.new_quantity(x)
        grid0.y = units.m.new_quantity(y)
        
        grid1 = datamodel.Grid(5,4,3)
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
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[0][0][0].x, 0.1 | units.m)
        self.assertEqual(grid[0][0][0].y, 0.125 | units.m)
        self.assertEqual(grid[0][0][0].z, 0.25 | units.m)
        self.assertEqual(grid[...,0,0].x, [0.1,0.3,0.5,0.7,0.9] | units.m)
        self.assertEqual(grid[0,0,...].z, [0.25, 0.75] | units.m)
        
        cellsize = grid.cellsize()
        self.assertAlmostRelativeEquals(cellsize[0], 0.2 | units.m)
        self.assertAlmostRelativeEquals(cellsize[1], 0.25 | units.m)
        self.assertAlmostRelativeEquals(cellsize[2], 0.5 | units.m)
    
    

    def test11(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        iarray,jarray,karray = grid.indices()
        for i in range(5):
            for j in range(4):
                for k in range(2):
                    self.assertEqual(iarray[i][j][k], i)
                    self.assertEqual(jarray[i][j][k], j)
                    self.assertEqual(karray[i][j][k], k)
        iarray,jarray,karray = grid.indices()
        i = 0
        for j in range(4):
            for k in range(2):
                self.assertEqual(iarray[i][j][k], i)
                self.assertEqual(jarray[i][j][k], j)
                self.assertEqual(karray[i][j][k], k)
       
        iarray,jarray,karray = grid[...,0,0].indices()
        j = 0
        k = 0
        for i in range(5):
            self.assertEqual(iarray[i], i)
            self.assertEqual(jarray[i], j)
            self.assertEqual(karray[i], k)
        iarray,jarray,karray = grid[3,2,...].indices()
        i = 3
        j = 2
        for k in range(2):
            self.assertEqual(iarray[k], i)
            self.assertEqual(jarray[k], j)
            self.assertEqual(karray[k], k)
        iarray,jarray,karray = grid[2,...,1].indices()
        i = 2
        k = 1
        for j in range(4):
            self.assertEqual(iarray[j], i)
            self.assertEqual(jarray[j], j)
            self.assertEqual(karray[j], k)
    
    
    def test12(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[grid.x > 0.4| units.m].x.shape, (24,))
        self.assertEqual(grid[grid.x>0.4| units.m].x, grid.x[grid.x>0.4|units.m])
        iarray,jarray,karray = grid.indices()
        self.assertEqual(grid[grid.x>0.4| units.m].indices()[0],  iarray[grid.x>0.4| units.m])
        self.assertEqual(grid[grid.x>0.4| units.m].indices()[1],  jarray[grid.x>0.4| units.m])
        self.assertEqual(grid[grid.x>0.4| units.m].indices()[2],  karray[grid.x>0.4| units.m])  
        
    def test13(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[0].shape, (4,2))
        self.assertEqual(grid[0][0].shape, (2,))
        self.assertEqual(grid[...,2,1].shape, (5,))
        
    def test14(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[0].x, grid.x[0])
        self.assertEqual(grid[0][1].x, grid.x[0][1])
        self.assertEqual(grid[1][2][1].x, grid.x[1][2][1])
        self.assertEqual(grid[...,2,1].x, grid.x[...,2,1])
        self.assertEqual(grid[1,...,1].x, grid.x[1,...,1])
        self.assertEqual(grid[1,2,...].x, grid.x[1,2,...])
        self.assertEqual(grid[...,1].x, grid.x[...,1])
        self.assertEqual(grid[2,...].x, grid.x[2,...])
        self.assertEqual(grid[:,3,:].x, grid.x[:,3,:])
        self.assertEqual(grid[:,3,:].y, grid.y[:,3,:])
        self.assertEqual(grid[:,3,:].z, grid.z[:,3,:])
    
    def test15(self):
        
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        nk = nj = ni =0
        for plane1 in grid:
            nk += 1
            for plane2 in plane1:
                nj += 1
                for plane3 in plane2:
                    ni += 1
        self.assertEqual(nk, 5)
        self.assertEqual(nj, 4 * 5)
        self.assertEqual(ni, 2 * 4 * 5)
                    
    def test16(self):
        grid = datamodel.new_regular_grid((5,4,2,4), [1.0 | units.m, 1.0 | units.m, 1.0 | units.m, 1.0 | units.s], ('x', 'y', 'z', 't') )
        self.assertEqual(grid.shape, (5,4,2,4))
        self.assertEqual(grid.x.shape, (5,4,2,4))
        self.assertAlmostRelativeEquals( grid[1][2][1].x, ([0.3] * 4) | units.m)
        self.assertAlmostRelativeEquals( grid[1][2][1].t, [0.125, 0.375, 0.625, 0.875] | units.s)
        self.assertEqual(grid[0].x, grid.x[0])
        self.assertEqual(grid[0][1].x, grid.x[0][1])
        self.assertEqual(grid[1][2][1].x, grid.x[1][2][1])
        self.assertEqual(grid[1][2][1][2].x, grid.x[1][2][1][2])

    def test17(self):
        grid = datamodel.new_regular_grid((4,2), [1.0 | units.m, 1.0 | units.m])
        self.assertEqual(grid.shape, (4,2) )
        self.assertEqual(grid.x.shape, (4,2))
        self.assertAlmostRelativeEquals( grid[1].x, ([0.375] * 2) | units.m)
        self.assertAlmostRelativeEquals( grid[1][1].y, 0.75 | units.m)
        self.assertEqual(grid[0].x, grid.x[0])
        self.assertEqual(grid[0][1].x, grid.x[0][1])
        
    def test18(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid.shape, (5,4,2) )
        self.assertEqual(grid[1:2,...].x.shape, (1,4,2 ))
        self.assertEqual(grid[1:2,...].shape, (1,4,2) )
        self.assertEqual(grid[1:2,...].x, grid.x[1:2,...])
        self.assertEqual(grid[1:3,...].x.shape, (2,4,2) )
        self.assertEqual(grid[1:3,...].shape, (2,4,2) )
        self.assertEqual(grid[1:3,...].x, grid.x[1:3,...])

    def test19(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[1:3,...].x, grid.x[1:3,...])
        self.assertEqual(grid[1:3,2:3,...].x, grid.x[1:3,2:3,...])
        self.assertEqual(grid[1:3,2:3,0:1].x, grid.x[1:3,2:3,0:1])
        self.assertEqual(grid[1:3,...,0:1].x, grid.x[1:3,...,0:1])
        self.assertEqual(grid[...,0:1].x, grid.x[...,0:1])
        self.assertEqual(grid[...,2:3,0:1].x, grid.x[...,2:3,0:1])
        
    def test20(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[1:3,:,:].x, grid.x[1:3,:,:])
        self.assertEqual(grid[1:3,2:3,:].x, grid.x[1:3,2:3,:])
        self.assertEqual(grid[1:3,2:3,0:1].x, grid.x[1:3,2:3,0:1])
        self.assertEqual(grid[1:3,:,0:1].x, grid.x[1:3,:,0:1])
        self.assertEqual(grid[:,:,0:1].x, grid.x[:,:,0:1])
        self.assertEqual(grid[:,2:3,0:1].x, grid.x[:,2:3,0:1])
        
        
    def test21(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid[1:3,:,:].copy().x, grid.x[1:3,:,:])
        self.assertEqual(grid[1:3,:,:].copy().shape, (2,4,2))
        self.assertEqual(grid[1:3,2:3,:].copy().x, grid.x[1:3,2:3,:])
        self.assertEqual(grid[1:3,2:3,:].copy().shape, (2,1,2))
        self.assertEqual(grid[1:3,2:3,0:1].copy().x, grid.x[1:3,2:3,0:1])
        self.assertEqual(grid[1:3,2:3,0:1].copy().shape, (2,1,1))
        self.assertEqual(grid[1:3,:,0:1].copy().x, grid.x[1:3,:,0:1])
        self.assertEqual(grid[1:3,:,0:1].copy().shape, (2,4,1))
        self.assertEqual(grid[:,:,0:1].copy().x, grid.x[:,:,0:1])
        self.assertEqual(grid[:,:,0:1].copy().shape, (5,4,1))
        self.assertEqual(grid[:,2:3,0:1].copy().x, grid.x[:,2:3,0:1])
        self.assertEqual(grid[:,2:3,0:1].copy().shape, (5,1,1))
        
    
    def test22(self):
        grid1 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        slice1 = grid1[1:3,:,:].copy()
        slice2 = grid2[1:3,:,:]
        slice1.x = -10 | units.m
        channel = slice1.new_channel_to(slice2)
        channel.copy()
        self.assertEqual(grid2.x[1:3], -10 | units.m)
        self.assertEqual(grid2.x[4],grid1.x[4])
    
    
    def test23(self):
        grid1 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual(grid1[1:3,...].shape,(2,4,2))
        slice1 = grid1[1:3,...].copy()
        slice2 = grid2[1:3,...]
        slice1.x = -10 | units.m
        channel = slice1.new_channel_to(slice2)
        channel.copy()
        self.assertEqual(grid2.x[1:3], -10 | units.m)
        self.assertEqual(grid2.x[4],grid1.x[4])
        
    
    def test24(self):
        particle = datamodel.Particle()
        particle.mass = 10 | units.kg
        
        grid = datamodel.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        grid.nounit = 10
        self.assertEqual(grid.nounit[0][1][2], 10)
        self.assertEqual(grid[0][1][2].nounit, 10)
        self.assertEqual(len(grid.nounit), 5)
        #grid[0][1][0].particle = particle
        #self.assertEquals(grid.mass[0][1][2], 2.0 | units.kg)
        #self.assertEquals(grid[0][1][0].particle, particle)
        #self.assertEquals(grid[0][1][1].particle, None)


    def test25(self):
        grid = datamodel.Grid(5,4,3)
        grid.mass = 2.0 | units.kg
        for cell in grid.iter_cells():
            self.assertEqual(cell.mass, 2.0 | units.kg)
    
    def test26(self):
        
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
    
        xv, yv, zv = numpy.mgrid[0:5,0:4,0:2]
        xv = xv.flatten()
        yv = yv.flatten()
        zv = zv.flatten()
        i = 0
        for cell in grid.iter_cells():
            expected_position = grid[xv[i], yv[i], zv[i]].position
            self.assertEqual(cell.position, expected_position)
            i += 1

    def test27(self):
        grid = datamodel.new_regular_grid((3,3), [1.0, 1.0] | units.m)
        subgrid1=grid[0:1,0:2]
        subgrid2=grid[0:1,0:2]
        subgrid3=grid[0:2,0:3][0:1,0:2]
        subgrid4=grid[0:1,0:3]
                
        self.assertTrue(subgrid1==subgrid2)
        self.assertTrue(subgrid1==subgrid3)
        self.assertTrue(subgrid2==subgrid3)
        self.assertFalse(subgrid1==subgrid4)
        self.assertFalse(subgrid2==subgrid4)
        self.assertFalse(subgrid3==subgrid4)

    def test28(self):
        grid = datamodel.Grid(200,400)
        subgrid1=grid[1:-1,1:-1]
        subgrid2=subgrid1[3:5,1:399]
        self.assertEqual(subgrid2.shape,(2,397))

    def test29(self):
        grid = datamodel.Grid(200)
        subgrid1=grid[1:-1]
        subgrid2=subgrid1[3:5]
        self.assertEqual(subgrid2.shape,(2,))

    def test30(self):
        grid = datamodel.Grid(200)
        subgrid1=grid[1:199]
        subgrid2=subgrid1[3:5]
        self.assertEqual(subgrid2.shape,(2,))

    def test31(self):
        grid = datamodel.Grid(20,20,20)
        a=numpy.zeros((20,20,20))
        self.assertEqual(a[1].shape,grid[1].shape) 
        self.assertEqual(a[1,2].shape,grid[1,2].shape) 
        self.assertEqual(a[1:5].shape,grid[1:5].shape) 
        self.assertEqual(a[2,1:5].shape,grid[2,1:5].shape)
        self.assertEqual(a[2,...].shape,grid[2,...].shape)       
        self.assertEqual(a[...,3,:].shape,grid[...,3,:].shape)       
        self.assertEqual(a[...,3:5,:].shape,grid[...,3:5,:].shape)       
        self.assertEqual(a[...,3:5,:].shape,grid[...,3:5,:].shape)
        self.assertEqual(a[::,::2].shape,grid[::,::2].shape)        

    def test32(self):
        grid = datamodel.Grid(200)
        grid.mass=numpy.arange(200)

        self.assertEqual(grid[1].mass,grid.mass[1]) 
        self.assertEqual(grid[1:-1].mass,grid.mass[1:-1])
        self.assertEqual(grid[-1:1:-1].mass,grid.mass[-1:1:-1])
        self.assertEqual(grid[-1:1:-1][1:-1].mass,grid.mass[-1:1:-1][1:-1])
        self.assertEqual(grid[-1:1:-1][-1:1].mass,grid.mass[-1:1:-1][-1:1])
        self.assertEqual(grid[-1:1:-1][-1:1:-3].mass,grid.mass[-1:1:-1][-1:1:-3])
        self.assertEqual(grid[300:1:-2][-1:5:-3].mass,grid.mass[300:1:-2][-1:5:-3])
        self.assertEqual(grid[100::-2][::3].mass,grid.mass[100::-2][::3])
        self.assertEqual(grid[100::-2][::-3].mass,grid.mass[100::-2][::-3])


    def test32b(self):
        grid = datamodel.Grid(200)
        grid.mass=numpy.arange(200)

        self.assertEqual(grid[::-1].mass,grid.mass[::-1])
        self.assertEqual(grid[10::-1].mass,grid.mass[10::-1])
        self.assertEqual(grid[:100:2].mass,grid.mass[:100:2])        
        self.assertEqual(grid[-1::-1].mass,grid.mass[-1::-1])
        self.assertEqual(grid[-1:-300:-1].mass,grid.mass[-1:-300:-1])
        self.assertEqual(grid[300:-300:-1].mass,grid.mass[300:-300:-1])        
        self.assertEqual(grid[300:-300:-7].mass,grid.mass[300:-300:-7])        
        
        
        
    def test33(self):
        grid = datamodel.Grid(20)
        grid.mass=numpy.zeros((20,5))
        self.assertEqual(grid[1].mass,numpy.zeros(5))
        grid.mass=numpy.arange(5)
        self.assertEqual(grid[-1].mass,numpy.arange(5))
        subgrid=grid[::2]
        self.assertEqual(subgrid[-1].mass,numpy.arange(5))
        subgrid[1].mass=5-numpy.arange(5)
        self.assertEqual(subgrid[1].mass,5-numpy.arange(5))
        self.assertEqual(grid[2].mass,5-numpy.arange(5))

    def test34(self):
        grid = datamodel.Grid(20)
        grid.mass=numpy.zeros((20,5)) | units.kg
        self.assertEqual(grid[1].mass,numpy.zeros(5) | units.kg)
        grid.mass=numpy.arange(5) | units.kg
        self.assertEqual(grid[-1].mass,numpy.arange(5)| units.kg)
        subgrid=grid[::2]
        self.assertEqual(subgrid[-1].mass,numpy.arange(5)| units.kg)
        subgrid[1].mass=(5-numpy.arange(5))| units.kg
        self.assertEqual(subgrid[1].mass,(5-numpy.arange(5))| units.kg)
        self.assertEqual(grid[2].mass,(5-numpy.arange(5))| units.kg)

    def test35(self):
        grid=datamodel.Grid(10,5)
        grid.mass=numpy.zeros((10,5,3))
        self.assertEqual(grid[2,2].mass,numpy.zeros(3))
        grid[::3,::2].mass=numpy.arange(3)
        self.assertEqual(grid[3,2].mass,numpy.arange(3))

    def test36(self):
        grid=datamodel.Grid(10)
        grid.mass=numpy.zeros((10,5,3))
        self.assertEqual(grid[2].mass,numpy.zeros((5,3)))
        grid[::3].mass=numpy.ones((5,3))
        self.assertEqual(grid[3].mass,numpy.ones((5,3)))

    def test37(self):
        grid = datamodel.Grid(20)
        grid.mass=numpy.zeros((20,5))
        grid.mass=numpy.arange(5)
        self.assertEqual(grid[-1].mass,numpy.arange(5))
        self.assertEqual(grid.mass.shape,(20,5))
        subgrid=grid[::2]
        self.assertEqual(subgrid[-1].mass,numpy.arange(5))
        subgrid[1].mass=5-numpy.arange(5)
        self.assertEqual(subgrid[1].mass,5-numpy.arange(5))
        self.assertEqual(grid[2].mass,5-numpy.arange(5))
        cp=subgrid.copy()
        self.assertEqual(cp[1].mass,5-numpy.arange(5))
        self.assertEqual(cp.mass.shape,(10,5))       
        cp=grid.copy()
        self.assertEqual(cp.mass.shape,(20,5))
        self.assertEqual(cp[2].mass,5-numpy.arange(5))
        self.assertEqual(cp[-1].mass,numpy.arange(5))

    def test38(self):
        grid=datamodel.new_cartesian_grid((10,),1)
        sub=grid[::2]
        self.assertEqual(sub[0].x,0.5)
        self.assertEqual(sub[(0,)].x,0.5)

        grid=datamodel.new_cartesian_grid((10,10),1)
        sub=grid[::2,::]
        self.assertEqual(sub[0,0].x,0.5)
        self.assertEqual(sub[(0,1)].y,1.5)

    def test39(self):
        grid=datamodel.new_cartesian_grid((10,10),1)
        sub=grid[3:6,5:8]
        self.assertEqual(sub[0:-1,0:-1].x,sub.x[0:-1,0:-1])
        self.assertEqual(sub[0:-1,-1].x,sub.x[0:-1,-1])
        self.assertEqual(sub[-1,-1].x,sub.x[-1,-1])
        self.assertEqual(sub[-1,-2].x,sub.x[-1,-2])

    def test40(self):
        grid1 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        grid1.m1 = 1
        grid1.m2 = 2
        channel = grid1.new_channel_to(grid2)
        channel.transform(["m3"],lambda x,y: (x,),["m1","m2"])
        self.assertEqual(grid2.m3, 1)
        channel.transform(["m3"],lambda x,y: (y,),["m1","m2"])
        self.assertEqual(grid2.m3, 2)
        channel.transform(["m3"],lambda x,y: (x+y,),["m1","m2"])
        self.assertEqual(grid2.m3, 3)
        channel.transform(["m3","m4"],lambda x,y: (x+y,2*x-y),["m1","m2"])
        self.assertEqual(grid2.m3, 3)
        self.assertEqual(grid2.m4, 0)

    def test40b(self):
        grid = datamodel.new_regular_grid((50,), [1.0] | units.m)
        for index in [ [0], [0,3,4], [1,2,2],[[2,3]],[[0,1],[2,3]],list(range(50)) ]:
          i=numpy.array(index)
          self.assertEqual(grid[i].x, grid.x[ i ])
          self.assertEqual(grid[i][1:].x, grid.x[ i ][1:])
          self.assertEqual(grid[i][1::2].x, grid.x[ i ][1::2])

    def test41(self):
        grid = datamodel.new_regular_grid((10,10), [1.0,2.0] | units.m)
        for _i,_j in [ ([0,1],[2,3]) ]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(grid[i,j].x, grid.x[ i,j ])
          self.assertEqual(grid[i,j][1:].x, grid.x[ i,j ][1:])
          self.assertEqual(grid[i,j][1::2].x, grid.x[ i,j ][1::2])

    def test42(self):
        grid = datamodel.new_regular_grid((3,4,5,6), [1.0,2.0,3.0,4.0],axes_names="abcd")
        for _i,_j in [ ([0],[3]),([0,1],[2,3]) ]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(grid[i,j].a, grid.a[ i,j ])
          self.assertEqual(grid[0,1,i,j].a, grid.a[0,1, i,j ])
          self.assertEqual(grid[i,j,0,0].a, grid.a[ i,j,0,0 ])
          self.assertEqual(grid[i,1,2,j].a, grid.a[ i,1,2,j])

    def test43(self):
        grid = datamodel.new_regular_grid((3,4,5,6), [1.0,2.0,3.0,4.0],axes_names="abcd")
        for _i,_j in [ ([0,1],[2,3]) ]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          subgrid=grid[:,0,:,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])
          subgrid=grid[:,0,2,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])
          subgrid=grid[:,0,-2,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])

    def test44(self):
        grid = datamodel.new_regular_grid((3,4,5,6), [1.0,2.0,3.0,4.0],axes_names="abcd")
        for _i,_j in [ ([-2],[-4]),([0,-1],[-2,3]) ]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          subgrid=grid[:,0,:,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])
          subgrid=grid[:,0,2,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])
          subgrid=grid[:,0,-2,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])

    def test45(self):
        grid = datamodel.new_regular_grid((3,4,5,6), [1.0,2.0,3.0,4.0],axes_names="abcd")
        for _i,_j in [ ([-1],[-4]) ]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          subgrid=grid[::2,0,:,:]
          self.assertEqual(subgrid[i,j].a, subgrid.a[ i,j ])

    def test46(self):
        grid = datamodel.new_regular_grid((6,), [1.0],axes_names="abcd")
        subgrid=grid[::2]
        self.assertEqual(subgrid[-2].a, subgrid.a[ -2 ])

        grid = datamodel.new_regular_grid((7,), [1.0],axes_names="abcd")
        subgrid=grid[::2]
        self.assertEqual(subgrid[-2].a, subgrid.a[ -2 ])

        grid = datamodel.new_regular_grid((7,), [1.0],axes_names="abcd")
        subgrid=grid[6:0:-2]
        self.assertEqual(subgrid[-2].a, subgrid.a[ -2 ])

    def test47(self):
        grid = datamodel.Grid()
        grid.mass=12.
        self.assertEqual(grid[...].mass,12.)
        self.assertEqual(grid.mass,12.)

    def test48(self):
        p=datamodel.Grid(3)

        p.a1=1.
        p.a2=1. | units.rad
        p.a3=1. | units.deg
        
        # test all combinations:
        
        p[0].a1=2.
        p[0].a2=2.
        p[0].a3=2.
                
        p[1].a1=2. | units.rad
        p[1].a2=2. | units.rad
        p[1].a3=2. | units.rad
        
        p[2].a1=2. | units.deg
        p[2].a2=2. | units.deg
        p[2].a3=2. | units.deg
        
        self.assertEqual( p.a1, [2.,2., (2. | units.deg).value_in(units.none)])
        self.assertEqual( p.a2, [2.,2., (2. | units.deg).value_in(units.none)])
        self.assertEqual( p.a3, [(2. | units.rad).in_(units.deg),
                                  (2. | units.rad).in_(units.deg) , 2. | units.deg])


    def test50(self):
        grid = datamodel.Grid(3,2)
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            grid.savepoint((i + 1) * 1.0 | units.s)

        dens = grid.get_timeline_of_attribute("density")
        self.assertEqual(len(dens), 3)
        self.assertEqual(dens[0][1].shape, (3,2))
        self.assertEqual(dens[1][1].shape, (3,2))
        self.assertEqual(dens[2][1].shape, (3,2))
        self.assertEqual(dens[0][0], 1.0 | units.s)
        self.assertEqual(dens[1][0], 2.0 | units.s)
        self.assertEqual(dens[2][0], 3.0 | units.s)

    def test51(self):
        grid = datamodel.Grid(3,2)
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            grid.savepoint((i + 1) * 1.0 | units.s)

        dens = grid[0:1,:].get_timeline_of_attribute("density")
        self.assertEqual(len(dens), 3)
        self.assertEqual(dens[0][1].shape, (1,2))
        self.assertEqual(dens[1][1].shape, (1,2))
        self.assertEqual(dens[2][1].shape, (1,2))
        self.assertEqual(dens[0][0], 1.0 | units.s)
        self.assertEqual(dens[1][0], 2.0 | units.s)
        self.assertEqual(dens[2][0], 3.0 | units.s)

    def test51(self):
        grid = datamodel.Grid(3,2)
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            grid.savepoint((i + 1) * 1.0 | units.s)

        time,dens = grid[0:1,:].get_timeline_of_attribute_as_vector("density")
        self.assertEqual(dens.shape, (3,1,2))
        self.assertEqual(dens[0].shape, (1,2))
        self.assertEqual(dens[1].shape, (1,2))
        self.assertEqual(dens[2].shape, (1,2))
        self.assertEqual(time[0], 1.0 | units.s)
        self.assertEqual(time[1], 2.0 | units.s)
        self.assertEqual(time[2], 3.0 | units.s)

    def test52(self):
        grid = datamodel.Grid(3,2)
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            grid.savepoint((i + 1) * 1.0 | units.s)

        time,dens = grid[2,1].get_timeline_of_attribute_as_vector("density")
        self.assertEqual(dens.shape, (3,))
        self.assertEqual(dens[0], 0 | units.kg/units.m**3)
        self.assertEqual(dens[1], 1 | units.kg/units.m**3)
        self.assertEqual(dens[2], 2 | units.kg/units.m**3)
        self.assertEqual(time[0], 1.0 | units.s)
        self.assertEqual(time[1], 2.0 | units.s)
        self.assertEqual(time[2], 3.0 | units.s)

    def test53(self):
        grid = datamodel.Grid(3,2)
        subgrid=grid[0:1,:]
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            subgrid.savepoint((i + 1) * 1.0 | units.s)

        time,dens = subgrid.get_timeline_of_attribute_as_vector("density")
        self.assertEqual(dens.shape, (3,1,2))
        self.assertEqual(dens[0], 0 | units.kg/units.m**3)
        self.assertEqual(dens[1], 1 | units.kg/units.m**3)
        self.assertEqual(dens[2], 2 | units.kg/units.m**3)
        self.assertEqual(time[0], 1.0 | units.s)
        self.assertEqual(time[1], 2.0 | units.s)
        self.assertEqual(time[2], 3.0 | units.s)

    def test54(self):
        """
        illustrates getting subgrid/gridpoint with history from subgrid with history
        """
        grid = datamodel.Grid(3,2)
        subgrid=grid[0:1,:]
        for i in range(3):
            grid.density = (i * 1.0) | units.kg/units.m**3
            subgrid.savepoint((i + 1) * 1.0 | units.s)
        # if the gridpoint derives directly from subgrid, its defined on the original
        # grid (which has no history...)
        subsub=subgrid.savepoint((i + 1) * 1.0 | units.s)
        time,dens = subsub[0,1].get_timeline_of_attribute_as_vector("density")
        self.assertEqual(dens.shape, (3,))
        self.assertEqual(dens[0], 0 | units.kg/units.m**3)
        self.assertEqual(dens[1], 1 | units.kg/units.m**3)
        self.assertEqual(dens[2], 2 | units.kg/units.m**3)
        self.assertEqual(time[0], 1.0 | units.s)
        self.assertEqual(time[1], 2.0 | units.s)
        self.assertEqual(time[2], 3.0 | units.s)

class TestGridFactories(amusetest.TestCase):
    def test1(self):
        grid1 = datamodel.new_cartesian_grid( (4,5), 1.0 | units.m)
        grid2 = datamodel.new_regular_grid( (4,5), [4.0,5.0] | units.m)
        grid3 = datamodel.new_rectilinear_grid( (4,5), [numpy.arange(5.) | units.m,numpy.arange(6.) | units.m])
        
        self.assertEqual(grid1.position,grid2.position)
        self.assertEqual(grid2.position,grid3.position)

    def test2(self):
        grid=datamodel.new_rectilinear_grid((10,),(1.*numpy.arange(11),))
        self.assertEqual(grid._axes_cell_boundaries,1.*numpy.arange(11))
        grid=datamodel.new_regular_grid((10,),[10.])
        self.assertEqual(grid._lengths,[10.])
        grid=datamodel.new_cartesian_grid((10,),1.)
        self.assertEqual(grid._cellsize,1.)
        grid=datamodel.new_regular_grid((10,20,),[10.,15.])
        self.assertEqual(grid._lengths,[10.,15.])

    def test3(self):
        N=10
        x,y=numpy.indices((N+1,N+1))
        grid=datamodel.new_structured_grid((N,N),[x,y])
        self.assertEqual(grid.shape,(N,N))
        x,y=numpy.indices((N,N))
        x=x+0.5
        y=y+0.5
        self.assertEqual(grid.x,x)
        self.assertEqual(grid.y,y)
        
    def test4(self):
        N=2
        x,y,z=numpy.indices((N+1,N+1,2*N+1))
        grid=datamodel.new_structured_grid((N,N,2*N),[x,y,z])
        self.assertEqual(grid.shape,(N,N,2*N))
        x,y,z=numpy.indices((N,N,2*N))
        x=x+0.5
        y=y+0.5
        z=z+0.5
        self.assertEqual(grid.x,x)
        self.assertEqual(grid.y,y)
        self.assertEqual(grid.z,z)



class TestGridAttributes(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertAlmostRelativeEquals(grid.get_minimum_position(),  ([0.0, 0.0, 0.0] | units.m) )
        self.assertAlmostRelativeEquals(grid.get_maximum_position(),  [1.0, 1.0, 1.0] | units.m)
        self.assertAlmostRelativeEquals(grid.get_volume(),  1.0 | units.m ** 3)
        self.assertTrue(grid.contains([0.5,0.5,0.5] | units.m))
        
    def test2(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertTrue(numpy.all(grid.contains([[0.5,0.5,0.5] , [0.1,0.1,0.1]]| units.m)))
        self.assertFalse(numpy.all(grid.contains([[1.1,0.5,0.5] , [0.1,1.1,0.1]]| units.m)))
        
    def test3(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        points = grid.points()
        self.assertEqual(points.shape, (6,5,3, 3))
        self.assertAlmostRelativeEquals(points[0][0][0],  ([0.0,0.0, 0.0] | units.m) )
        self.assertAlmostRelativeEquals(points[1][0][0] ,  ([0.2,0.0, 0.0] | units.m) )
        self.assertAlmostRelativeEquals(points[1][1][1],  [0.2,0.25, 0.5] | units.m )
        self.assertAlmostRelativeEquals(points[0][-1][-1] ,  ([0.0, 1.0, 1.0] | units.m) )
        self.assertAlmostRelativeEquals(points[-1][0][-1] ,  ([1.0, 0.0, 1.0] | units.m) )
        self.assertAlmostRelativeEquals(points[-1][-1][0] ,  ([1.0, 1.0, 0.0] | units.m) )
        self.assertAlmostRelativeEquals(points[-1][-1][-1] ,  ([1.0,1.0, 1.0] | units.m) )
        
    
    def test4(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        points = grid.points().reshape([6*5*3,3])
        connectivity = grid.connectivity()
        
        self.assertEqual(connectivity.shape, (5,4,2, 8))
        first_cell = connectivity[0][0][0]
        self.assertAlmostRelativeEquals(points[first_cell[0]], [0,0,0 ] | units.m)
        self.assertAlmostRelativeEquals(points[first_cell[1]], [0.2,0,0] | units.m)
        self.assertAlmostRelativeEquals(points[first_cell[2]] ,  ([0,0.25,0] | units.m))
        self.assertAlmostRelativeEquals(points[first_cell[3]] ,  ([0.2,0.25,0] | units.m))
        self.assertAlmostRelativeEquals(points[first_cell[4]] ,  ([0.0,0.0,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[first_cell[5]] ,  ([0.2,0.0,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[first_cell[6]] ,  ([0.0,0.25,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[first_cell[7]] ,  ([0.2,0.25,0.5] | units.m))
        
        self.assertEqual(connectivity[0][0][0], [ 0,15,  3, 18,  1, 16, 4, 19])
        
        
    
    def test5(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        points = grid.points().reshape([6*5*3,3])
        connectivity = grid.connectivity()
        
        self.assertEqual(connectivity.shape, (5,4,2, 8))
        cell = connectivity[0][0][1]
        self.assertAlmostRelativeEquals(points[cell[0]] ,  ([0.0,0.0,0.5 ] | units.m))
        self.assertAlmostRelativeEquals(points[cell[1]] ,  ([0.2,0,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[2]] ,  ([0,0.25,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[3]] ,  ([0.2,0.25,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[4]] ,  ([0.0,0.0,1.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[5]] ,  ([0.2,0.0,1.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[6]] ,  ([0.0,0.25,1.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[7]] ,  ([0.2,0.25,1.0] | units.m))
        
        self.assertEqual(connectivity[0][0][0], [ 0,15,  3, 18,  1, 16, 4, 19])
        
    
    
    def test6(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        points = grid.points().reshape([6*5*3,3])
        connectivity = grid.connectivity()
        
        self.assertEqual(connectivity.shape, (5,4,2, 8))
        cell = connectivity[1][1][1]
        self.assertAlmostRelativeEquals(points[cell[0]], ([0.2, 0.25, 0.5]|units.m) + ([0.0,0.0,0.0 ] | units.m))
        self.assertAlmostRelativeEquals(points[cell[1]], ([0.2, 0.25, 0.5]|units.m) + ([0.2,0,0.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[2]], ([0.2, 0.25, 0.5]|units.m) + ([0,0.25,0.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[3]], ([0.2, 0.25, 0.5]|units.m) + ([0.2,0.25,0.0] | units.m))
        self.assertAlmostRelativeEquals(points[cell[4]], ([0.2, 0.25, 0.5]|units.m) + ([0.0,0.0,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[5]], ([0.2, 0.25, 0.5]|units.m) + ([0.2,0.0,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[6]], ([0.2, 0.25, 0.5]|units.m) + ([0.0,0.25,0.5] | units.m))
        self.assertAlmostRelativeEquals(points[cell[7]], ([0.2, 0.25, 0.5]|units.m) + ([0.2,0.25,0.5] | units.m))
        
        self.assertEqual(connectivity[0][0][0], [ 0,15,  3, 18,  1, 16, 4, 19])
    
    def test7(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        
        self.assertAlmostRelativeEquals(grid[1][2][3].position, [3,5,7] |units.m)
        
        grid[1][2][3].position = [7,5,3] |units.m
        self.assertAlmostRelativeEquals(grid[1][2][3].position, [7,5,3] |units.m)
        
        
        grid[1][2][3].position += [1,2,3] |units.m
        self.assertAlmostRelativeEquals(grid[1][2][3].position, [8,7,6] |units.m)

    def test8(self):
        grid = datamodel.new_regular_grid((5,4), [1.0, 1.0] | units.m)
        self.assertAlmostRelativeEquals(grid.get_minimum_position(),  ([0.0, 0.0] | units.m) )
        self.assertAlmostRelativeEquals(grid.get_maximum_position(),  [1.0, 1.0] | units.m)
        self.assertAlmostRelativeEquals(grid.get_volume(),  1.0 | units.m ** 2)
        self.assertTrue(grid.contains([0.5,0.5] | units.m))

    def test9(self):
        grid = datamodel.new_regular_grid((5,4,2), [1.0, 1.0, 1.0] | units.m)
        self.assertEqual((0,0,0),grid.get_minimum_index())
        self.assertEqual((4,3,1),grid.get_maximum_index())

    def test10(self):
        grid1 = datamodel.new_regular_grid((5,4), [1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((5,4), [.1, .1] | units.m)
        grid3 = datamodel.new_regular_grid((5,4), [.1, .1] | units.m,offset=[0.5,0.6] | units.m)
        self.assertTrue(grid1.overlaps(grid2))
        self.assertTrue(grid1.overlaps(grid3))
        self.assertFalse(grid2.overlaps(grid3))
        self.assertTrue(grid2.overlaps(grid1))
        self.assertTrue(grid3.overlaps(grid1))
        self.assertFalse(grid3.overlaps(grid2))

    def test11(self):
        grid1 = datamodel.new_regular_grid((4,4), [1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((4,4), [1.0, 1.0] | units.m,offset=[-0.5,-0.5] | units.m)
        self.assertTrue(grid1.overlaps(grid2))
        overlap=grid1.get_overlap_with(grid2)
        self.assertEqual(overlap.position,grid1[0:3,0:3].position)

    def test12(self):
        grid1 = datamodel.new_regular_grid((4,4), [1.0, 1.0] | units.m)
        grid2 = datamodel.new_regular_grid((4,4), [1.0, 1.0] | units.m,offset=[-0.5,-0.5] | units.m)
        self.assertTrue(grid1.overlaps(grid2))
        overlap=grid1.get_overlap_with(grid2,eps=grid2.cellsize()[0]*1.e-6)
        self.assertEqual(overlap.position,grid1[0:2,0:2].position)


class TestGridSampling(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        self.assertEqual(sample.index , [1,1,1])
        sample = grid.samplePoint([2.5,2.5,2.5]| units.m,method="interpolation")
        self.assertEqual(sample.index , [1,1,1])
        sample = grid.samplePoint([3.5,3.5,3.5]| units.m,method="interpolation")
        self.assertEqual(sample.index , [1,1,1])
        
        for x in range(0,200):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m,method="interpolation")
            self.assertEqual(sample.index , [0,2,3])
            
        for x in range(200,400):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m,method="interpolation")
            self.assertEqual(sample.index , [1,3,4])

    def test2(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        self.assertEqual(sample.index_for_000_cell , [1,1,1])
        sample = grid.samplePoint([2.5,2.5,2.5]| units.m,method="interpolation")
        self.assertEqual(sample.index_for_000_cell , [0,0,0])
        sample = grid.samplePoint([3.5,3.5,3.5]| units.m,method="interpolation")
        self.assertEqual(sample.index_for_000_cell , [1,1,1])
        sample = grid.samplePoint([4.5,4.5,4.5]| units.m,method="interpolation")
        self.assertEqual(sample.index_for_000_cell , [1,1,1])
        self.assertEqual(sample.index , [2,2,2])
        
        for x in range(0,100):
            
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m,method="interpolation")
            self.assertEqual(sample.index_for_000_cell , [-1,1,2])
        for x in range(100,300):
            
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m,method="interpolation")
            self.assertEqual(sample.index_for_000_cell , [0,2,3])
            
        for x in range(300,400):
            sample = grid.samplePoint([0.0 + (x/100.0),4.0+(x/100.0),6.0+(x/100.0)]| units.m,method="interpolation")
            self.assertEqual(sample.index_for_000_cell , [1,3,4])
            
    
    def test3(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        self.assertEqual(sample.index_for_000_cell , [1,1,1])
        self.assertEqual(sample.surrounding_cell_indices , [
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
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        self.assertEqual(sample.surrounding_cells[0].position , [3.0,3.0,3.0] | units.m )    
        self.assertEqual(sample.surrounding_cells[1].position , [5.0,3.0,3.0] | units.m )   
        self.assertEqual(sample.surrounding_cells[-1].position , [5.0,5.0,5.0] | units.m )        

    
    def test5(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        masses = sample.get_values_of_attribute("mass")
        self.assertEqual(masses[0] , 3.0 | units.kg ) 
        self.assertEqual(masses[1] , 5.0 | units.kg ) 
        self.assertEqual(masses[2] , 3.0 | units.kg )  
        self.assertEqual(masses[3] , 3.0 | units.kg )  
        self.assertEqual(masses[4] , 5.0 | units.kg )  
        self.assertEqual(masses[5] , 3.0 | units.kg )  
        self.assertEqual(masses[6] , 5.0 | units.kg )  
        self.assertEqual(masses[7] , 5.0 | units.kg ) 
        factors = sample.weighing_factors
        self.assertEqual(factors[0] , 1.0 | units.none ) 
        self.assertEqual(factors[1] , 0.0 | units.none ) 
        self.assertEqual(factors[2] , 0.0 | units.none )  
        self.assertEqual(factors[3] , 0.0 | units.none )  
        self.assertEqual(factors[4] , 0.0 | units.none )  
        self.assertEqual(factors[5] , 0.0 | units.none )  
        self.assertEqual(factors[6] , 0.0 | units.none )  
        self.assertEqual(factors[7] , 0.0 | units.none ) 
        
        self.assertAlmostRelativeEquals(sample.mass , 3.0 | units.kg ) 
            
    def test6(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        for xpos in numpy.arange(3.0,5.0,0.1):
            sample = grid.samplePoint([xpos,3.0,3.0]| units.m,method="interpolation")
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,3.0,3.0]| units.m,method="interpolation")
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,5.0,3.0]| units.m,method="interpolation")
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
          
            sample = grid.samplePoint([xpos,3.0,5.0]| units.m,method="interpolation")
            self.assertAlmostRelativeEquals(sample.mass , (3.0 | units.kg) + ((2.0 * (xpos - 3.0) / 2.0) | units.kg) ) 
        

            sample = grid.samplePoint([4.0,4.0,4.0]| units.m,method="interpolation")
            self.assertAlmostRelativeEquals(sample.mass , (4.0 | units.kg)) 

    def test7(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m,method="interpolation")
        self.assertTrue(sample.isvalid)
        sample = grid.samplePoint([11.0,3.0,3.0]| units.m,method="interpolation")
        self.assertFalse(sample.isvalid)
        sample = grid.samplePoint([3.0,-1.0,3.0]| units.m,method="interpolation")
        self.assertFalse(sample.isvalid)
        
    
    def test8(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m, method="nearest")    
        self.assertEqual(sample.position, [3.0,3.0,3.0]| units.m)    
        self.assertEqual(sample.mass, 3.0 | units.kg)
        sample = grid.samplePoint([3.5,3.0,3.0]| units.m, method="nearest")    
        self.assertEqual(sample.position, [3.0,3.0,3.0]| units.m)  
        self.assertEqual(sample.mass, 3.0 | units.kg)
        
    def test9(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        sample = grid.samplePoint([3.0,3.0,3.0]| units.m, method="linear")    
        self.assertEqual(sample.position, [3.0,3.0,3.0]| units.m)    
        self.assertEqual(sample.mass, 3.0 | units.kg)
        sample = grid.samplePoint([3.5,3.0,3.0]| units.m, method="linear")    
        self.assertEqual(sample.position, [3.5,3.0,3.0]| units.m)  
        self.assertEqual(sample.mass, 3.5 | units.kg)
        
    
class TestGridSamplingMultiplePoints(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        samples = grid.samplePoints([[3.0,3.0,3.0], [4.0,3.0,3.0]]| units.m, method="linear")
        self.assertEqual(len(samples), 2)
        self.assertEqual(samples.position[0] , [3.0,3.0,3.0]| units.m)
        self.assertEqual(samples.position[0] , samples[0].position)
        self.assertEqual(samples.position[1] , samples[1].position)
        self.assertEqual(samples.mass , [3.0, 4.0] | units.kg)
    
    def test2(self):
        grid = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid.mass = grid.x.value_in(units.m) | units.kg
        samples = grid.samplePoints([[3.5,3.0,3.0], [4.5,3.0,3.0]]| units.m, method="linear")
        self.assertEqual(len(samples), 2)
        self.assertEqual(samples.mass , [3.5, 4.5] | units.kg)
        
    def test3(self):
        grid1 = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid1.mass = grid1.x.value_in(units.m) | units.kg
        grid2 = datamodel.new_regular_grid((5,5,5), [10.0, 10.0, 10.0] | units.m)
        grid2.position += (10.0,0,0) | units.m
        grid2.mass = grid2.x.value_in(units.m) | units.kg
        samples = SamplePointsOnMultipleGrids((grid1, grid2), [[3.0,3.0,3.0], [4.0,3.0,3.0], [13,3,3]]| units.m)
        self.assertEqual(len(samples), 3)
        self.assertEqual(samples.mass , [3.0, 4.0, 13.0] | units.kg)
