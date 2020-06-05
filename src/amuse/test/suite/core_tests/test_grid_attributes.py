from amuse.test import amusetest

from amuse.datamodel.grids import *

import numpy

class TestGridAttributes(amusetest.TestCase):
    
    def test1(self):
        grid=new_cartesian_grid((10,20),1.)
        self.assertEqual(grid.cellsize(),[1.,1.])
        self.assertEqual(grid.get_minimum_index(), (0,0))
        self.assertEqual(grid.get_maximum_index(), (9,19))
        self.assertEqual(grid.get_minimum_position(), [0.,0.])
        self.assertEqual(grid.get_maximum_position(), [10.,20.])
        self.assertEqual(grid.get_volume(), 10*20.)
        self.assertEqual(grid.contains(numpy.asarray([[1,1],[-10,-10]])), [True,False])
        points=grid.points()
        x=points[:,:,0]
        y=points[:,:,1]
        self.assertEqual([x,y],numpy.indices((11,21)))

        grid2=new_cartesian_grid((10,20),1., offset=(5,5))
        self.assertTrue(grid.overlaps(grid2))
        grid3=new_cartesian_grid((10,20),1., offset=(15,25))
        self.assertFalse(grid.overlaps(grid3))
        overlap=grid.get_overlap_with(grid2)
        self.assertEqual(overlap.shape,(5,15))
    
    def test2(self):
        grid=new_regular_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])
        self.assertEqual(grid.get_minimum_index(), (0,0))
        self.assertEqual(grid.get_maximum_index(), (9,19))
        self.assertEqual(grid.get_minimum_position(), [0.,0.])
        self.assertEqual(grid.get_maximum_position(), [10.,10.])
        self.assertEqual(grid.get_volume(), 10*10.)
        self.assertEqual(grid.contains(numpy.asarray([[1,1],[-10,-10]])), [True,False])
        points=grid.points()
        x=points[:,:,0]
        y=points[:,:,1]
        xp,yp=numpy.indices((11,21))
        self.assertEqual(x,xp)
        self.assertEqual(y,yp*0.5)

        grid2=new_cartesian_grid((10,20),1., offset=(5,5))
        self.assertTrue(grid.overlaps(grid2))
        grid3=new_cartesian_grid((10,20),1., offset=(15,25))
        self.assertFalse(grid.overlaps(grid3))
        overlap=grid.get_overlap_with(grid2)
        self.assertEqual(overlap.shape,(5,10))


    def test3(self):
        x=numpy.arange(11)/10.
        y=numpy.arange(21)/20.
        grid=new_rectilinear_grid((10,20),(x,y**2))
        self.assertRaises(Exception, grid.cellsize,
          expected_message="a RectilinearGrid does not have a constant cellsize, use the cellsizes method instead")

    def xtest4(self):
        grid=new_structured_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])

    def xtest5(self):
        grid=new_unstructured_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])

    def test6(self):
        grid=new_regular_grid((10,10),(10,10))
        self.assertEqual(grid.get_index((6.2,3.7)), [6,3])
        self.assertEqual(grid.get_index(x=[6.2],y=[3.7]), [6,3])
        self.assertEqual(grid.get_index(y=[6.2],x=[3.7]), [3,6])
        self.assertEqual(grid.get_index(y=6.2,x=3.7)[0], 3)
        self.assertEqual(grid.get_index(y=6.2,x=3.7)[1], 6)
        
    def test7(self):
        grid=new_regular_grid((10,10),[20,10] | units.m,axes_names="ab")
        self.assertEqual(grid.get_index([16.2,3.7] | units.m), [8,3])
        self.assertEqual(grid.get_index(a=16.2 | units.m,b=3.7 | units.m), [8,3])
        self.assertEqual(grid.get_index(a=[16.2, 4.5] | units.m,b=[3.7,4.2] | units.m), [[8,3],[2,4]])

