from amuse.test import amusetest

from amuse.datamodel.grids import *

import numpy

class TestGridAttributes(amusetest.TestCase):
    
    def test1(self):
        grid=new_cartesian_grid((10,20),1.)
        self.assertEqual(grid.cellsize(),[1.,1.])
        self.assertEqual(grid.get_minimum_index(), [0,0])
        self.assertEqual(grid.get_maximum_index(), [9,19])
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

    def xtest3(self):
        grid=new_rectilinear_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])

    def xtest4(self):
        grid=new_structured_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])

    def xtest5(self):
        grid=new_unstructured_grid((10,20),(10,10))
        self.assertEqual(grid.cellsize(),[1.,0.5])
