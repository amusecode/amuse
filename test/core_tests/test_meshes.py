from amuse.test import amusetest

from amuse.datamodel import grids
from amuse.datamodel import meshes



class TestGridWithLevel(amusetest.TestCase):
    
    def test1(self):
        grid = meshes.GridWithLevel(None, 1, (1,2))
        print grid.position_at_level(1)
        self.assertEquals(grid.position_at_level(1), (1,2))
        self.assertEquals(grid.position_at_level(2), (2,4))
        self.assertEquals(grid.position_at_level(3), (4,8))
        self.assertEquals(grid.position_at_level(4), (8,16))
        
    def test2(self):
        grid = meshes.GridWithLevel(None, 1, (0,0))
        self.assertEquals(grid.position_at_level(1), (0,0))
        self.assertEquals(grid.position_at_level(2), (0,0))
        self.assertEquals(grid.position_at_level(3), (0,0))
        
    def test3(self):
        grid = meshes.GridWithLevel(None, 3, (4,8))
        self.assertEquals(grid.position_at_level(1), (1,2))
        self.assertEquals(grid.position_at_level(2), (2,4))
        self.assertEquals(grid.position_at_level(3), (4,8))
        self.assertEquals(grid.position_at_level(4), (8,16))

    def test4(self):
        grid = meshes.GridWithLevel(None, 4, (8,16))
        self.assertEquals(grid.position_at_level(1), (1,2))
        self.assertEquals(grid.position_at_level(2), (2,4))
        self.assertEquals(grid.position_at_level(3), (4,8))
        self.assertEquals(grid.position_at_level(4), (8,16))
        
    def test5(self):
        grid = meshes.GridWithLevel(None, 4, (80,))
        self.assertEquals(grid.position_at_level(1), (10,))
    
    def test6(self):
        grid = meshes.GridWithLevel(grids.Grid(3,4), 1, (1,2))
        self.assertEquals(grid.shape_at_level(1), (3,4))
        self.assertEquals(grid.shape_at_level(2), (6,8))
        self.assertEquals(grid.shape_at_level(3), (12,16))
        self.assertEquals(grid.shape_at_level(4), (24,32))
    
    def test7(self):
        grid = meshes.GridWithLevel(None, 1, (1,2), factor_per_level=4)
        self.assertEquals(grid.position_at_level(1), (1,2))
        self.assertEquals(grid.position_at_level(2), (4,8))
        self.assertEquals(grid.position_at_level(3), (16,32))
        self.assertEquals(grid.position_at_level(4), (64,128))


class TestMultilevelBlockStructuredMesh(amusetest.TestCase):
    
    def test1(self):
        mesh = meshes.MultilevelBlockStructuredMesh((4,6), 1)
        
        self.assertEquals(mesh.shape, (4,6))
        self.assertEquals(mesh[0].shape, (6,))
        self.assertEquals(mesh[1].shape, (6,))
        
    def test2(self):
        mesh = meshes.MultilevelBlockStructuredMesh((4,6), 1)
        for i,j in (
                (0,0),
                (2,0),
                (0,3),
                (2,3)
            ):
            grid = grids.Grid(2,3)
            grid.i = i,
            grid.j = j,
            mesh.add_grid(grid, 1, i, j)
        self.assertEquals(mesh[1][2].i, 0)
        self.assertEquals(mesh[2][2].i, 2)
        self.assertEquals(mesh[2][4].i, 2)
        self.assertEquals(mesh[2][3].i, 2)
        self.assertEquals(mesh[2][3].j, 3)
        self.assertEquals(mesh[2][4].j, 3)
