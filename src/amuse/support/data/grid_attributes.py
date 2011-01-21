
from amuse.support.data import base
from amuse.support.data import grids

import numpy

grids.Grid.add_global_vector_attribute("position", ["x","y","z"])
grids.Grid.add_global_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])


        

@grids.AbstractGrid.function_for_set
def cellsize(grid, dimensions = 3):
    """Returns the lenght of each direction in the grid.
    Works on 3d, 2d and 1d grids, although the grid cells
    must have 3d positions (x,y and z)
    """
    cell1 = grid[0][0][0]
    cell2 = grid[1][0][0]
    dx = cell2.x - cell1.x
    
    result = [0.0,0.0,0.0] | dx.unit
    result[0] = dx
   
    if dimensions > 1:
        cell2 = grid[0][1][0]
        result[1] = cell2.y - cell1.y
    else:
        return result[0:1]
        
    if dimensions > 2:
        cell2 = grid[0][0][1]
        result[2] = cell2.z - cell1.z
        return result
    else:
        return result[0:2]




@grids.AbstractGrid.function_for_set
def get_minimum_position(grid):
    return grid[0,0,0].position - (0.5 * grid.cellsize())
    
@grids.AbstractGrid.function_for_set
def get_maximum_position(grid):
    return grid[-1,-1,-1].position + (0.5 * grid.cellsize())
    
@grids.AbstractGrid.function_for_set
def get_volume(grid):
    maximum_position = grid.get_maximum_position()
    minimum_position = grid.get_minimum_position()
    delta = maximum_position - minimum_position
    return delta.prod()
    

@grids.AbstractGrid.function_for_set
def contains(grid, points):
    return numpy.logical_and(
        numpy.all(points >= grid.get_minimum_position(), axis=len(points.shape)-1),
        numpy.all(points < grid.get_maximum_position(), axis=len(points.shape)-1)
    )
   
