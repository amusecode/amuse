from amuse.support.data import base
from amuse.support.data import grids
from amuse.units.quantities import zero

import numpy

grids.AbstractGrid.add_global_vector_attribute("position", ["x","y","z"])
grids.AbstractGrid.add_global_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])
grids.AbstractGrid.add_global_vector_attribute("magnetic_field", ["B1i","B2i","B3i"])


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
        if grid.shape[1] > 1:
            cell2 = grid[0][1][0]
            result[1] = cell2.y - cell1.y
    else:
        return result[0:1]
        
    if dimensions > 2:
        if grid.shape[2] > 1:
            cell2 = grid[0][0][1]
            result[2] = cell2.z - cell1.z
    else:
        return result[0:2]
    
    return result



@grids.AbstractGrid.function_for_set
def get_minimum_index(grid):
    return [0,0,0]
    
@grids.AbstractGrid.function_for_set
def get_maximum_index(grid):
    return grid.shape - numpy.asarray([1,1,1])
    

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
   

@grids.AbstractGrid.function_for_set
def points(grid):
    cellcenters = grid.position
    nx,ny,nz = grid.shape
         
    shape_with_boundary = numpy.asarray(cellcenters.shape) + 1
    shape_with_boundary[-1] -= 1
    result = cellcenters.zeros(shape_with_boundary, cellcenters.unit)
    dx = grid.cellsize() 
    half_dx = grid.cellsize() / 2.0
    
    
    result[1:,1:,1:] = grid.position + (half_dx * [1,1,1])
    result[:-1,1:,1:] = grid.position + (half_dx * [-1,1,1])
    result[1:,:-1,1:] = grid.position + (half_dx * [1,-1,1])
    result[:-1,:-1,1:] = grid.position + (half_dx * [-1,-1,1])
    result[1:,1:,:-1] = grid.position + (half_dx * [1,1,-1])
    result[:-1,1:,:-1] = grid.position + (half_dx * [-1,1,-1])
    result[1:,:-1,:-1] = grid.position + (half_dx * [1,-1,-1])
    result[:-1,:-1,:-1] = grid.position + (half_dx * [-1,-1,-1])
    
    return result
    
@grids.AbstractGrid.function_for_set
def connectivity(grid):
    cellcenters = grid.position
    shape = numpy.asarray(cellcenters.shape)
    shape[-1] = 8
    shape_with_boundary = numpy.asarray(cellcenters.shape) + 1
    shape_with_boundary = shape_with_boundary[:3]
    indices = numpy.arange(0, numpy.prod(shape_with_boundary), dtype = numpy.int).reshape(shape_with_boundary)
    result = numpy.zeros(shape, dtype = numpy.int)

    result[...,...,...,0] = indices[ :-1, :-1, :-1]
    result[...,...,...,1] = indices[1:  , :-1, :-1]
    result[...,...,...,2] = indices[ :-1,1:  , :-1]
    result[...,...,...,3] = indices[1:  ,1:  , :-1]
    
    result[...,...,...,4] = indices[ :-1, :-1,1:  ]
    result[...,...,...,5] = indices[1:  , :-1,1:  ]
    result[...,...,...,6] = indices[ :-1,1:  ,1:  ]
    result[...,...,...,7] = indices[1:  ,1:  ,1:  ]
    return result


