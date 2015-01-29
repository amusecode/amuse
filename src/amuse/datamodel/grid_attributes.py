from amuse.units.quantities import zero

import numpy

from amuse.datamodel import base
from amuse.datamodel import grids
grids.AbstractGrid.add_global_vector_attribute("position", ["x","y","z"])
grids.AbstractGrid.add_global_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])
grids.AbstractGrid.add_global_vector_attribute("magnetic_field", ["B1i","B2i","B3i"])


@grids.AbstractGrid.caching_function_for_set
def cellsize(grid):
    """Returns the lenght of each direction in the grid.
    Works for regular cartesion grids.
    """
    result = grid[tuple(grid.get_minimum_index())].position*0.
    Ndim=len(grid.shape)
    cell1 = grid[(0,)*Ndim]
    for i in range(len(result)):
        if grid.shape[i] > 1:
            cell2=grid[(0,)*i+(1,)+(0,)*(Ndim-1-i)]
            result[i:i+1]=(cell2.position-cell1.position)[i]      
    return result

@grids.AbstractGrid.caching_function_for_set
def get_minimum_index(grid):
    return numpy.zeros_like(grid.shape)
    
@grids.AbstractGrid.caching_function_for_set
def get_maximum_index(grid):
    return grid.shape - numpy.ones_like(grid.shape)
    
@grids.AbstractGrid.caching_function_for_set
def get_minimum_position(grid):
    return grid[tuple(grid.get_minimum_index())].position - 0.5 * grid.cellsize()
    
@grids.AbstractGrid.caching_function_for_set
def get_maximum_position(grid):
    return grid[tuple(grid.get_maximum_index())].position + 0.5 * grid.cellsize()
    
@grids.AbstractGrid.caching_function_for_set
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

@grids.AbstractGrid.function_for_set
def overlaps(grid, grid1,eps=None):
    """simple test for overlap
       optional keyword parameter:
       
       [eps]: size of buffer (to ignore just touching regions)
    """
    minp=grid.get_minimum_position()
    maxp=grid.get_maximum_position()
    minp1=grid1.get_minimum_position()
    maxp1=grid1.get_maximum_position()
    if eps is not None:
      minp+=eps
      maxp-=eps
      minp1+=eps
      maxp1-=eps
    if (maxp<=minp1).sum()>0 or (minp>=maxp1).sum()>0:
      return False
    return True
    
@grids.AbstractGrid.function_for_set
def get_overlap_with(grid, grid1,eps=None):
    """return overlapping subgrid"""
    if not grid.overlaps(grid1,eps):
      return None
    minindex=grid.get_minimum_index()
    maxindex=grid.get_maximum_index()
    cellsize=grid.cellsize()
        
    minp=grid.get_minimum_position()
    maxp=grid.get_maximum_position()
    minp1=grid1.get_minimum_position()
    maxp1=grid1.get_maximum_position()
    if eps is not None:
      minp1+=eps
      maxp1-=eps
    
    index_of_minp1=numpy.maximum( numpy.array((minp1-minp)/cellsize,'int'), minindex[:len(cellsize)])
    index_of_maxp1=numpy.minimum( numpy.array((maxp1-minp)/cellsize,'int'), maxindex[:len(cellsize)])
    slices=()
    for i,j in zip(index_of_minp1,index_of_maxp1):
      slices+=(slice(i,j+1),)
    if len(slices)!= len(minindex): slices+=(Ellipsis,)
    return grid[slices]


#@grids.AbstractGrid.function_for_set
#def select_fully_inside(grid, cellsizes=(), coordinates=()):
#    """returns boolean array with cells with cellsizes centered on 
#    coordinates fully inside the grid """
#    gridminx,gridminy=sys.grid.get_minimum_position()
#    gridmaxx,gridmaxy=sys.grid.get_maximum_position()

