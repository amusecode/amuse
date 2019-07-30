from amuse.units.quantities import zero, as_vector_quantity, column_stack

import numpy

from amuse.datamodel import base
from amuse.datamodel import grids

# maintain for backwards compatibility, these should go..
grids.Grid.add_global_vector_attribute("position", ["x","y","z"])
grids.Grid.add_global_vector_attribute("momentum", ["rhovx","rhovy","rhovz"])
grids.Grid.add_global_vector_attribute("magnetic_field", ["B1i","B2i","B3i"])

@grids.BaseGrid.caching_function_for_set
def cellsize(grid):
    raise Exception("a {0} does not have a constant cellsize, use the cellsizes method instead".format(grid.__class__.__name__))

@grids.RegularBaseGrid.caching_function_for_set
def cellsize(grid):
    """Returns the lenght of each direction in the grid.
    Works for regular and cartesian grids.
    """
    result = []
    Ndim= len(grid.shape)
    cell1 = grid[(0,)*Ndim]
    for i in range(len( grid[grid.get_minimum_index()].position )): # shape of position not necessarily same as shape of the grid?
        if grid.shape[i] > 1:
            cell2=grid[(0,)*i+(1,)+(0,)*(Ndim-1-i)]
            result.append((cell2.position-cell1.position)[i])      
    return as_vector_quantity(result)

@grids.BaseGrid.caching_function_for_set
def get_minimum_index(grid):
    raise Exception("not implemented")

@grids.StructuredBaseGrid.caching_function_for_set
def get_minimum_index(grid):
    return tuple(numpy.zeros_like(grid.shape))
    
@grids.BaseGrid.caching_function_for_set
def get_maximum_index(grid):
    raise Exception("not implemented")

@grids.StructuredBaseGrid.caching_function_for_set
def get_maximum_index(grid):
    return tuple(grid.shape - numpy.ones_like(grid.shape))

@grids.BaseGrid.caching_function_for_set
def get_minimum_position(grid):
    raise Exception("not implemented")

@grids.RectilinearBaseGrid.caching_function_for_set
def get_minimum_position(grid):
    return grid[grid.get_minimum_index()].position - 0.5 * grid.cellsize()

@grids.BaseGrid.caching_function_for_set
def get_maximum_position(grid):
    raise Exception("not implemented")

@grids.RectilinearBaseGrid.caching_function_for_set
def get_maximum_position(grid):
    return grid[grid.get_maximum_index()].position + 0.5 * grid.cellsize()

@grids.BaseGrid.caching_function_for_set
def get_volume(grid):
    raise Exception("not implemented")
    
@grids.RectilinearBaseGrid.caching_function_for_set
def get_volume(grid):
    maximum_position = grid.get_maximum_position()
    minimum_position = grid.get_minimum_position()
    delta = maximum_position - minimum_position
    return delta.prod()
    
@grids.BaseGrid.function_for_set
def contains(grid, points):
    raise Exception("not implemented")

@grids.RectilinearBaseGrid.function_for_set
def contains(grid, points):
    return numpy.logical_and(
        numpy.all(points >= grid.get_minimum_position(), axis=len(points.shape)-1),
        numpy.all(points < grid.get_maximum_position(), axis=len(points.shape)-1)
    )
   
@grids.BaseGrid.function_for_set
def points(grid):
    raise Exception("not implemented")

@grids.RegularBaseGrid.function_for_set
def points(grid):
    shape=grid.shape
    dx = grid.cellsize()/2
    cell_centers=grid.position

    shape_with_boundary = numpy.asarray(cell_centers.shape) + 1
    shape_with_boundary[-1] -= 1
    result = numpy.zeros(shape_with_boundary)* cell_centers.flat[0]
    
    for i in range(2**len(shape)):
        slicing=[]
        offset=[]
        for j in range(len(shape)):
            if i & 2**j:
                slicing.append(slice(1,None)) 
                offset.append(1)
            else:
                slicing.append(slice(None,-1))
                offset.append(-1)
                        
        result[slicing]=cell_centers+dx*numpy.asarray(offset)    

    return result
    
@grids.BaseGrid.function_for_set
def connectivity(grid):
    raise Exception("not implemented")

@grids.RegularBaseGrid.function_for_set
def connectivity(grid):
    cellcenters = grid.position
    shape = numpy.asarray(cellcenters.shape)
    dim=len(shape)-1
    shape[-1] = 2**dim
    shape_with_boundary = numpy.asarray(cellcenters.shape) + 1
    shape_with_boundary = shape_with_boundary[:dim]
    indices = numpy.arange(0, numpy.prod(shape_with_boundary), dtype = numpy.int).reshape(shape_with_boundary)
    result = numpy.zeros(shape, dtype = numpy.int)

    for i in range(2**dim):
        slicing1=[]
        slicing2=[]
        for j in range(dim):
            if i & 2**j:
                slicing2.append(slice(1,None)) 
                if len(slicing1) == 0 or slicing1[-1] is not Ellipsis:
                    slicing1.append(Ellipsis)
            else:
                slicing2.append(slice(None,-1))
                if len(slicing1) == 0 or slicing1[-1] is not Ellipsis:
                    slicing1.append(Ellipsis)

        slicing1.append(i)

        result[slicing1]=indices[slicing2]    

    return result

@grids.BaseGrid.function_for_set
def overlaps(grid, grid1,eps=None):
    raise Exception("not implemented")

@grids.RectilinearBaseGrid.function_for_set
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
    
@grids.BaseGrid.function_for_set
def get_overlap_with(grid, grid1,eps=None):
    raise Exception("not implemented")

@grids.RegularBaseGrid.function_for_set
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

@grids.BaseGrid.function_for_set
def get_index(grid, pos=None, **kwargs):
    raise Exception("not implemented for a {0} grid".format(grid.__class__.__name__))

@grids.RegularBaseGrid.function_for_set
def get_index(grid, pos=None, **kwargs):
    pos=grid._get_array_of_positions_from_arguments(pos=pos,**kwargs)
    offset = pos - grid.get_minimum_position()
    indices = (offset / grid.cellsize())
    return numpy.floor(indices).astype(numpy.int)

@grids.BaseGrid.function_for_set
def _get_array_of_positions_from_arguments(grid, **kwargs):
    return grids._get_array_of_positions_from_arguments(grid.get_axes_names(), **kwargs)

           
      

