from amuse.support.data import values
from amuse.support.data.values import Quantity, new_quantity, zero
from amuse.support.units import constants
from amuse.support.units import units
from amuse.support.units import generic_unit_system
from amuse.support.core import CompositeDictionary, late

from amuse.support import exceptions
from amuse.support.data.base import *
from amuse.support.data.memory_storage import *
from amuse.support.data import indexing
from amuse.support.data import values

class AbstractGrid(AbstractSet):
    
    GLOBAL_DERIVED_ATTRIBUTES = {}
    
    def _get_value_of_attribute(self, key, attribute):
        if attribute in self._derived_attributes:
            return self._derived_attributes[attribute].get_value_for_entity(self, key)
        else:
            return self._convert_to_entities_or_quantities(self._get_values(key, [attribute])[0])
        
    def _set_value_of_attribute(self, key, attribute, value):
        if attribute in self._derived_attributes:
            return self._derived_attributes[attribute].set_value_for_entity(self, key, value)
        else:
            return self._set_values(key, [attribute], value)

    def _get_values_for_entity(self, key, attributes):
        return self._get_values(key, attributes)
        
    def _set_values_for_entity(self, key, attributes, values):
        return self._set_values(key, attributes, values)
    
    def _get_particle(self, index):
        return GridPoint(index, self._original_set())
    
    
    def previous_state(self):
        return self._private.previous
        
        
    def savepoint(self, timestamp=None):
        instance = type(self)()
        instance._private.attribute_storage = self._private.attribute_storage.copy()
        instance._private.timestamp = timestamp
        instance._private.previous = self._private.previous
        self._private.previous = instance
        return instance
    
    
    def new_channel_to(self, other):
        return GridInformationChannel(self, other)
    
    
    def copy_to_memory(self):
        attributes = self._get_attribute_names()
        values = self._get_values(None, attributes)
        result = Grid(*self.shape)
        result._set_values(None, attributes, values)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
       
        return result
        
    def samplePoint(self, position, must_return_values_on_cell_center = False):
        if must_return_values_on_cell_center:
            return SamplePointOnCellCenter(self, position)
        else:
            return SamplePointWithIntepolation(self, position)
        
        
class Grid(AbstractGrid):
    def __init__(self, number_of_points_in_x_direction = 1, number_of_points_in_y_direction = 1, number_of_points_in_z_direction = 1, storage = None):
        AbstractGrid.__init__(self)
        
        if storage is None:
            self._private.attribute_storage = InMemoryGridAttributeStorage(
                number_of_points_in_x_direction,
                number_of_points_in_y_direction,
                number_of_points_in_z_direction
            )
        else:
            self._private.attribute_storage = storage
            
    @classmethod
    def create(cls, shape, lengths):
        """Returns a grid with cells between 0 and lengths.
        """
        result = cls(shape[0], shape[1], shape[2])
        ix,iy,iz=(numpy.indices(shape)+0.5)
        
        def positions(indices, length, n):
            return length * (indices/n)
    
        result.x = positions(ix, lengths[0], shape[0])
        result.y = positions(iy, lengths[1], shape[1])
        result.z = positions(iz, lengths[2], shape[2])
        
        return result
            
    def _get_values(self, indices, attributes):
        result = self._private.attribute_storage._get_values(indices, attributes)
        return result
        
    def _set_values(self, indices, attributes, values):
        self._private.attribute_storage._set_values(indices, attributes, values)

    def _get_attribute_names(self):
        return self._private.attribute_storage._get_attribute_names()
        
    def _get_writeable_attribute_names(self):
        return self._private.attribute_storage._get_writeable_attribute_names()


    def _original_set(self):
        return self
        
    def _get_keys(self):
        return self._private.attribute_storage._get_keys()
    
    def __getitem__(self, index):
        if indexing.number_of_dimensions_after_index(3, index) == 0:
            return GridPoint(index, self)
        else:
            return SubGrid(self, index)
            
    def number_of_dimensions(self):
        return 3
        
    @property
    def shape(self):
        return self._private.attribute_storage.storage_shape()
        
    @property
    def size(self):
        return numpy.prod(self.shape)
        
    def indices(self):
        return numpy.indices(self.shape)
    
    
class SubGrid(AbstractGrid):
    def __init__(self, grid, indices):
        AbstractGrid.__init__(self, grid)
        
        self._private.grid = grid
        self._private.indices = indices
    
    def _original_set(self):
        return self._private.grid
        
    def _get_values(self, indices, attributes):
        combined_index = indexing.combine_indices(self._private.indices, indices)
        result = self._private.grid._get_values(combined_index, attributes)
        return result
    
    def _set_values(self, indices, attributes, values):
        combined_index = indexing.combine_indices(self._private.indices, indices)
        self._private.grid._set_values(combined_index, attributes, values)
            
    def _get_keys(self):
        return None

    def number_of_dimensions(self):
        return indexing.number_of_dimensions_after_index(3, self._private.indices)
        
    def __getitem__(self, index):
        combined_index = indexing.combine_indices(self._private.indices, index)
        if indexing.number_of_dimensions_after_index(3, combined_index) == 0:
            return GridPoint(combined_index, self._original_set())
        else:
            return SubGrid(self._original_set(), combined_index)
            
    
    def _get_attribute_names(self):
        return self._private.grid._get_attribute_names()
        
    def _get_writeable_attribute_names(self):
        return self._private.grid._get_writeable_attribute_names()
    
        


    def indices(self):
        return [x[self._private.indices] for x in self._original_set().indices()]   
    
class GridPoint(object):

    def __init__(self, index, grid):
        object.__setattr__(self,"index",index)
        object.__setattr__(self,"grid",grid)    
    
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.grid._set_value_of_attribute(self.index, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise AttributeError("Can only assign quantities or other particles to an attribute.")
    
    def __getattr__(self, name_of_the_attribute):
        return self.grid._get_value_of_attribute(self.index, name_of_the_attribute)
        

class GridInformationChannel(object):
    """
    A channel to copy attributes from one grid to another.
    For each dimension copies cells from 0 - min(grid0.size, grid1.size).
    """
    
    def __init__(self, source, target):
        self.source = source
        self.target = target
        self._reindex()
        
    def _reindex(self):
        source_shape = self.source.shape
        target_shape = self.target.shape
        if len(source_shape) != len(target_shape):
            raise exceptions.AmuseException("The source and target grids do not have the same dimensions, cannot use this channel")
        index = [numpy.s_[0:min(x,y)] for x,y in zip(source_shape, target_shape)]
        index = tuple(index)
        
        self.index = index
        
        
    def copy_attributes(self, attributes):
        data = self.source._get_values(self.index, attributes)
        self.target._set_values(self.index, attributes, data)
    
    def copy(self):
        self.copy_attributes(self.target._get_writeable_attribute_names())
    
class SamplePointWithIntepolation(object):
    """
    Vxyz =	 V000 (1 - x) (1 - y) (1 - z) +
V100 x (1 - y) (1 - z) + 
V010 (1 - x) y (1 - z) + 
V001 (1 - x) (1 - y) z +
V101 x (1 - y) z + 
V011 (1 - x) y z + 
V110 x y (1 - z) + 
V111 x y z
    """
    
    def __init__(self, grid, point):
        self.grid = grid
        self.point = point
        
    @late
    def index(self):
        offset = self.point - self.grid.get_minimum_position()
        indices = (offset / self.grid.cellsize()).value_in(units.none)
        return numpy.floor(indices)
        
    @late
    def index_for_000_cell(self):
        offset = self.point - self.grid[0,0,0].position
        indices = (offset / self.grid.cellsize()).value_in(units.none)
        return numpy.floor(indices).astype(numpy.int)

    @late
    def surrounding_cell_indices(self):
        cell000 = self.index_for_000_cell
        translations = [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 0],
            [1, 1, 1],
        ]
        return cell000 + translations
        
    @late
    def weighing_factors(self):
        x0,y0,z0 = self.grid[tuple(self.index_for_000_cell)].position
        x1,y1,z1 = self.grid[tuple(self.index_for_000_cell + [1,1,1])].position
        x,y,z = self.point
        
        
        dx1 = (x1 - x) / (x1 - x0)
        dy1 = (y1 - y) / (y1 - y0)
        dz1 = (z1 - z) / (z1 - z0)
        dx0 = (x - x0) / (x1 - x0)
        dy0 = (y - y0) / (y1 - y0)
        dz0 = (z - z0) / (z1 - z0)
        
        result = values.AdaptingVectorQuantity()
        result.extend([
            dx1 * dy1 * dz1,
            dx0 * dy1 * dz1,
            dx1 * dy0 * dz1,
            dx1 * dy1 * dz0,
            dx0 * dy1 * dz0,
            dx1 * dy0 * dz0,
            dx0 * dy0 * dz1,
            dx0 * dy0 * dz0
        ] )
        return result 
    
    @late
    def surrounding_cells(self):
        return [self.grid[tuple(x)] for x in self.surrounding_cell_indices]
        
    def get_values_of_attribute(self, name_of_the_attribute):
        result = values.AdaptingVectorQuantity()
        for x in self.surrounding_cells:
            result.append(getattr(x, name_of_the_attribute))
        return result
    
    def __getattr__(self, name_of_the_attribute):
        values = self.get_values_of_attribute(name_of_the_attribute)
        print values, self.weighing_factors
        return (values * self.weighing_factors).sum()
        
        
        
