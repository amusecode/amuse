from amuse.support.data import values
from amuse.support.data.values import Quantity, new_quantity, zero
from amuse.support.units import constants
from amuse.support.units import units
from amuse.support.core import CompositeDictionary

from amuse.support.data.base import *
from amuse.support.data.memory_storage import *
from amuse.support.data import indexing

class AbstractGrid(AbstractSet):
    
    GLOBAL_DERIVED_ATTRIBUTES = {}
    


class Grid(AbstractGrid):
    def __init__(self, number_of_points_in_x_direction, number_of_points_in_y_direction, number_of_points_in_z_direction, storage = None):
        AbstractGrid.__init__(self)
        
        if storage is None:
            self._private.attribute_storage = InMemoryGridAttributeStorage(
                number_of_points_in_x_direction,
                number_of_points_in_y_direction,
                number_of_points_in_z_direction
            )
        else:
            self._private.attribute_storage = storage

    def _get_values(self, indices, attributes):
        result = self._private.attribute_storage._get_values(indices, attributes)
        return result
        
    def _set_values(self, indices, attributes, values):
        self._private.attribute_storage._set_values(indices, attributes, values)

    def _get_attribute_names(self):
        return self._private.attribute_storage._get_attribute_names()

    def _original_set(self):
        return self
        
    def _get_keys(self):
        return self._private.attribute_storage._get_keys()
    
    def __getitem__(self, index):
        return SubGrid(self, index)
        
class GridPoint(object):

    def __init__(self, index, grid):
        self.index = index
        self.grid = grid
    
    
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.grid._set_value_of_attribute(self.index, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise AttributeError("Can only assign quantities or other particles to an attribute.")
    
    def __getattr__(self, name_of_the_attribute):
        return self.particles_set._get_value_of_attribute(self.index, name_of_the_attribute)
                
class SubGrid(AbstractGrid):
    def __init__(self, grid, indices):
        AbstractParticleSet.__init__(self, particles)
        
        self._private.grid = grid
        self._private.indices = indices
    
    
    def _original_set(self):
        return self._private.grid
        
    def __getitem__(self, index):
        return SubGrid(self, indexing.combine_indices(self._private.indices, index))
        
    def _get_values(self, indices, attributes):
        result = self._private.attribute_storage._get_values(indices, attributes)
        return result
    
    def _set_values(self, indices, attributes, values):
        self._private.attribute_storage._set_values(indices, attributes, values)
            
    def _get_keys(self):
        return self._private_indices
