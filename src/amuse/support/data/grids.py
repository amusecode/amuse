from amuse.support.data import values
from amuse.support.data.values import Quantity, new_quantity, zero
from amuse.support.units import constants
from amuse.support.units import units
from amuse.support.core import CompositeDictionary

from amuse.support.data.base import *
from amuse.support.data.memory_storage import *

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
        return SubGrid(sel, index)
        
        
class SubGrid(AbstractGrid):
    def __init__(self, grid, indices):
        AbstractParticleSet.__init__(self, particles)
        
        self._private.grid = grid
        self._private.indices = indices
    
    
    def _original_set(self):
        return self._private.grid
        
    def __getitem__(self, index):
        pass
        
    def _get_values(self, indices, attributes):
        result = self._private.attribute_storage._get_values(indices, attributes)
        return result
    
    def _set_values(self, indices, attributes, values):
        self._private.attribute_storage._set_values(indices, attributes, values)
            
    def _get_keys(self):
        return self._private_indices
