from amuse.support.core import CompositeDictionary, late
from amuse.units import constants
from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import quantities
from amuse.units.quantities import Quantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import zero
from amuse.support import exceptions
from amuse.datamodel.base import *
from amuse.datamodel.memory_storage import *
import numpy

from amuse.datamodel import indexing

class AbstractGrid(AbstractSet):
    
    GLOBAL_DERIVED_ATTRIBUTES = {}
        
    def _get_value_of_attribute(self, particle, index, attribute):
        if attribute in self._derived_attributes:
            return self._derived_attributes[attribute].get_value_for_entity(self, particle, index)
        else:
            return self._convert_to_entities_or_quantities(self.get_values_in_store(index, [attribute])[0])
            
    def _set_value_of_attribute(self, key, attribute, value):             
        if attribute in self._derived_attributes: 	 	 
            return self._derived_attributes[attribute].set_value_for_entity(self, key, value) 	 	 
        else:
            return self.set_values_in_store(key, [attribute], [value])
            
    def _get_values_for_entity(self, key, attributes):
        return self.get_values_in_store(key, attributes)
        
    def _set_values_for_entity(self, key, attributes, values):
        return self.set_values_in_store(key, attributes, values)
    
    def _get_particle(self, index):
        return GridPoint(index, self._original_set())
    
    
    def previous_state(self):
        return self._private.previous
        
        
    def savepoint(self, timestamp=None, **attributes):
        instance = type(self)()
        instance._private.attribute_storage = self._private.attribute_storage.copy()
        instance.collection_attributes.timestamp = timestamp
        
        for name, value in attributes.iteritems():
            setattr(instance.collection_attributes, name, value)
            
        instance._private.previous = self._private.previous
        self._private.previous = instance
        return instance
    
    
    def get_timestamp(self):
        return self.collection_attributes.timestamp
        
    def new_channel_to(self, other):
        return GridInformationChannel(self, other)
    
    def copy(self, memento = None, keep_structure = False):
        attributes = self.get_attribute_names_defined_in_store()
        values = self.get_values_in_store(None, attributes)
        result = self._factory_for_new_collection()(*self.shape)
        
        if memento is None:
            memento = {}
        memento[id(self._original_set())] = result
        
        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy(memento, keep_structure))
            else:
                converted.append(x)
        result.set_values_in_store(None, attributes, converted)
        
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        return result
        
    
    def _factory_for_new_collection(self):
        return Grid
        
    def empty_copy(self):
        result = self._factory_for_new_collection()(*self.shape)
        result.set_values_in_store(None, [],[])
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        return result
        
    def samplePoint(self, position, must_return_values_on_cell_center = False):
        if must_return_values_on_cell_center:
            return SamplePointOnCellCenter(self, position)
        else:
            return SamplePointWithIntepolation(self, position)
        
    def samplePoints(self, positions, must_return_values_on_cell_center = False):
        if must_return_values_on_cell_center:
            return SamplePointsOnGrid(self, positions, SamplePointOnCellCenter)
        else:
            return SamplePointsOnGrid(self, positions, SamplePointWithIntepolation)
    
    def __len__(self):
        return self.shape[0]
    
    def __iter__(self):
        for i in range(self.shape[0]):
            yield self[i]
            
    def get_all_indices_in_store(self):
        return self.get_all_keys_in_store()
        
    def can_extend_attributes(self):
        return self._original_set().can_extend_attributes()
        
    def __str__(self):
        dimensionstr = ' x '.join(([str(x) for x in self.shape]))
        
        return "{0} ({1})".format(
            self.__class__.__name__, 
            dimensionstr
        )
        
class Grid(AbstractGrid):
    DEFAULT_AXES_NAMES = ('x', 'y', 'z')
    
    
    def __init__(self, *args, **kwargs):
        AbstractGrid.__init__(self)
        
        if "storage" in kwargs:
            self._private.attribute_storage = kwargs['storage']
        else:
            self._private.attribute_storage = InMemoryGridAttributeStorage(*args)
        
        
        if "axes_names" in kwargs:
            self._private.axes_names = kwargs['axes_names']
        else:
            self._private.axes_names = self.DEFAULT_AXES_NAMES
            
        self._private.previous = None
        self.collection_attributes.timestamp = None
        self.add_vector_attribute("position", self._private.axes_names[0:len(self.shape)])
        
    def can_extend_attributes(self):
        return self._private.attribute_storage.can_extend_attributes()
    
    @classmethod
    def create(cls, shape, lengths, axes_names = DEFAULT_AXES_NAMES):
        """Returns a grid with cells between 0 and lengths.
        """
        result = cls(*shape, axes_names = axes_names)
    
        all_indices = numpy.indices(shape)+0.5
        
        def positions(indices, length, n):
            return length * (indices/n)
    
        for indices, length, n, axis_name in zip(all_indices, lengths, shape, axes_names):
            setattr(result, axis_name, positions(indices, length, n))
       
        return result
            
    def get_values_in_store(self, indices, attributes, by_key = True):
        result = self._private.attribute_storage.get_values_in_store(indices, attributes)
        return result
        
    def set_values_in_store(self, indices, attributes, values, by_key = True):
        self._private.attribute_storage.set_values_in_store(indices, attributes, values)

    def get_attribute_names_defined_in_store(self):
        return self._private.attribute_storage.get_defined_attribute_names()
        
    def _get_writeable_attribute_names(self):
        return self._private.attribute_storage._get_writeable_attribute_names()


    def _original_set(self):
        return self
        
    def get_all_keys_in_store(self):
        return self._private.attribute_storage.get_all_keys_in_store()
    
    def __getitem__(self, index):
        if indexing.number_of_dimensions_after_index(self.number_of_dimensions(), index) == 0:
            return GridPoint(index, self)
        else:
            return SubGrid(self, index)
    
    def _get_gridpoint(self, index):
        return GridPoint(index, self)
        
    def number_of_dimensions(self):
        return len(self.shape)
        
    @property
    def shape(self):
        return self._private.attribute_storage.storage_shape()
        
    @property
    def size(self):
        return numpy.prod(self.shape)
        
        
    def indices(self):
        return numpy.indices(self.shape)
    
    def iter_history(self):
        current = self._private.previous
        while not current is None:
            yield current
            current = current._private.previous
    
    @property
    def history(self):
        return reversed(list(self.iter_history()))
        
class SubGrid(AbstractGrid):
    def __init__(self, grid, indices):
        AbstractGrid.__init__(self, grid)
        
        self._private.grid = grid
        self._private.indices = indices
    
    def _original_set(self):
        return self._private.grid
        
    def get_values_in_store(self, indices, attributes, by_key = True):
        combined_index = indexing.combine_indices(self._private.indices, indices)
        result = self._private.grid.get_values_in_store(combined_index, attributes)
        return result
    
    def set_values_in_store(self, indices, attributes, values, by_key = True):
        combined_index = indexing.combine_indices(self._private.indices, indices)
        self._private.grid.set_values_in_store(combined_index, attributes, values)
            
    def get_all_keys_in_store(self):
        return None

    def number_of_dimensions(self):
        return indexing.number_of_dimensions_after_index(self._original_set().number_of_dimensions(), self._private.indices)
        
    def __getitem__(self, index):
        combined_index = indexing.combine_indices(self._private.indices, index)
        if indexing.number_of_dimensions_after_index(self._original_set().number_of_dimensions(), combined_index) == 0:
            return GridPoint(combined_index, self._original_set())
        else:
            return SubGrid(self._original_set(), combined_index)
            
    
    def get_attribute_names_defined_in_store(self):
        return self._private.grid.get_attribute_names_defined_in_store()
        
    def _get_writeable_attribute_names(self):
        return self._private.grid._get_writeable_attribute_names()
    
    @property
    def shape(self):
        return indexing.shape_after_index(self._private.grid.shape, self._private.indices )

    def indices(self):
        return [x[self._private.indices] for x in self._original_set().indices()]   
    
class GridPoint(object):

    def __init__(self, index, grid):
        object.__setattr__(self,"index",index)
        object.__setattr__(self,"grid",grid)    
    
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        try:
            self.grid._set_value_of_attribute(self.index, name_of_the_attribute, new_value_for_the_attribute)
        except Exception as ex:
            raise
            raise AttributeError("Could not assign to attribute {0}.".format(name_of_the_attribute))
    
    def __getattr__(self, name_of_the_attribute):
        return self.grid._get_value_of_attribute(self, self.index, name_of_the_attribute)
        
    def __eq__(self, other):
        return isinstance(other, type(self)) and other.index == self.index and other.grid == self.grid

    def __ne__(self, other):
        return not(isinstance(other, type(self)) and other.index == self.index and other.grid == self.grid)
    
    def get_containing_set(self):
        return self.grid
        
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
        values = self.source.get_values_in_store(self.index, attributes)
        
        
        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy_with_link_transfer(self.source, self.target))
            else:
                converted.append(x)
                
        self.target.set_values_in_store(self.index, attributes, converted)
    
        
    def copy(self):
        if not self.target.can_extend_attributes():
            self.copy_overlapping_attributes()
        else:
            self.copy_all_attributes()
        
    def copy_all_attributes(self):
        names_to_copy = self.source.get_attribute_names_defined_in_store()
        self.copy_attributes(list(names_to_copy))    

    def copy_overlapping_attributes(self):

        from_names = self.source.get_attribute_names_defined_in_store()
        to_names = self.target._get_writeable_attribute_names()
        names_to_copy = set(from_names).intersection(set(to_names))
             
        self.copy_attributes(list(names_to_copy))  

class SamplePointOnCellCenter(object):
    def __init__(self, grid, point):
        self.grid = grid
        self.point = point
        
    @late
    def position(self):
        return self.cell.position
    
    @late
    def index(self):
        offset = self.point - self.grid.get_minimum_position()
        indices = (offset / self.grid.cellsize())
        return numpy.floor(indices).astype(numpy.int)
        
    @late
    def isvalid(self):
        return numpy.logical_and(
            numpy.all(self.index >= self.grid.get_minimum_index()),
            numpy.all(self.index <= self.grid.get_maximum_index())
        )
    
    @late
    def cell(self):
        return self.grid[tuple(self.index)]
    
    def get_value_of_attribute(self, name_of_the_attribute):
        return getattr(self.cell, name_of_the_attribute)
    
    def __getattr__(self, name_of_the_attribute):
        return self.get_value_of_attribute(name_of_the_attribute)

class SamplePointWithIntepolation(object):
    """
    Vxyz =
    V000 (1 - x) (1 - y) (1 - z) +
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
    def position(self):
        return self.point
    
    @late
    def index(self):
        offset = self.point - self.grid.get_minimum_position()
        indices = (offset / self.grid.cellsize())
        return numpy.floor(indices)
        
    @late
    def index_for_000_cell(self):
        offset = self.point - self.grid[0,0,0].position
        indices = (offset / self.grid.cellsize())
        return numpy.floor(indices).astype(numpy.int)

    @late
    def index_for_111_cell(self):
        return self.index_for_000_cell + [1,1,1]
    
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
        
        result = numpy.asarray([
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
    
    @late
    def isvalid(self):
        return numpy.logical_and(
            numpy.all(self.index_for_000_cell >= self.grid.get_minimum_index()),
            numpy.all(self.index_for_111_cell <= self.grid.get_maximum_index())
        )
        
    def get_values_of_attribute(self, name_of_the_attribute):
        result = quantities.AdaptingVectorQuantity()
        for x in self.surrounding_cells:
            result.append(getattr(x, name_of_the_attribute))
        return result
    
    def __getattr__(self, name_of_the_attribute):
        values = self.get_values_of_attribute(name_of_the_attribute)
        return (values * self.weighing_factors).sum()
        
        

class SamplePointsOnGrid(object):
    
    def __init__(self, grid, points, samples_factory = SamplePointWithIntepolation):
        self.grid = grid
        self.samples = [samples_factory(grid, x) for x in points]
        self.samples = [x for x in self.samples if x.isvalid ]
        
    @late
    def indices(self):
        for x in self.samples:
            yield x.index
            
    @late
    def positions(self):
        for x in self.samples:
            yield x.position
        
    def __getattr__(self, name_of_the_attribute):
        result = quantities.AdaptingVectorQuantity()
        for x in self.samples:
            result.append(getattr(x, name_of_the_attribute))
        return result
    
    def __iter__(self):
        for x in len(self):
            yield self[x]
    
    def __getitem__(self, index):
        return self.samples[index]
    
    def __len__(self):
        return len(self.samples)

class SamplePointsOnMultipleGrids(object):
    
    def __init__(self, grids, points, samples_factory = SamplePointWithIntepolation, index_factory = None):
        self.grids = grids
        self.points = points
        self.samples_factory = samples_factory
        if index_factory is None:
            self.index = None
        else:
            self.index = index_factory(self.grids)
        
    def _grid_for_point(self, point):
        if self.index is None:
            for grid in self.grids:
                if (numpy.all(point >= grid.get_minimum_position()) and
                     numpy.all(point < grid.get_maximum_position())):
                    return grid
            return None 
        else:
            return self.index.grid_for_point(point)
    
    def filterout_duplicate_indices(self):
        previous_grid = None
        previous_index = None  
        filteredout = [] 
        for x in self.samples:
            if x.grid is previous_grid and numpy.all(x.index == previous_index):
                pass
            else:
                previous_grid= x.grid
                previous_index = x.index
                filteredout.append(x)
        self.samples = filteredout
        
    def get_samples(self):
        result = []
        for x in self.points:
            grid = self._grid_for_point(x)
            if grid is None:
                continue
            sample = self.samples_factory(grid, x)
            if not sample.isvalid:
                continue
            result.append(sample)
        return result
        
    @late
    def samples(self):
        result = []
        for x in self.points:
            grid = self._grid_for_point(x)
            if grid is None:
                continue
            sample = self.samples_factory(grid, x)
            if not sample.isvalid:
                continue
            result.append(sample)
        return result
        
    @late
    def indices(self):
        for x in self.samples:
            yield x.index
            
    @late
    def positions(self):
        for x in self.samples:
            yield x.position
        
    def __getattr__(self, name_of_the_attribute):
        self.get_samples()
        result = quantities.AdaptingVectorQuantity()
        for x in self.samples:
            result.append(getattr(x, name_of_the_attribute))
        return result
    
    def __iter__(self):
        for x in len(self):
            yield self[x]
    
    def __getitem__(self, index):
        return self.samples[index]
    
    def __len__(self):
        return len(self.samples)
        
        

class NonOverlappingGridsIndexer(object):
        
    def __init__(self, grids):
        self.grids = grids
        self.setup_index()
    
    @late
    def minimum_position(self):
        result = self.grids[0].get_minimum_position()
        for x in self.grids[1:]:
            minimum = x.get_minimum_position()
            result = result.minimum(minimum)
        print result
        return result
        
    def setup_index(self):
        smallest_boxsize = None
        for x in self.grids:
            boxsize = x.get_maximum_position() - x.get_minimum_position()
            if smallest_boxsize is None:
                smallest_boxsize = boxsize
            else:
                smallest_boxsize = boxsize.minimum(smallest_boxsize)
            
        self.smallest_boxsize = smallest_boxsize
        max_index = [0,0,0]
        
        for x in self.grids:
            index = (x.get_maximum_position() / smallest_boxsize)
            index = numpy.floor(index).astype(numpy.int)
            max_index = numpy.where(index > max_index, index, max_index)
            
        self.grids_on_index = numpy.zeros(max_index, 'int')
        
        for index,x in enumerate(self.grids):
            bottom_left = x.get_minimum_position()
            index_of_grid = (bottom_left / smallest_boxsize)
            size = ((x.get_maximum_position() - x.get_minimum_position()) / smallest_boxsize)
            i,j,k = numpy.floor(index_of_grid).astype(numpy.int)
            ni,nj,nk = numpy.floor(size).astype(numpy.int)
            self.grids_on_index[i:i+ni,j:j+nj,k:k+nk] = index
        
        
    def grid_for_point(self, position):
        index = ((position - self.minimum_position) / self.smallest_boxsize)
        index = numpy.floor(index).astype(numpy.int)
        index_of_grid = self.grids_on_index[tuple(index)]
        return self.grids[index_of_grid]
        
    def grids_for_points(self, points):
        index = ((points - self.minimum_position) / self.smallest_boxsize)
        index = numpy.floor(index).astype(numpy.int)
        print indices
        index_of_grid = self.grids_on_index[tuple(index)]
        return self.grids[index_of_grid]
