from amuse.support.core import CompositeDictionary, late
from amuse.units import constants
from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import quantities
from amuse.units.quantities import Quantity
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import zero
from amuse.units.quantities import column_stack
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
        try:
            instance = type(self)()
            instance._private.attribute_storage = self._private.attribute_storage.copy()
        except:
            instance=self.copy() # for the case of subgrid, maybe always ok

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
    def new_remapping_channel_to(self, other, remapper):
        return GridRemappingChannel(self, other, remapper)
    
    def copy(self, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        attributes = self.get_attribute_names_defined_in_store()
        attributes = [x for x in attributes if filter_attributes(self, x)]
        
        values = self.get_values_in_store(Ellipsis, attributes)
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
        result.set_values_in_store(Ellipsis, attributes, converted)
        
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        return result
        
    
    def _factory_for_new_collection(self):
        return self.__class__
        
    def empty_copy(self):
        result = self._factory_for_new_collection()(*self.shape)
        result.set_values_in_store(None, [],[])
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        return result
        
    def samplePoint(self, position=None, method="nearest", **kwargs):
        if method in ["nearest"]:
            return SamplePointOnCellCenter(self, position=position, **kwargs)
        elif method in ["interpolation", "linear"]:
            return SamplePointWithInterpolation(self, position=position, **kwargs)
        else:
            raise Exception("unknown sample method")
        
    def samplePoints(self, positions=None, method="nearest", **kwargs):
        if method in ["nearest"]:
            return SamplePointsOnGrid(self, positions, SamplePointOnCellCenter, **kwargs)
        elif method in ["interpolation", "linear"]:
            return SamplePointsOnGrid(self, positions, SamplePointWithInterpolation, **kwargs)
        else:
            raise Exception("unknown sample method")

    
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
        attributes=self.get_attribute_names_defined_in_store()
        settable=self.get_defined_settable_attribute_names()
        strings=[a if a in settable else a+" (ro)" for a in attributes]
        attrstr= ', '.join(strings)
        return "{0}({1}) ( {2} )".format(
            self.__class__.__name__, 
            dimensionstr,
            attrstr
        )

    def iter_history(self):
        raise Exception("not implemented")
        
    @property
    def history(self):
        return reversed(list(self.iter_history()))

    def get_timeline_of_attribute(self, attribute):
        timeline = []
        for x in self.history:
            timeline.append((x.collection_attributes.timestamp, getattr(x,attribute)))
        return timeline

    def get_timeline_of_attribute_as_vector(self, attribute):
        timestamps = AdaptingVectorQuantity()
        timeline = AdaptingVectorQuantity()
        for x in self.history:
            timestamps.append(x.collection_attributes.timestamp)
            timeline.append(getattr(x,attribute))
        return timestamps,timeline
        
class BaseGrid(AbstractGrid):
    def __init__(self, *args, **kwargs):
        AbstractGrid.__init__(self)
        
        if "storage" in kwargs:
            self._private.attribute_storage = kwargs['storage']
        else:
            self._private.attribute_storage = InMemoryGridAttributeStorage(*args)
        
        self._private.previous = None
        self.collection_attributes.timestamp = None
        
    def can_extend_attributes(self):
        return self._private.attribute_storage.can_extend_attributes()
                
    def get_values_in_store(self, indices, attributes, by_key = True):
        result = self._private.attribute_storage.get_values_in_store(indices, attributes)
        return result
        
    def set_values_in_store(self, indices, attributes, values, by_key = True):
        self._private.attribute_storage.set_values_in_store(indices, attributes, values)
        
    def set_values_in_store_async(self, indices, attributes, values, by_key = True):
        return self._private.attribute_storage.set_values_in_store_async(indices, attributes, values)

    def get_attribute_names_defined_in_store(self):
        return self._private.attribute_storage.get_defined_attribute_names()
        
    def get_defined_settable_attribute_names(self):
        return self._private.attribute_storage.get_defined_settable_attribute_names()


    def _original_set(self):
        return self
        
    def get_all_keys_in_store(self):
        return self._private.attribute_storage.get_all_keys_in_store()
    
    def __getitem__(self, index):
        return new_subgrid_from_index(self, index)
    
    def iter_cells(self):
        shape = numpy.asarray(self.shape)
        
        index = 0 * shape
        
        while index[0] < shape[0]:
            yield self._get_gridpoint(tuple(index))
            
            index[-1] += 1
            for i in range(len(self.shape) - 1, 0, -1):
                if index[i] >= shape[i]:
                    index[i] = 0
                    index[i-1] += 1
            
                
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
    
    @classmethod
    def create(cls,*args,**kwargs):
        print ("Grid.create deprecated, use new_regular_grid instead")
        return new_regular_grid(*args,**kwargs)

    def get_axes_names(self):
        if "position" in self.GLOBAL_DERIVED_ATTRIBUTES:
            result=self.GLOBAL_DERIVED_ATTRIBUTES["position"].attribute_names
        elif "position" in self._derived_attributes:
            result=self._derived_attributes["position"].attribute_names
        else:
            try:
              result=self._axes_names
            except:
              raise Exception("do not know how to find axes_names")
        return list(result)

class UnstructuredGrid(BaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(BaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class StructuredBaseGrid(BaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(BaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class StructuredGrid(StructuredBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(StructuredBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class RectilinearBaseGrid(StructuredBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(StructuredBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class RectilinearGrid(RectilinearBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(RectilinearBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class RegularBaseGrid(RectilinearBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(RectilinearBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class RegularGrid(RegularBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(RegularBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class CartesianBaseGrid(RegularBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(RegularBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
class CartesianGrid(CartesianBaseGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(CartesianBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)

# maintains compatibility with previous def.
class Grid(RegularGrid):
    GLOBAL_DERIVED_ATTRIBUTES=CompositeDictionary(RegularBaseGrid.GLOBAL_DERIVED_ATTRIBUTES)
    
def new_cartesian_grid(shape, cellsize, axes_names = "xyz",offset=None):
        """Returns a cartesian grid with cells of size cellsize.
        """
        if len(axes_names)<len(shape):
          raise Exception("provide enough axes names")

        result = CartesianGrid(*shape)
    
        all_indices = numpy.indices(shape)+0.5
        
        if offset is None:
            offset=[0.*cellsize]*len(shape)
        
        def positions(indices):
            return cellsize * indices
    
        for indices, n, axis_name, of in zip(all_indices, shape, axes_names,offset):
            setattr(result, axis_name, positions(indices)+of)
                 
        result.add_vector_attribute("position", axes_names[0:len(shape)])

        object.__setattr__(result,"_grid_type","cartesian") # for now for convenience, eventually to be converted in seperate classes
        object.__setattr__(result,"_cellsize",cellsize)
       
        return result

def new_regular_grid(shape, lengths, axes_names = "xyz",offset=None):
        """Returns a regular grid with cells between 0 and lengths.
        """
        if len(axes_names)<len(shape):
          raise Exception("provide enough axes names")
        if len(lengths)!=len(shape):
          raise Exception("shape and lengths do not conform")

        result = RegularGrid(*shape)
    
        all_indices = numpy.indices(shape)+0.5
        
        if offset is None:
            offset=[0.*l for l in lengths]
        
        def positions(indices, length, n):
            return length * (indices/n)
    
        for indices, length, n, axis_name, of in zip(all_indices, lengths, shape, axes_names,offset):
            setattr(result, axis_name, positions(indices, length, n)+of)
                 
        result.add_vector_attribute("position", axes_names[0:len(shape)])

        object.__setattr__(result,"_grid_type","regular")
        object.__setattr__(result,"_lengths",lengths)
       
        return result

def new_rectilinear_grid(shape, axes_cell_boundaries=None, cell_centers=None, axes_names = "xyz",offset=None):
        """Returns a rectilinear grid with cells at positions midway given cell boundaries.
        """
        if len(axes_names)<len(shape):
            raise Exception("provide enough axes names")
        if not (axes_cell_boundaries or cell_centers):
            raise Exception("provide cell boundaries or cell_centers")
        if axes_cell_boundaries and len(axes_cell_boundaries)!=len(shape):
            raise Exception("length of shape and axes positions do not conform")
        if axes_cell_boundaries:
            for s,b in zip(shape,axes_cell_boundaries):
                if len(b)!=s+1:
                    raise Exception("number of cell boundary arrays error (must be {0} instead of {1})".format(s+1,len(b)))
        if cell_centers and len(cell_centers)!=len(shape):
            raise Exception("length of shape and axes positions do not conform")
        if cell_centers:
            for s,b in zip(shape,cell_centers):
                if len(b)!=s:
                    raise Exception("number of cell_center arrays error (must be {0} instead of {1})".format(s+1,len(b)))

        result = RectilinearGrid(*shape)

        all_indices = numpy.indices(shape)
    
        #~ axes_cell_boundaries=[numpy.sort(b) for b in axes_cell_boundaries]
    
        if axes_cell_boundaries:
            positions=[(b[1:]+b[:-1])/2 for b in axes_cell_boundaries]
        if cell_centers:
            positions=cell_centers
        
        if offset is None:
            offset=[0.*l[0] for l in positions]
            
        for indices, axis_pos, axis_name, of in zip(all_indices, positions, axes_names, offset):
            setattr(result, axis_name, axis_pos[indices]+of)
       
        result.add_vector_attribute("position", axes_names[0:len(shape)])
        
        object.__setattr__(result,"_grid_type","rectilinear")
        object.__setattr__(result,"_axes_cell_boundaries",axes_cell_boundaries)
        object.__setattr__(result,"_cell_centers",cell_centers)
       
        return result

def new_structured_grid(shape, cell_corners, cell_positions=None, axes_names = "xyz", offset=None):
        """Returns a structured grid with cells with given corners and cell_positions.
           if not present, cell positions default to average of corner positions.
        """
        if len(axes_names)<len(shape):
            raise Exception("provide enough axes names")
        if len(cell_corners)!=len(shape):
            raise Exception("dimensions of shape and cell_boundaries do not conform")
        for c in cell_corners:
            if not numpy.all([s1==s2+1 for s1,s2 in zip(c.shape,shape)]):
                shape1=[s+1 for s in shape]
                raise Exception("size of cell_corner arrays must be {0} instead of {1}".format(shape1.__str__(),c.shape.__str__()))        

        if cell_positions is None:
              cell_positions=[]
              for cc in cell_corners:
                  cp=numpy.zeros(shape) * cc.flat[0]
                  for i in range(2**len(shape)):
                      slicing=[]
                      for j in range(len(shape)):
                          if i & 2**j:
                              slicing.append(slice(1,None)) 
                          else:
                              slicing.append(slice(None,-1))
                      cp=cp+cc[slicing]              
                  cell_positions.append(cp/2**len(shape))          

        if len(cell_positions)!=len(shape):
            raise Exception("dimensions of shape and cell_positions do not conform")
        for c in cell_positions:
            if not numpy.all([s1==s2 for s1,s2 in zip(c.shape,shape)]):
                raise Exception("size of cell_position arrays must be {0} instead of {1}".format(shape1.__str__(),c.shape.__str__()))        

        if offset is None:
            offset=[0.*l.flat[0] for l in cell_positions]

        result = StructuredGrid(*shape)

        for axis_name, pos, of in zip(axes_names, cell_positions, offset):
            setattr(result, axis_name, pos + of)

        result.add_vector_attribute("position", axes_names[0:len(shape)])
        
        object.__setattr__(result,"_grid_type","structured")
        object.__setattr__(result,"_cell_corners", cell_corners)
       
        return result

def new_unstructured_grid(size, num_corners, cell_corners, cell_positions=None, axes_names="xyz", offset=None):
        """Returns an unstructured grid with cells with given corners and cell_positions.
           if not present, cell positions default to average of corner positions.
        """
        dimensions = cell_corners.size / (num_corners * size)
        if len(axes_names)<dimensions:
            raise Exception("provide enough axes names")
        if len(cell_corners.shape) != 2:
            raise Exception("incorrect shape for cell_corners, the number of dimensions of the array should be exactly three (dimensions, size, corners)")
        if cell_corners.shape[0] != dimensions:
            raise Exception("incorrect shape for cell_corners, first dimension should equal the number of dimensions of the space in which the grid is defined")
        if cell_corners.shape[1] != size * num_corners:
            raise Exception("incorrect shape for cell_corners, second dimension should equal the grid size times the number of corners per cell")

        if cell_positions is None:
            cell_positions=[]
            for cc in cell_corners:
                c = cc.reshape(size, num_corners)
                cp=numpy.zeros(size)
                for i in range(size):
                    cp[i] = c[i].sum() / num_corners

                cell_positions.append(cp)

        cell_positions = numpy.array(cell_positions)

        if len(cell_positions.shape) != 2:
            raise Exception("incorrect shape for cell_positions, the number of dimensions of the array should be exactly two (dimensions, size)")
        if cell_positions.shape[0] != dimensions:
            raise Exception("dimensions of cell_positions and cell_corners do not conform")
        if cell_positions.shape[1] != size:
            raise Exception("size of cell_positions and size do not conform")

        if offset is None:
            offset=[0.*l.flat[0] for l in cell_positions]

        result = UnstructuredGrid(size)

        for axis_name, pos, of in zip(axes_names, cell_positions, offset):
            setattr(result, axis_name, pos + of)

        result.add_vector_attribute("position", axes_names[0:dimensions])
        
        object.__setattr__(result,"_grid_type","unstructured")
        object.__setattr__(result,"_num_corners", num_corners)
        object.__setattr__(result,"_cell_corners", cell_corners)
       
        return result

class SubGrid(AbstractGrid):
    def __init__(self, grid, indices):
        AbstractGrid.__init__(self, grid)
        
        self._private.previous=None
        self._private.grid = grid
        self._private.indices = indexing.normalize_slices(grid.shape,indices)
        self._private.collection_attributes=grid.collection_attributes
        
    def _original_set(self):
        return self._private.grid
    
    def previous_state(self):
        previous=self._private.previous
        if previous:
            return previous
        previous=self._private.grid.previous_state()
        if previous:
            return previous[self._private.indices]
        return previous
        
    def get_values_in_store(self, indices, attributes, by_key = True):
        normalized_indices = indexing.normalize_slices(self.shape,indices)
        combined_index = indexing.combine_indices(self._private.indices, normalized_indices)
        result = self._private.grid.get_values_in_store(combined_index, attributes)
        return result
    
    def set_values_in_store(self, indices, attributes, values, by_key = True):
        normalized_indices = indexing.normalize_slices(self.shape,indices)
        combined_index = indexing.combine_indices(self._private.indices, normalized_indices)
        self._private.grid.set_values_in_store(combined_index, attributes, values)
        
    def set_values_in_store_async(self, indices, attributes, values, by_key = True):
        normalized_indices = indexing.normalize_slices(self.shape,indices)
        combined_index = indexing.combine_indices(self._private.indices, normalized_indices)
        return self._private.grid.set_values_in_store_async(combined_index, attributes, values)
            
    def get_all_keys_in_store(self):
        return Ellipsis

    def number_of_dimensions(self):
        return indexing.number_of_dimensions_after_index(self._original_set().number_of_dimensions(), self._private.indices)
        
    def __getitem__(self, index):
        normalized_index= indexing.normalize_slices(self.shape,index)
        combined_index = indexing.combine_indices(self._private.indices, normalized_index)
        return new_subgrid_from_index(self._original_set(), combined_index)
                
    def get_attribute_names_defined_in_store(self):
        return self._private.grid.get_attribute_names_defined_in_store()
        
    def get_defined_settable_attribute_names(self):
        return self._private.grid.get_defined_settable_attribute_names()
    
    @property
    def shape(self):
        return indexing.shape_after_index(self._private.grid.shape, self._private.indices )

    def indices(self):
        return [x[self._private.indices] for x in self._original_set().indices()]
        
    def __eq__(self, other):
        if self._private.grid!=other._private.grid:
          return False
        else:
          if numpy.all(numpy.array(self.indices())==numpy.array(other.indices())):
            return True
          else:
            return False
        
    def __ne__(self,other):
        return not(self==other)

    def _factory_for_new_collection(self):
        return Grid

    def iter_history(self):
        if self._private.previous:
            current = self._private.previous
            while not current is None:
                yield current
                current = current._private.previous
            return

        current = self._original_set().previous_state()
        while not current is None:
            yield current[self._private.indices]
            current = current.previous_state()
        
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

    def iter_history(self):
        current = self.get_containing_set().previous_state()
        while not current is None:
            yield current[self.index]
            current = current.previous_state()

    @property
    def history(self):
        return reversed(list(self.iter_history()))

    def get_timeline_of_attribute(self, attribute):
        timeline = []
        for x in self.history:
            timeline.append((x.grid.collection_attributes.timestamp, getattr(x,attribute)))
        return timeline

    def get_timeline_of_attribute_as_vector(self, attribute):
        timestamps = AdaptingVectorQuantity()
        timeline = AdaptingVectorQuantity()
        for x in self.history:
            timestamps.append(x.grid.collection_attributes.timestamp)
            timeline.append(getattr(x,attribute))
        return timestamps,timeline


        
def new_subgrid_from_index(grid, index):
    if indexing.number_of_dimensions_after_index(grid.number_of_dimensions(), index) == 0:
        return GridPoint(index, grid)
    else:
        return SubGrid(grid, index)

class GridRemappingChannel(object):
    """
    A channel to remap attributes from one grid to another.
    """
    
    def __init__(self, source, target, remapper):
        self.source = source
        self.target = target
        if callable(remapper):
            self.remapper = remapper( source, target)
        else:
            self.remapper = remapper

    def get_overlapping_attributes(self):
        from_names = self.source.get_attribute_names_defined_in_store()
        to_names = self.target.get_defined_settable_attribute_names()
        names_to_copy = set(from_names).intersection(set(to_names))
        return list(names_to_copy)

    def copy_attributes(self, attributes):
        self.remapper.forward_mapping(attributes)
                
    def copy(self):
        if not self.target.can_extend_attributes():
            self.copy_overlapping_attributes()
        else:
            self.copy_all_attributes()
        
    def copy_all_attributes(self):
        names_to_copy = self.source.get_attribute_names_defined_in_store()
        self.copy_attributes(list(names_to_copy))    

    def copy_overlapping_attributes(self):
        names_to_copy = self.get_overlapping_attributes()
        self.copy_attributes(names_to_copy)  


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
        
    def get_values(self, attributes):
        values = self.source.get_values_in_store(self.index, attributes)
        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy_with_link_transfer(self.source, self.target))
            else:
                converted.append(x)
        return converted

    def get_overlapping_attributes(self):
        from_names = self.source.get_attribute_names_defined_in_store()
        to_names = self.target.get_defined_settable_attribute_names()
        names_to_copy = set(from_names).intersection(set(to_names))
        return list(names_to_copy)
    
    def copy_attributes(self, attributes):
        converted=self.get_values(attributes)        
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
        names_to_copy = self.get_overlapping_attributes()
        self.copy_attributes(names_to_copy)

    def transform_values(self, attributes, f):
        values = self.source.get_values_in_store(self.index, attributes)
        return f(*values)
        
    def transform(self, target, function, source):
        """ Copy and transform values of one attribute from the source set to the target set.

        :argument target: name of the attributes in the target set
        :argument function: function used for transform, should return tuple
        :argument source: name of the attribute in the source set

        >>> from amuse.datamodel import Grid
        >>> grid1 = Grid(2)
        >>> grid2 = Grid(2)
        >>> grid1.attribute1 = 1
        >>> grid1.attribute2 = 2
        >>> channel = grid1.new_channel_to(grid2)
        >>> channel.transform(["attribute3","attribute4"], lambda x,y: (y+x,y-x), ["attribute1","attribute2"])
        >>> print grid2.attribute3
        [3 3]
        >>> print grid2.attribute4
        [1 1]

        """
        if function is None:
            function=lambda *x : x
        
        if not self.target.can_extend_attributes():
            target_attributes = self.target.get_defined_settable_attribute_names()
            if not set(target).issubset(set(target_attributes)):
                raise Exception("trying to set unsettable attributes {0}".format(
                                list(set(target)-set(target_attributes))) )
        converted=self.transform_values(source, function)
        if len(converted) != len(target):
            raise Exception("function {0} returns {1} values while target attributes are {2} of length {3}".format(
                            function.__name__, len(converted), target, len(target)))
        self.target.set_values_in_store(self.index, target, converted)        

class SamplePointOnCellCenter(object):
    def __init__(self, grid, point=None, **kwargs):
        self.grid = grid
        self.point = self.grid._get_array_of_positions_from_arguments(pos=point, **kwargs)
        
    @late
    def position(self):
        return self.cell.position
    
    @late
    def index(self):
        return self.grid.get_index(self.point)
        
    @late
    def isvalid(self):
        return numpy.logical_and(
            numpy.all(self.index >= self.grid.get_minimum_index()[:len(self.index)]),
            numpy.all(self.index <= self.grid.get_maximum_index()[:len(self.index)])
        )
    
    @late
    def cell(self):
        return self.grid[tuple(self.index)]
    
    def get_value_of_attribute(self, name_of_the_attribute):
        return getattr(self.cell, name_of_the_attribute)
    
    def __getattr__(self, name_of_the_attribute):
        return self.get_value_of_attribute(name_of_the_attribute)

class SamplePointWithInterpolation(object):
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
    
    def __init__(self, grid, point=None, **kwargs):
        self.grid = grid
        self.point = self.grid._get_array_of_positions_from_arguments(pos=point, **kwargs)
        
    @late
    def position(self):
        return self.point
    
    @late
    def index(self):
        return self.grid.get_index(self.point)
                
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
            numpy.all(self.index_for_000_cell >= self.grid.get_minimum_index()[:len(self.index)]),
            numpy.all(self.index_for_111_cell <= self.grid.get_maximum_index()[:len(self.index)])
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
    
    def __init__(self, grid, points=None, samples_factory = SamplePointWithInterpolation, **kwargs):
        self.grid = grid
        points=self.grid._get_array_of_positions_from_arguments(pos=points,**kwargs)
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
    
    def __init__(self, grids, points, samples_factory = SamplePointWithInterpolation, index_factory = None):
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
        index_of_grid = self.grids_on_index[tuple(index)]
        return self.grids[index_of_grid]


# convenience function to convert input arguments to positions (or vector of "points")
def _get_array_of_positions_from_arguments(axes_names, **kwargs):
    if kwargs.get('pos',None):
        return kwargs['pos']
    if kwargs.get('position',None):
        return kwargs['position']
    
    coordinates=[kwargs[x] for x in axes_names]
    if numpy.ndim(coordinates[0])==0:
      return VectorQuantity.new_from_scalar_quantities(*coordinates)
    return column_stack(coordinates)
