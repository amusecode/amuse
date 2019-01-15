"""This module defines the classes needed to map
functions defined by the codes into particle sets
and grids.

The attribute values of Particles or Gridpoints are 
stored in Particle Sets or Grids. These sets or grids
manage:

1. The storage allocation (deletion and removal of particles)
2. The attribute access (getting or setting the value(s) of attribute(s))
3. Queries or selections of particles (selection of subsets of particles)

All 3 functions can be provided by a code. The classes in this 
module provide a mapping from the functions in a code to
the datamodel used in AMUSE. 

When a code manages a particular set all the data of that set is 
stored in the memory space of that code. The code needs to provide
functions to acces the data in the set.

.. note::

    Most codes already implement a particle set or a grid. 
    The only extra requirement for AMUSE is to provide functions to access
    this set. When a code does not have any knowlegde of sets or grids, the
    management will take place in AMUSE and only some data transfer code
    is needed

All incode storage is build on mapping attributes to functions. These
mappings are provided by a number of helper classes:

**setter/getter**

    :py:class:`ParticleGetAttributesMethod`
        Given particle indices or gridpoints (i,j,k) return a vector quantity
        for each attribute
        
    :py:class:`ParticleSetAttributesMethod`
        Send values to the code given particle indices or gridpoints (i,j,k)
        and a vector quantities for each attribute.
    
**new/delete**

    :py:class:`NewParticleMethod`
        Given vector quantities for attributes return the indices
        of newly allocated particles

**function**

    :py:class:`ParticleMethod`
        Given particle indices or gridpoints (i,j,k) and optional arguments 
        return one or more vector quantities

**selection**
        
    :py:class:`ParticleSpecificSelectMethod`
        Given a particle return a subset of particles. For links between
        particles (nearest neighbor, subparticle)

    :py:class:`ParticleQueryMethod`
        Retrieve indices from the code and return a subset of particles.
        For selection of a limited number of particles by the code (get
        the escaper)

    :py:class:`ParticleSetSelectSubsetMethod`
        Like ParticleQueryMethod but can handle larger subsets of
        particles, the code can provide a special function
        the return the number of particles in the set.
        

The InCode storage system is based on a number of classes:

:py:class:`AbstractInCodeAttributeStorage`
    Handle attribute set/get functionality but no particle or
    grid management
    
:py:class:`InCodeAttributeStorage`
    Subclass of AbstractInCodeAttributeStorage, manages particles
    
:py:class:`InCodeGridAttributeStorage`
    Subclass of AbstractInCodeAttributeStorage, manages grids

    

"""

from amuse.support.methods import AbstractCodeMethodWrapper
from amuse.units import nbody_system
from amuse.units import units
from amuse.units import quantities
from amuse.units.quantities import is_quantity
from amuse.support.core import late
from amuse.support import exceptions

import numpy
import inspect

from amuse.datamodel import parameters
from amuse.datamodel import base
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.datamodel import ParticleInformationChannel
from amuse.datamodel import Particle
from amuse.datamodel import Grid
from amuse.datamodel import AttributeStorage

from amuse.rfi.async_request import ASyncRequestSequence

class ParticleMappingMethod(AbstractCodeMethodWrapper):
    def __init__(self, method, attribute_names = None):
        AbstractCodeMethodWrapper.__init__(self, method)
        
        if attribute_names is None:
            self._attribute_names = []
        else:
            self._attribute_names = attribute_names
        
    @late
    def name_of_the_indexing_parameter(self):
        return 'index_of_the_particle'
        

class ParticleGetAttributesMethod(ParticleMappingMethod):
    """
    Instances wrap other methods and provide mappings
    from attribute names to results.
    
    Simple attribute getter methods take an array of indices
    and return a tuple with arrays of result values.
    
    .. code-block:: python
    
        x, y, z = instance.get_xyz(indices)
    
    Instances of this class make it possible to access the 
    return values by their attribute names.
    
    For this it employs two strategies:
    
    1. It uses the provided array of names and
       maps each name to the positional output.
    
    2. If no array of names is provided it asks the wrapped 
       method for all the names of the output parameters 
       (this scheme only works for legacy 
       functions or for wrapped legacy functions)
    
    """
    def __init__(self, method, attribute_names = None):
        ParticleMappingMethod.__init__(self, method, attribute_names)
    
    @late
    def attribute_names(self):
        if self._attribute_names:
            return self._attribute_names
        else:
            result = []
            for x in self.method_output_argument_names:
                if x == self.name_of_the_indexing_parameter:
                    continue
                else:
                    result.append(x)
            return result
    
    def check_arguments(self, storage, attributes_to_return, *indices):
        if len(indices[0]) > 1: 
            if self.method_is_legacy and not (self.method.specification.can_handle_array or self.method.specification.must_handle_array):
                raise Exception(
                    "getter method {0} cannot handle arrays".format(self.method)
                )
            elif self.method_is_code:
                if not self.method.legacy_specification is None:
                    if not (self.method.legacy_specification.can_handle_array or self.method.legacy_specification.must_handle_array):
                        raise exceptions.AmuseException(
                            "getter method {0} cannot handle arrays".format(self.method)
                        )
    
    def convert_return_value(self, return_value, storage, attributes_to_return):
        if len(self.attribute_names) == 1:
            return_value = (return_value,)
        
        set_of_attributes_to_return = set(attributes_to_return)
        
        result = {}
        
        if self.index_output_attributes:
            index_output_attributes = self.index_output_attributes
        else:
            index_output_attributes = [False] * len(return_value)
        
        for value, attribute, isindex in zip(return_value, self.attribute_names, index_output_attributes):
            if attribute in set_of_attributes_to_return:
                if isindex:
                    result[attribute] = quantities.new_quantity(storage._get_keys_for_indices_in_the_code(value), units.object_key)
                else:
                    result[attribute] = value
                    
        return result
    
    def get_attribute_values(self, storage, attributes_to_return, *indices):
        
        self.check_arguments(storage, indices, attributes_to_return)
        
        try:
            return_value = self.method(*indices, **storage.extra_keyword_arguments_for_getters_and_setters)
        except:
            print self.method
            raise
        return self.convert_return_value(return_value, storage, attributes_to_return)
    

    def get_attribute_values_async(self, storage, attributes_to_return, *indices):
        
        self.check_arguments(storage, indices, attributes_to_return)
        
        def result_handler(inner):
            return self.convert_return_value(inner(), storage, attributes_to_return)
            
        async_request = self.method.asynchronous(*indices, **storage.extra_keyword_arguments_for_getters_and_setters)
        async_request.add_result_handler(result_handler)
        return async_request
    
    
class ParticleGetGriddedAttributesMethod(ParticleGetAttributesMethod):
    
    def __init__(self, method, get_range_method, attribute_names = None):
        ParticleGetAttributesMethod.__init__(self, method, attribute_names)
        self.get_range_method = get_range_method
        
    def get_attribute_values(self, storage, attributes_to_return, *indices):
        
        self.check_arguments(storage, indices, attributes_to_return)
        
        minmax_per_dimension = self.get_range_method( **storage.extra_keyword_arguments_for_getters_and_setters)
        
        result = [slice(0, len(indices[0]))]
        gridshape=[len(indices[0])]
        for ind in indices[1:]:
            result.append(slice(0, 1))
        for i in range(0, len(minmax_per_dimension), 2):
            minval = minmax_per_dimension[i]
            maxval = minmax_per_dimension[i+1]
            result.append(slice(minval, maxval+1))
            gridshape.append(maxval+1-minval)
        
        grid_indices = numpy.mgrid[tuple(result)]
        u=grid_indices[0].copy()
        for i, ind in enumerate(indices):
            grid_indices[i] = numpy.asarray(ind)[u]
        one_dimensional_arrays_of_indices = [x.reshape(-1) for x in grid_indices]
        try:
            return_value = self.method(*one_dimensional_arrays_of_indices, **storage.extra_keyword_arguments_for_getters_and_setters)
        except:
            print self.method
            raise
            
        mapping_from_name_to_value = self.convert_return_value(return_value, storage, attributes_to_return)
        for key, value in mapping_from_name_to_value.iteritems():
            mapping_from_name_to_value[key] = value.reshape(gridshape)
        return mapping_from_name_to_value
        
class ParticleSetAttributesMethod(ParticleMappingMethod):
    """
    Instances wrap other methods and provide mappings
    from attribute names to input parameters.
    
    Simple attribute setter methods take an array of indices
    and one or more arrays of new values.
    
    .. code-block:: python
    
       instance.set_xyz(indices, x, y, z)
    
    Instances of this class make it possible to access the 
    possitional parameters with attribute names.
    
    .. Note::
        the index argument is assumed to always come first!
    
    For this it employs two strategies:
    
    1. It uses the provided array of names and
       maps each name to the positional output.
    
    2. If no array of names is provided it asks the wrapped 
       method for all the names of the input parameters 
       (this scheme works for legacy 
       functions and sometimes for python native functions (if
       they have named arguments))
    
    """
    def __init__(self, method,  attribute_names = None):
        ParticleMappingMethod.__init__(self, method,  attribute_names)
    
    @late
    def attribute_names(self):
        if self._attribute_names:
            return self._attribute_names
        else:
            result = []
            for x in self.method_input_argument_names:
                if x == self.name_of_the_indexing_parameter:
                    continue
                else:
                    result.append(x)
            return result
            
    @late
    def optional_attribute_names(self):
        if hasattr(self.method, 'optional_method_input_argument_names'):
            return self.method.optional_method_input_argument_names
        else:
            return []
        
    @late
    def names_to_index(self):
        result = {}
        for index, name in enumerate(self.attribute_names):
            result[name] = index
        return result
        
    def set_attribute_values(self, storage, attributes, values, *indices):
        list_arguments = list(indices)
        list_args, keyword_args = self.convert_attributes_and_values_to_list_and_keyword_arguments(attributes, values)
        list_arguments.extend(list_args)
        keyword_args.update(storage.extra_keyword_arguments_for_getters_and_setters)
        self.method(*list_arguments, **keyword_args)
            
        
    def set_attribute_values_async(self, storage, attributes, values, *indices, **extra_keyword_arguments_for_getters_and_setters):
        keyword_args = {}
        keyword_args.update(storage.extra_keyword_arguments_for_getters_and_setters)
        keyword_args.update(extra_keyword_arguments_for_getters_and_setters)
        
        list_arguments = list(indices)
        list_args, keyword_args2 = self.convert_attributes_and_values_to_list_and_keyword_arguments(attributes, values)
        keyword_args.update(keyword_args2)
        list_arguments.extend(list_args)
        async_request = self.method.asynchronous(*list_arguments, **keyword_args)
        return async_request
    
    
    def convert_attributes_and_values_to_list_and_keyword_arguments(self, attributes, values):
        not_set_marker = object()
        list_arguments = [not_set_marker] * (len(self.attribute_names))
        
        names_to_index = self.names_to_index
        for attribute, quantity in zip(attributes, values):
            if attribute in names_to_index:
                index = names_to_index[attribute]
                list_arguments[index] = quantity
        
        default_argument_found = False
        missing_attributes = []
        dict_arguments = {}
        for index, x in enumerate(list_arguments):
            if x is not_set_marker:
                name_of_attribute = self.attribute_names[index]
                if not name_of_attribute in self.optional_attribute_names:
                    missing_attributes.append(name_of_attribute)
                else:
                    default_argument_found = True
            elif default_argument_found:
                name_of_attribute = self.attribute_names[index]
                if not name_of_attribute in self.optional_attribute_names:
                    raise exceptions.AmuseException("Optional before required arguments")
                dict_arguments[name_of_attribute] = x
                list_arguments[index] = not_set_marker
        
        if len(missing_attributes) > 0:
            if len(missing_attributes) == 1:
                missing_attributes_string = "{0!r} attribute".format(missing_attributes[0])
            else:
                missing_attributes_string = "{0!r} and {1!r} attributes".format(", ".join(missing_attributes[:-1]), missing_attributes[-1])
            raise exceptions.MissingAttributesAmuseException(
                missing_attributes,
                "To add particles to this code you need to specify the {0}".format(missing_attributes_string))
        
        list_arguments = [x for x in list_arguments if not x is not_set_marker]
        return list_arguments, dict_arguments
        
class ParticleSetGriddedAttributesMethod(ParticleSetAttributesMethod):
    
    def __init__(self, method, get_range_method, attribute_names = None):
        ParticleSetAttributesMethod.__init__(self, method, attribute_names)
        self.get_range_method = get_range_method
        
    
    def set_attribute_values(self, storage, attributes, values, *indices):
        list_args, keyword_args = self.convert_attributes_and_values_to_list_and_keyword_arguments(attributes, values)

        minmax_per_dimension = self.get_range_method( **storage.extra_keyword_arguments_for_getters_and_setters)
        
        result = [slice(0, len(indices[0]))]
        gridshape=[len(indices[0])]
        for ind in indices[1:]:
            result.append(slice(0, 1))
        for i in range(0, len(minmax_per_dimension), 2):
            minval = minmax_per_dimension[i]
            maxval = minmax_per_dimension[i+1]
            result.append(slice(minval, maxval+1))
            gridshape.append(maxval+1-minval)
        
        grid_indices = numpy.mgrid[tuple(result)]
        u=grid_indices[0].copy()
        for i, ind in enumerate(indices):
            grid_indices[i] = numpy.asarray(ind)[u]
        one_dimensional_arrays_of_indices = [x.reshape(-1) for x in grid_indices]
        
        list_arguments = list(grid_indices)
        list_arguments.extend(list_args)
        one_dimensional_arrays_of_args = [x.reshape(-1) for x in list_arguments]
        
        for key, value in keyword_args.iteritems():
            keyword_args[key] = value.reshape(-1)
        
        keyword_args.update(storage.extra_keyword_arguments_for_getters_and_setters)
        self.method(*one_dimensional_arrays_of_args, **keyword_args)
        
        
        
    def set_attribute_values_async(self, storage, attributes, values, *indices):
        keyword_args = {}
        keyword_args.update(storage.extra_keyword_arguments_for_getters_and_setters)
        list_args, keyword_args2 = self.convert_attributes_and_values_to_list_and_keyword_arguments(attributes, values)
        keyword_args.update(keyword_args2)
#        print keyword_args2
        minmax_per_dimension = self.get_range_method( **storage.extra_keyword_arguments_for_getters_and_setters)
        
        result = [slice(0, len(indices[0]))]
        gridshape=[len(indices[0])]
        for ind in indices[1:]:
            result.append(slice(0, 1))
        for i in range(0, len(minmax_per_dimension), 2):
            minval = minmax_per_dimension[i]
            maxval = minmax_per_dimension[i+1]
            result.append(slice(minval, maxval+1))
            gridshape.append(maxval+1-minval)
        
        grid_indices = numpy.mgrid[tuple(result)]
        u=grid_indices[0].copy()
        for i, ind in enumerate(indices):
            grid_indices[i] = numpy.asarray(ind)[u]
        one_dimensional_arrays_of_indices = [x.reshape(-1) for x in grid_indices]

        list_arguments = list(grid_indices)
        list_arguments.extend(list_args)
        one_dimensional_arrays_of_args = [x.reshape(-1) for x in list_arguments]
        
        for key, value in keyword_args.iteritems():
            keyword_args[key] = value.reshape(-1)
        
        async_request = self.method.asynchronous(*one_dimensional_arrays_of_args, **keyword_args)
        return async_request
        

class NewParticleMethod(ParticleSetAttributesMethod):
    """
    Instances wrap a method to create particles. The method may
    take attributes values to set initial values on
    the created particles. 
    
    The new particle functions work a lot like 
    the set attribute methods, only the new particle 
    function is supposed to return an array
    of the indices of the created particles.
    
    .. code-block:: python
    
       indices = instance.new_particle(x, y, z)
       
    """
    def __init__(self,  method, attribute_names = None):
        ParticleSetAttributesMethod.__init__(self, method, attribute_names)

    def add_entities(self, attributes, values):
        list_arguments,keyword_arguments = self.convert_attributes_and_values_to_list_and_keyword_arguments(attributes, values)
        indices = self.method(*list_arguments, **keyword_arguments)
        return indices
        
class ParticleQueryMethod(object):
    """
    Instances wrap a function that can take one or more arguments
    and returns an index (or a list of indices, if the arguments are
    lists). This method is most useful to select one particle form
    all particles in the set
    
    .. code-block:: python
    
        index = instance.get_escaper()
    
    The index or indices are converted to a particle subset.
    """
    def __init__(self, method, names = (), public_name = None, query_superset=False):
        self.method = method
        self.name_of_the_out_parameters = names
        self.public_name = public_name
        if query_superset:
            self.apply = self.apply_for_superset
        else:
            self.apply = self.apply_normal

    def apply_normal(self, particles, *args, **kwargs):
        indices = self.method(*args, **kwargs)
        keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices)
        return particles._subset(keys)
    
    def apply_for_superset(self, particles, *args, **kwargs):
        indices = self.method(*args, **kwargs)
        subset_results = []
        for subset in particles._private.particle_sets:
            keys = []
            for index in indices:
                if index in subset._private.attribute_storage.mapping_from_index_in_the_code_to_particle_key:
                    keys.append(subset._private.attribute_storage.mapping_from_index_in_the_code_to_particle_key[index])
            subset_results.append(subset._subset(keys))
        return ParticlesSuperset(subset_results)
    

class ParticleSpecificSelectMethod(object):
    """
    Instances wrap a function that can take a particle index
    and returns one or more indices
    (but a limited and fixed number of indices). This method is most 
    useful to return links between particles (subparticles or
    nearest neighbors)
    
    .. code-block:: python
    
        output_index = instance.get_nearest_neigbord(input_index)
    
    The idex or indices are converted to a particle subset.
    """
    def __init__(self, method, names = (), public_name = None):
        self.method = method
        self.name_of_the_out_parameters = names
        self.public_name = public_name

    def apply_on_all(self, particles):
        
        all_indices = particles._private.attribute_storage.mapping_from_index_in_the_code_to_particle_key.keys()
        
        lists_of_indices = self.method(list(all_indices))
        
        lists_of_keys = []
        for indices in lists_of_indices:
            keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices)        
            lists_of_keys.append(keys)
        
        result = []
        for keys in zip(lists_of_keys):
            result.append(particles._subset(keys))
            
        return result
    
    def apply_on_one(self, set,  particle):
        
        index = set._private.attribute_storage.get_indices_of(particle.key)
        
        result = self.method(index)
        
        keys = set._private.attribute_storage._get_keys_for_indices_in_the_code(result)  
        
        result = []
        return particle.as_set()._subset(keys)
        
        
class ParticleMethod(AbstractCodeMethodWrapper):
    """
    Instances wrap a function that returns quanties given particle
    indices and optional arguments. Instances have a lot in common
    with attribute getters, but can take extra arguments.
    
    .. code-block:: python
    
        pressure = instance.get_pressure(index, gamma)
    """
    def __init__(self, method, public_name = None):
        AbstractCodeMethodWrapper.__init__(self, method)
        self.public_name = public_name

    def apply_on_all(self, particles, *list_arguments, **keyword_arguments):
        storage = particles._private.attribute_storage
        all_indices = list(storage.mapping_from_index_in_the_code_to_particle_key.keys())
        return self.method(all_indices, *list_arguments, **keyword_arguments)
    
    def apply_on_one(self, set,  particle, *list_arguments, **keyword_arguments):
        storage = particle.particles_set._private.attribute_storage
        index = storage.get_indices_of([particle.key])
        return self.method(index[0], *list_arguments, **keyword_arguments)
        
class ParticleSetSelectSubsetMethod(object):
    """
    Generic method to query and retrieve particles from the
    set. This selection can have up to tree stages:
    
    1. start the query given a number of optional arguments
    2. get the number of selected particles
    3. get the index of each particle 
    
    The pseudo-code for this selection is:
    
    .. code-block:: python
    
        set_selection_criteria(r = 10.0 | units.m)
        n = get_number_of_selected_particles()
        for i in range(n):
            particle_index = get_index_of_selected_particle(i)
    
    The first and second step are optional. If no number of 
    particles method is provided the class assumes the selection
    only returns 1 particle.
    
    Generalisation of ParticleQueryMethod
    """
    
    def __init__(self,  method, set_query_arguments_method = None, get_number_of_particles_in_set_method = None, public_name = None):
        self.method = method
        self.set_query_arguments_method = set_query_arguments_method
        self.get_number_of_particles_in_set_method = get_number_of_particles_in_set_method
        self.public_name = public_name

    def apply_on_all(self, particles, *list_arguments, **keyword_arguments):
        query_identifiers = None
        if not self.set_query_arguments_method is None:
            query_identifiers = self.set_query_arguments_method(*list_arguments, **keyword_arguments)
        
        if query_identifiers is None:
            query_identifiers = ()
        elif not hasattr(query_identifiers, '__iter__'):
            query_identifiers = (query_identifiers,)
            
        if not self.get_number_of_particles_in_set_method is None:
            number_of_particles_in_set = self.get_number_of_particles_in_set_method(*query_identifiers)
            indices = self.method(range(number_of_particles_in_set))
        else:
            index = self.method(*query_identifiers)
            indices = [index]
        
        query_identifiers = [ [x]*len(indices) for x in query_identifiers ]
        keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices, *query_identifiers)    
        
        return particles._subset(keys)


class ParticlesAddedUpdateMethod(object):
   
    
    def __init__(self,  get_number_of_particles_added_method = None, get_id_of_added_particles_method = None):
        self.get_number_of_particles_added_method = get_number_of_particles_added_method
        self.get_id_of_added_particles_method = get_id_of_added_particles_method
    

    def apply_on_all(self, particles, *list_arguments, **keyword_arguments):
        query_identifiers = None
        if not self.set_query_arguments_method is None:
            query_identifiers = self.set_query_arguments_method(*list_arguments, **keyword_arguments)
        
        if query_identifiers is None:
            query_identifiers = ()
        elif not hasattr(query_identifiers, '__iter__'):
            query_identifiers = (query_identifiers,)
            
        if not self.get_number_of_particles_in_set_method is None:
            number_of_particles_in_set = self.get_number_of_particles_in_set_method(*query_identifiers)
            indices = self.method(range(number_of_particles_in_set))
        else:
            index = self.method(*query_identifiers)
            indices = [index]
        
        query_identifiers = [ [x]*len(indices) for x in query_identifiers ]
        keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices, *query_identifiers)    
        
        return particles._subset(keys)

class ParticleGetIndexMethod(object):
    """
    Instances return the index of a particle in the code
    """
    ATTRIBUTE_NAME = "index_in_code"
    
    def __init__(self):
        pass
    
    @late
    def attribute_names(self):
        return [self.ATTRIBUTE_NAME]
    
    def get_attribute_values(self, storage, attributes_to_return, *indices):
        
        return {self.ATTRIBUTE_NAME : indices[0]}

class AbstractInCodeAttributeStorage(base.AttributeStorage):
    """
    Abstract base storage for incode attribute storage.
    It provides functions to handle getters and setters of 
    attributes but not for creating or deleting of particles as
    this differs between grids and particle sets.
    
    """
    
    def __init__(self, 
            code_interface, 
            setters,
            getters,
            extra_keyword_arguments_for_getters_and_setters = {},
    ):
        
        self.code_interface = code_interface
        
        self.getters = list(getters)
        self.setters = setters
        
        self.attributes = set([])
        for x in self.getters:
            self.attributes |= set(x.attribute_names)
        for x in self.setters:
            self.attributes |= set(x.attribute_names)
            
        self.writable_attributes = set([])
        for x in self.setters:
            self.writable_attributes |= set(x.attribute_names)
        
        
        self.extra_keyword_arguments_for_getters_and_setters = extra_keyword_arguments_for_getters_and_setters
        
    
    def select_getters_for(self, attributes):
        set_of_attributes = set(attributes)
        
        # first check for an exact match
        result = [getter for getter in self.getters if set(getter.attribute_names) == set_of_attributes]
        if result:
            return result
        
        # sort methods on attribute lengths, longest first
        sorted_getters = sorted(self.getters, key=lambda x : len(x.attribute_names), reverse = True)
        
        # next, select the longest fitting method(s), to minize the number of calls
        for access_method in sorted_getters:
            if set_of_attributes >= set(access_method.attribute_names):
                result.append(access_method)
                set_of_attributes -= set(access_method.attribute_names)
        
        # next, select the sortest method(s), to minimize the extra parameters
        if set_of_attributes:
            for access_method in reversed(sorted_getters):
                if set_of_attributes & set(access_method.attribute_names):
                    result.append(access_method)
                    set_of_attributes -= set(access_method.attribute_names)
                    
        if set_of_attributes:
            raise exceptions.AmuseException("Do not have attributes {0}".format(sorted(set_of_attributes)))
        
        return result
    
    def select_setters_for(self, attributes):
        set_of_attributes = set(attributes)
        result = []
        for access_method in self.setters:
            if set_of_attributes >= set(access_method.attribute_names):
                result.append(access_method)
                set_of_attributes -= set(access_method.attribute_names)
                
        if set_of_attributes:
            raise exceptions.AmuseException("Cannot set attributes {0}".format(sorted(set_of_attributes)))
            
        return result
    
    
    def get_defined_attribute_names(self):
        return sorted(self.attributes)
        
    def get_defined_settable_attribute_names(self):
        return sorted(self.writable_attributes)

    
class InCodeAttributeStorage(AbstractInCodeAttributeStorage):
    """
    Manages sets of particles stored in codes.
    
    Maps indices returned by the code to keys defined in AMUSE.
    """
    def __init__(self, 
            code_interface, 
            new_particle_method, 
            delete_particle_method, 
            number_of_particles_method, 
            setters,
            getters,
            name_of_the_index):
        
        
        for x in getters:
            x.name_of_the_indexing_parameter = name_of_the_index
            
        for x in setters:
            x.name_of_the_indexing_parameter = name_of_the_index
        
        getters = list(getters)
        
        AbstractInCodeAttributeStorage.__init__(self, code_interface, setters, getters)
    
        self.mapping_from_particle_key_to_index_in_the_code = {}
        self.mapping_from_index_in_the_code_to_particle_key = {}
        self.particle_keys = numpy.zeros(0)
        self.code_indices = numpy.zeros(0)
        
        self._get_number_of_particles = number_of_particles_method
        self.delete_particle_method = delete_particle_method
        self.new_particle_method = new_particle_method
        
        self.getters.append(ParticleGetIndexMethod())

    def __len__(self):
        return len(self.mapping_from_particle_key_to_index_in_the_code)

    def can_extend_attributes(self):
        return False
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        
        indices = self.new_particle_method.add_entities(attributes, values)
        
        if len(self.particle_keys) > 0:
            previous_length = len(self.particle_keys)
            self.particle_keys = numpy.concatenate((self.particle_keys, numpy.array(list(keys))))
            self.code_indices =  numpy.concatenate((self.code_indices, numpy.array(indices)))
            result = self.code_indices[previous_length:]
        else:
            self.particle_keys = numpy.array(keys)
            self.code_indices = numpy.array(indices)
            result = self.code_indices
            
        index = 0
        for key in keys:
            if key in self.mapping_from_particle_key_to_index_in_the_code:
                raise Exception("particle with same key added twice: {0}".format(key))
            self.mapping_from_particle_key_to_index_in_the_code[key] = indices[index]
            self.mapping_from_index_in_the_code_to_particle_key[indices[index]] = key
            index = index + 1
        
        return result

    def get_indices_of(self, keys):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
        
        notfoundkeys = []
        foundkeys = []
        for particle_key in keys:
            try:
                indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
                foundkeys.append(particle_key)
            except KeyError:
                notfoundkeys.append(particle_key)
          
        if not len(notfoundkeys) == 0:
            raise exceptions.KeysNotInStorageException(
                numpy.asarray(foundkeys), 
                numpy.asarray(indices_in_the_code), 
                numpy.asarray(notfoundkeys)
            )
        
        return numpy.asarray(indices_in_the_code)
        
   

    def get_key_indices_of(self, keys):
        result = []
        if keys is None:
            keys = self.particle_keys
        
        keys_set = set(keys)
        for index in range(len(self.particle_keys)):
            key = self.particle_keys[index]
            if key in keys_set:
                result.append(index)
          
        return result
         
        
    def get_positions_of_indices(self, indices):
        result = []
        if indices is None:
            indices = self.code_indices
        
        indices_set = set(indices)
        for index in range(len(self.code_indices)):
            index_in_code = self.code_indices[index]
            if index_in_code in indices_set:
                result.append(index)
          
        return result
        
    def get_value_of(self, index, attribute):
        return self.get_value_in_store(index, attribute)
   
    def get_values_in_store(self, indices_in_the_code, attributes):

        if indices_in_the_code is None:
            indices_in_the_code = self.code_indices
            
        if len(indices_in_the_code) == 0:
            return [[] for attribute in attributes]
             
        mapping_from_attribute_to_result = {}
        
        for getter in self.select_getters_for(attributes):
            result = getter.get_attribute_values(self, attributes, indices_in_the_code)
            mapping_from_attribute_to_result.update(result)
            
        results = []
        for attribute in attributes:
            results.append(mapping_from_attribute_to_result[attribute])
        return results
        
    

    def get_values_in_store_async(self, indices_in_the_code, attributes):
    
        if indices_in_the_code is None:
            indices_in_the_code = self.code_indices
            
        if len(indices_in_the_code) == 0:
            return [[] for attribute in attributes]
             
        mapping_from_attribute_to_result = {}
        
        getters = self.select_getters_for(attributes)
        if len(getters) > 1:
            def new_request(index, getters, attributes, indices_in_the_code):
                if index >= len(getters):
                    return None
                getter = getters[index]
                request = getter.get_attribute_values_async(self, attributes, indices_in_the_code)
                def result_handler(inner, mapping):
                    mapping.update(inner())
                request.add_result_handler(result_handler, (mapping_from_attribute_to_result,))
                return request

            request = ASyncRequestSequence(new_request, args=(getters,attributes, indices_in_the_code))
        else:
            for getter in getters:
                request = getter.get_attribute_values_async(self, attributes, indices_in_the_code)
                def result_handler(inner, mapping):
                    mapping.update(inner())
                request.add_result_handler(result_handler, (mapping_from_attribute_to_result,))
    
        def all_handler(inner, mapping):
            inner()
            results = []
            for attribute in attributes:
                results.append(mapping[attribute])
            return results
            
        request.add_result_handler(all_handler, (mapping_from_attribute_to_result,))
        return request
        
    def set_values_in_store(self, indices_in_the_code, attributes, values):
        if indices_in_the_code is None:
            indices_in_the_code = self.code_indices

        if len(indices_in_the_code) == 0:
            return
            
        for setter in self.select_setters_for(attributes):
            setter.set_attribute_values(self, attributes, values, indices_in_the_code)

    def set_values_in_store_async(self, indices_in_the_code, attributes, values):
        if indices_in_the_code is None:
            indices_in_the_code = self.code_indices

        if len(indices_in_the_code) == 0:
            return
        setters = self.select_setters_for(attributes)
        if len(setters) > 1:
            def new_request(index, setters, attributes, values, indices_in_the_code):
                if index >= len(setters):
                    return None
                setter = setters[index]
                request = setter.set_attribute_values_async(self, attributes, values, indices_in_the_code)
                return request

            request = ASyncRequestSequence(new_request, args=(setters,attributes, values, indices_in_the_code))
        else:
            for setter in setters:
                request = setter.set_attribute_values(self, attributes, values, indices_in_the_code)
        return request
    

    def remove_particles_from_store(self, indices_in_the_code):
        if indices_in_the_code is None:
            return
        self.delete_particle_method(indices_in_the_code)
        
        mapping_key = self.mapping_from_particle_key_to_index_in_the_code
        mapping_index = self.mapping_from_index_in_the_code_to_particle_key
        for i in indices_in_the_code:
            key = mapping_index[i]
            del mapping_index[i]
            del mapping_key[key]
        
        indices_to_delete = self.get_positions_of_indices(indices_in_the_code)
        
        self.particle_keys =  numpy.delete(self.particle_keys, indices_to_delete)
        self.code_indices =  numpy.delete(self.code_indices, indices_to_delete)
            
        
    def get_all_keys_in_store(self):
        return self.particle_keys

    def get_all_indices_in_store(self):
        return self.code_indices
        
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_key_to_index_in_the_code
        
    def _get_keys_for_indices_in_the_code(self, indices):
        result = []
        for i in indices:
            result.append(self.mapping_from_index_in_the_code_to_particle_key.get(i, 0))
        return result
        
    def _remove_indices(self, indices):
        keys = []
        for i in indices:
            if i in self.mapping_from_index_in_the_code_to_particle_key:
                key = self.mapping_from_index_in_the_code_to_particle_key[i]
                del self.mapping_from_index_in_the_code_to_particle_key[i]
                del self.mapping_from_particle_key_to_index_in_the_code[key]
                keys.append(key)
                
        indices_to_delete = self.get_key_indices_of(keys)
        self.particle_keys =  numpy.delete(self.particle_keys, indices_to_delete)
        self.code_indices =  numpy.delete(self.code_indices, indices_to_delete)
        
    
    def _add_indices(self, indices):
        keys = []
        for i in indices:
            if i in self.mapping_from_index_in_the_code_to_particle_key:
                raise exceptions.AmuseException("adding an index '{0}' that is already managed, bookkeeping is broken".format(i))
            newkey = base.UniqueKeyGenerator.next()
            self.mapping_from_index_in_the_code_to_particle_key[i] = newkey
            self.mapping_from_particle_key_to_index_in_the_code[newkey] = i
            
            keys.append(newkey)
        if len(self.particle_keys) > 0:
            self.particle_keys = numpy.concatenate((self.particle_keys, 
                numpy.asarray(list(keys), dtype=self.particle_keys.dtype)))
            self.code_indices =  numpy.concatenate((self.code_indices, 
                numpy.asarray(list(indices), dtype=self.code_indices.dtype)))
        else:
            self.particle_keys = numpy.array(keys)
            self.code_indices = numpy.array(indices)
                

class InCodeGridAttributeStorage(AbstractInCodeAttributeStorage):
    """
    Manages grids stored in codes. 
    """
    def __init__(self, 
            code_interface, 
            get_range_method,
            setters,
            getters,
            extra_keyword_arguments_for_getters_and_setters = {},
    ):
        AbstractInCodeAttributeStorage.__init__(self, code_interface, setters, getters, extra_keyword_arguments_for_getters_and_setters)
        self.get_range_method = get_range_method
        self._indices_grid = None
            
    def can_extend_attributes(self):
        return False
            
    def storage_shape(self):
        try:
            minmax_per_dimension = self.get_range_method(**self.extra_keyword_arguments_for_getters_and_setters)
            result = []
            for i in range(0, len(minmax_per_dimension), 2):
                minval = minmax_per_dimension[i]
                maxval = minmax_per_dimension[i+1]
                result.append(maxval - minval + 1)
            return tuple(result)
        except:
            import traceback
            traceback.print_exc()
            raise
        
    def add_particles_to_store(self, keys, attributes = [], quantities = []):
        raise exceptions.AmuseException("adding points to the grid is not implemented")
            
    def remove_particles_from_store(self, keys):
        raise exceptions.AmuseException("removing points from the grid is not implemented")
    
    def _to_arrays_of_indices(self, index):
        #imin, imax, jmin, jmax, kmin, kmax = self.get_range_method(**self.extra_keyword_arguments_for_getters_and_setters)
        if self._indices_grid is None:
            minmax_per_dimension = self.get_range_method(**self.extra_keyword_arguments_for_getters_and_setters)
            result = []
            for i in range(0, len(minmax_per_dimension), 2):
                minval = minmax_per_dimension[i]
                maxval = minmax_per_dimension[i+1]
                result.append(slice(minval, maxval+1))

            self._indices_grid  = numpy.mgrid[tuple(result)]
        indices = self._indices_grid 
        
        if index is None:
            return indices
        else:
            return [x[index] for x in indices]
        
    def get_values_in_store(self, indices, attributes):
        array_of_indices = self._to_arrays_of_indices(indices)
        mapping_from_attribute_to_result = {}    
        one_dimensional_array_of_indices = [x.reshape(-1) for x in array_of_indices]
        for getter in self.select_getters_for(attributes):
            result = getter.get_attribute_values(self, attributes, *one_dimensional_array_of_indices)
            mapping_from_attribute_to_result.update(result)
            
        results = []
        for attribute in attributes:
            returned_value = mapping_from_attribute_to_result[attribute]
            
            if len(array_of_indices)==0:
                value=returned_value
            elif len(array_of_indices[0].shape) == 0:
                value = returned_value[0]
            else:
                if len(returned_value)!=numpy.product(array_of_indices[0].shape):
                    raise Exception("unexpected mismatch of array shapes")
                if isinstance(returned_value,list):
                  returned_value=numpy.asarray(returned_value)
                value = returned_value.reshape(array_of_indices[0].shape+returned_value.shape[1:])
                
            results.append(value)
            
        return results
        
    def set_values_in_store(self,  indices, attributes, quantities):
        array_of_indices = self._to_arrays_of_indices(indices)
        one_dimensional_array_of_indices = [x.reshape(-1) for x in array_of_indices]
        if len(one_dimensional_array_of_indices)==0:
            one_dimensional_values = [x for x in quantities]
        else:
            one_dimensional_values = [(x.reshape(-1) if is_quantity(x) else numpy.asanyarray(x).reshape(-1)) for x in quantities]

        
        for setter in self.select_setters_for(attributes):
            setter.set_attribute_values(self, attributes, one_dimensional_values, *one_dimensional_array_of_indices)
     
        
    def set_values_in_store_async(self,  indices, attributes, quantities):
        array_of_indices = self._to_arrays_of_indices(indices)
    
        one_dimensional_array_of_indices = [x.reshape(-1) for x in array_of_indices]
        if len(one_dimensional_array_of_indices)==0:
            one_dimensional_values = [x for x in quantities]
        else:
            one_dimensional_values = [(x.reshape(-1) if is_quantity(x) else numpy.asanyarray(x).reshape(-1)) for x in quantities]
        selected_setters = list([setter for setter in self.select_setters_for(attributes)])
        
        def next_request(index, setters):
            if index < len(setters):
                setter = setters[index]
                return setter.set_attribute_values_async(self, attributes, one_dimensional_values, *one_dimensional_array_of_indices)
            else:
                return None
        
        request = ASyncRequestSequence(next_request, args = (selected_setters,))
        return request
        
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index
        
    def get_all_keys_in_store(self):
        return Ellipsis
        
    def __len__(self):
        shape = self.storage_shape()
        return shape[0] * shape[1] * shape[2]
        
    def copy(self):
        from .memory_storage import InMemoryGridAttributeStorage
        copy = InMemoryGridAttributeStorage()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy() 
        return copy
        
    
    def get_defined_attribute_names(self):
        return sorted(self.attributes)
        
    def _get_writeable_attribute_names(self):
        return self.writable_attributes
        
    def get_defined_settable_attribute_names(self):
        return sorted(self.writable_attributes)





class ParticleSpecificSelectSubsetMethod(object):
    """
    Instances wrap a function that can take a particle index, plus a list
    offset and returns one index. This method is most 
    useful to return links between particles (subparticles or
    nearest neighbors). Instances also need a function to get
    the number of links.
    
    .. code-block:: python
    
        output_index = instance.get_nearest_neigbors(index_of_the_particle, input_index)
    
    The index or indices are converted to a particle subset.
    """
    def __init__(self, method,  get_number_of_particles_in_set_method = None, public_name = None):
        self.method = method
        self.public_name = public_name
        self.get_number_of_particles_in_set_method = get_number_of_particles_in_set_method

    def apply_on_all(self, particles):
        raise Exception("Getting all links to other particles from all particles in a set is not implemented yet")
    
    def apply_on_one(self, set,  particle):
        
        from_indices = set._private.attribute_storage.get_indices_of([particle.key,])
        
        if not self.get_number_of_particles_in_set_method is None:
            number_of_particles_in_set = self.get_number_of_particles_in_set_method(from_indices)[0]
            indices = self.method([from_indices[0]] * number_of_particles_in_set, range(number_of_particles_in_set))
        else:
            index = self.method()
            indices = [index]
            
        keys = set._private.attribute_storage._get_keys_for_indices_in_the_code(indices)                          
     
        return particle.as_set()._subset(keys)

