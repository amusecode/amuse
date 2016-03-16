from amuse.units import units
from amuse.support import exceptions
import numpy

from amuse.datamodel.grids import *


class StaggeredGrid(object):
    def __init__(self, elements, nodes = None, get_corners=None):
        self.elements=elements
        self.nodes=nodes
        self._get_corners_func = get_corners
        self.get_corners()

    def new_remapping_channel_to(self, other, remapper):
        return GridRemappingChannel(self, other, remapper)

    def get_attribute_names_defined_in_store(self):
        elem = self.elements.get_attribute_names_defined_in_store()
        node = self.nodes.get_attribute_names_defined_in_store()
        return list(set(elem).union(set(node)))

    def get_defined_settable_attribute_names(self):
        elem = self.elements.get_defined_settable_attribute_names()
        node = self.nodes.get_defined_settable_attribute_names()
        return list(set(elem).union(set(node)))

    def get_corners(self):
        elements = self.elements
        nodes = self.nodes

        #if a replacing get_corners method is defined for this grid call it and return
        if not (self._get_corners_func is None):
            corners = self._get_corners_func()
            if not hasattr(elements, '_cell_corners'):
                object.__setattr__(elements, "_cell_corners", corners)
            self.corners = corners
            return corners            

        #element and node grids should be of the same type
        if type(elements) != type(nodes):
            raise Exception("elements and nodes grids should be of the same type")

        #get set of participating axes names
        dims = elements.number_of_dimensions()
        if dims > 2:
            raise Exception("Staggered grids with more than 2 dimensions are currently not supported")
            #this is mainly because we can only access the coordinates in
            #a dimensions by their attribute names 'lat' and 'lon'
            #if the _axes_names were defined properly for in-code stored grids I could do something like
            # xpos=to_quantity(getattr(nodes,self._axes_names[0])), but currently that's not possible right now

        #structured or unstructured
        if type(elements) is StructuredGrid:

            #structured grid assumes that the nodes are at the cell corners
            #of the elements grid and therefore the nodes grid should be
            #exactly 1 gridpoint larger in each dimension to fully encapsulate
            #the element grid
            if len(elements.shape) != len(nodes.shape):
                raise Exception("elements and nodes should have the same number of the dimensions")

            for i in range(len(elements.shape)):
                if elements.shape[i]+1 != nodes.shape[i]:
                    raise Exception("nodes grid should be exactly 1 grid point larger than element grid in each dimension")

            corners = numpy.zeros([dims] + list(nodes.shape), dtype=numpy.double)

            #use node positions as corner positions
            corners[0] = nodes.lon.value_in(units.rad)
            corners[1] = nodes.lat.value_in(units.rad)

        elif type(elements) == UnstructuredGrid:
            #the following bit of code tries to access 'n0', 'n1' up to 'n9' of the elements grid
            #a cleaner implementation would be to call get_element_nodes() directly, but we can not do that from here
            attributes = elements.all_attributes()
            max_corners = 10
            corner_indices = []
            for i in range(max_corners):
                node = 'n' + str(i)
                if node in attributes:
                    corner_indices.append(getattr(elements, node))

            self.num_corners = num_corners = len(corner_indices)
            object.__setattr__(elements,"_num_corners", num_corners)
            self.corner_indices = corner_indices
            size = elements.size
            self.inverse_mapping = inverse_mapping = [[] for i in range(nodes.size)]

            #only 2 dimensions supported currently
            corners = numpy.zeros((2, size*num_corners), dtype=numpy.double)
            node_lon = nodes.lon.value_in(units.rad)
            node_lat = nodes.lat.value_in(units.rad)
            for i in range(size):
                for d in range(num_corners):
                    n = corner_indices[d][i]
                    inverse_mapping[n].append(i)
                    nlon = node_lon[n]
                    corners[0][i*num_corners+d] = nlon
                    nlat = node_lat[n]
                    corners[1][i*num_corners+d] = nlat
            
        else:
            raise Exception("unknown grid type for elements: should be either StructuredGrid or UnstructuredGrid")



        if not hasattr(elements, '_cell_corners'):
            object.__setattr__(elements, "_cell_corners", corners)
        self.corners = corners
        return corners            

    def map_elements_to_nodes_structured_larger(self, elements, nodes, elem_values):
        node_values = numpy.zeros(self.nodes.shape, dtype=numpy.float)
        node_values[1:,1:] = elem_values[:,:]
        node_values[0,:] = 0.0
        #assume the grid is cyclic east-west
        node_values[1:,0] = elem_values[:,-1]
        return node_values

    def map_elements_to_nodes_structured_same_size(self, elements, nodes, elem_values):
        node_values = numpy.zeros(self.nodes.shape, dtype=numpy.float)
        node_values = elem_values[:]
        return node_values
    
    def map_elements_to_nodes(self, element_values):
        #currently very rough remapping schemes, more sophisticated methods will be added later

        if not hasattr(self, 'corners'):
            self.get_corners()
        elements = self.elements
        nodes = self.nodes

        element_values = element_values.reshape(elements.shape)

        if type(elements) is StructuredGrid:
            if len(elements.shape) != len(nodes.shape):
                raise Exception("elements and nodes should have the same number of the dimensions")

            if numpy.all([s1+1==s2 for s1,s2 in zip(elements.shape,nodes.shape)]):
                return self.map_elements_to_nodes_structured_larger(elements, nodes, element_values)
            if numpy.all([s1==s2 for s1,s2 in zip(elements.shape,nodes.shape)]):
                return self.map_elements_to_nodes_structured_same_size(elements, nodes, element_values)
            else:
                raise Exception("nodes grid should have either exactly same shape or 1 grid point larger than element grid in each dimension")

        elif type(elements) == UnstructuredGrid:
            if (len(element_values) != self.elements.size):
                raise Exception("number of values passed does not match size of elements grid")
            #do a simple average value of the elements around the node
            node_values = numpy.zeros(self.nodes.size, dtype=numpy.float)
            for i in range(len(node_values)):
                num_neighbors = len(self.inverse_mapping[i])
                value = 0.0
                #add value of neighboring element
                for neighbor in self.inverse_mapping[i]:
                    value += element_values[neighbor]
                #store result
                node_values[i] = value / (1.0*num_neighbors)
        else:
            raise Exception("unknown grid type for elements: should be either StructuredGrid or UnstructuredGrid")
        return node_values


    #this method is for structured grids where the nodes grid is exactly 1 grid point larger in each dimension
    def map_nodes_to_elements_structured_larger(self, elements, nodes, node_values):
        #do simple translation/shift of the values from the north-east corners of each grid cell to the cell centers
        elem_values = numpy.zeros(self.elements.shape, dtype=numpy.float)
        elem_values = node_values[1:,1:]
        return elem_values

    #this method is for structured grids where the nodes grid is of the same size as the elements grid, if so 
    #the grid is assumed to be cyclic east-west
    def map_nodes_to_elements_structured_same_size(self, elements, nodes, node_values):
        #do simple translation/shift of the values from the north-east corners of each grid cell to the cell centers
        elem_values = numpy.zeros(self.elements.shape, dtype=numpy.float)
        elem_values = node_values.flatten()
        return elem_values


    def map_nodes_to_elements(self, node_values):

        if not hasattr(self, 'corners'):
            self.get_corners()
        elements = self.elements
        nodes = self.nodes

        node_values = node_values.reshape(nodes.shape)

        if type(elements) is StructuredGrid:
            if len(elements.shape) != len(nodes.shape):
                raise Exception("elements and nodes should have the same number of the dimensions")
            if numpy.all([s1+1==s2 for s1,s2 in zip(elements.shape,nodes.shape)]):
                return self.map_nodes_to_elements_structured_larger(elements, nodes, node_values)
            if numpy.all([s1==s2 for s1,s2 in zip(elements.shape,nodes.shape)]):
                return self.map_nodes_to_elements_structured_same_size(elements, nodes, node_values)
            else:
                raise Exception("nodes grid should have either exactly same shape or 1 grid point larger than element grid in each dimension")


        elif type(elements) == UnstructuredGrid:
            if (len(node_values) != self.nodes.size):
                raise Exception("number of values passed does not match size of nodes grid")

            elem_values = numpy.zeros(self.elements.size, dtype=numpy.float)

            #do a simple average value of the nodes around the element
            for i in range(len(elem_values)):
                value = 0.0
                for c in range(self.num_corners):
                    index = self.corner_indices[c][i]
                    value += node_values[index]            
                elem_values[i] = value / (1.0*self.num_corners)

        else:
            raise Exception("unknown grid type for elements: should be either StructuredGrid or UnstructuredGrid")
        return elem_values




