from amuse.support.data import core
from amuse.support.core import late
from amuse.support.units import units
from amuse.support.io import base

import numpy

class VtkStructuredGrid(base.FileFormatProcessor):
    """
    Process a text file containing a table of values separated by a predefined character
    """
    
    provided_formats = ['txt']
    
    def __init__(self, filename = None, stream = None, set = None, format = None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        
        self.filename = filename
        self.stream = stream
        self.set = set
    
    def store(self):
        
        if self.stream is None:
            self.stream = open(self.filename, "w")
            close_function = self.stream.close 
        else:
            close_function = lambda : None
            
        try:
            return self.store_on_stream()
        finally:
            close_function()
            
    def store_on_stream(self):
        self.write_header()
        self.write_grid()
        self.write_footer()
        
    @late
    def quantities(self):
        if self.set is None:
            return []
        else:
            return map(lambda x:getattr(self.set, x),self.attribute_names)

    @base.format_option
    def attribute_names(self):
        "list of the names of the attributes to load or store"
        if self.set is None:
            return map(lambda x : "col({0})".format(x), range(len(self.quantities)))
        else:
            all_attributes = self.set.stored_attributes()
            return [ x for x in all_attributes if x not in set(['x','y','z'])]
        
    @base.format_option
    def attribute_types(self):
        "list of the types of the attributes to store"
        quantities = self.quantities
        if self.quantities:
            return map(lambda x : x.unit.to_simple_form(), quantities)
        elif self.set is None:
            return map(lambda x : units.none, self.attribute_names)
    

    @base.format_option
    def float_format_string(self):
        "format specification string to convert numbers to strings, see format_spec in python documentation"
        return ".{0}e".format(self.precision_of_number_output)

    @base.format_option
    def precision_of_number_output(self):
        "The precision is a decimal number indicating how many digits should be displayed after the decimal point"
        return 12
        
        
    @base.format_option
    def extent(self):
        "The number of points of the grid in the x, y and z direction, array with 6 float xmin, xmax, ymin, ymax, zmin, zmax"
        quantities = self.quantities
        if self.quantities:
            first_quantity = self.quantities[0]
            nx,ny,nz = first_quantity.shape[0:3]
            return (0,nx-1,0,ny-1,0,nz-1)
        elif self.set is None:
            nx,ny,nz = self.set.shape
            return (0,nx-1,0,ny-1,0,nz-1)
    
            
    @base.format_option
    def minmax(self):
        "The extent of the grid, array with 6 float xmin, xmax, ymin, ymax, zmin, zmax"
        if not self.set is None:
            xmin = self.set.x.min()
            xmax = self.set.x.max()
            ymin = self.set.y.min()
            ymax = self.set.y.max()
            zmin = self.set.z.min()
            zmax = self.set.z.max()
            result = values.AdaptingVectorQuantity
            result.append(xmin)
            result.append(xmax)
            result.append(ymin)
            result.append(ymax)
            result.append(zmin)
            result.append(zmax)
            return result.append(zmax)
        
        if 'x' in self.attribute_names and 'y' in  self.attribute_names and 'z' in self.attribute_names:
            pass
        
    @base.format_option
    def points(self):
        "The position vector of all grid cells"
        if not self.set is None:
            return self.set.position
        else:
            pass
            
    @base.format_option
    def length_unit(self):
        if not self.set is None:
            return self.set.position.unit
        else:
            pass
            
    def write_header(self):
        self.stream.write('<?xml version="1.0" encoding="utf-8"?>\n')
        self.stream.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        
    def write_grid(self):
        self.stream.write('<StructuredGrid WholeExtent="')
        self.stream.write(" ".join(map(str,self.extent)))
        self.stream.write('">')

        self.stream.write('<Piece Extent="')
        self.stream.write(" ".join(map(str,self.extent)))
        self.stream.write('">')
        self.stream.write('<PointData>')
        for name, quantity, unit in zip(self.attribute_names, self.quantities, self.attribute_types):
            self.write_float64_data(quantity.value_in(unit), name = name)
        self.stream.write('</PointData>')
        self.stream.write('<CellData></CellData>')
        self.stream.write('<Points>')
        self.write_float64_data(self.points.value_in(self.length_unit))
        self.stream.write('</Points>')
        self.stream.write('</Piece>')

        self.stream.write('</StructuredGrid>')
    
    def write_footer(self):
        self.stream.write('</VTKFile>')
    
    def write_float64_data(self, array, name = None):
        number_of_components = numpy.prod(array.shape[3:])
        if len(array.shape[3:]) == 0:
            number_of_components = 1
            array = numpy.transpose(array, (2,1,0,)).flatten()
        else:
            x,y,z = numpy.split(array, 3, axis = array.ndim - 1)
            array = numpy.transpose(array, (2,1,0,) + tuple(array.shape[3:])).flatten()
        self.stream.write('<DataArray type="Float64" NumberOfComponents="{0}"'.format(number_of_components))
        if not name is None:
            self.stream.write(' Name="{0}"'.format(name))
        
        self.stream.write('>')
        self.stream.write(' '.join(map(self.convert_number_to_string, iter(array.flatten()))));
        self.stream.write('</DataArray>\n')
        
    def convert_number_to_string(self, number):
        return str(number)
    
    
