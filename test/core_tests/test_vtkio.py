from amuse.support import io
from amuse.support.io import vtk
from amuse.support.units import units
from amuse.support.units import generic_unit_system
from amuse.support.data import core
from amuse.test import amusetest
import StringIO
import textwrap
import os
import numpy

from amuse.support.units  import generic_unit_system
        
class VtkStructuredGridTests(amusetest.TestCase):
    
    def test1(self):
        grid = core. Grid.create([2,3,4], [1,1,1] | generic_unit_system.length)
        grid.rho = grid.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        data_file = StringIO.StringIO()
        instance = vtk.VtkStructuredGrid("test.vts", data_file, grid)
        instance.store()
        
        contents = data_file.getvalue()
        self.assertTrue(contents.find('WholeExtent="0 2 0 3 0 4"')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="3">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="rho">')> 0)
    

    def test2(self):
        grid = core. Grid.create([4,5,6], [1,1,1] | generic_unit_system.length)
        grid.mass = generic_unit_system.density(numpy.random.rand(4,5,6))
        data_file = StringIO.StringIO()
        instance = vtk.VtkStructuredGrid("test.vts", data_file, grid)
        instance.store()
        
        contents = data_file.getvalue()
        print contents
        self.assertTrue(contents.find('<Piece Extent="0 4 0 5 0 6">')> 0)
        self.assertTrue(contents.find('<CellData>')> 0)
        self.assertTrue(contents.find('<PointData>')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="3">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="mass">')> 0)
    
