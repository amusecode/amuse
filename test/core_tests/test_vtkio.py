from amuse.test import amusetest
from StringIO import StringIO
import textwrap
import os
import numpy


from amuse import io
from amuse.io import vtk
from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import generic_unit_system
from amuse import datamodel

class VtkStructuredGridTests(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel. Grid.create([2,3,4], [1,1,1] | generic_unit_system.length)
        grid.rho = grid.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        data_file = StringIO()
        instance = vtk.VtkStructuredGrid("test.vts", data_file, grid)
        instance.store()
        
        contents = data_file.getvalue()
        self.assertTrue(contents.find('WholeExtent="0 2 0 3 0 4"')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="3">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="rho">')> 0)
    

    def test2(self):
        grid = datamodel. Grid.create([4,5,6], [1,1,1] | generic_unit_system.length)
        grid.mass = generic_unit_system.density(numpy.random.rand(4,5,6))
        data_file = StringIO()
        instance = vtk.VtkStructuredGrid("test.vts", data_file, grid)
        instance.store()
        
        contents = data_file.getvalue()
        self.assertTrue(contents.find('<Piece Extent="0 4 0 5 0 6">')> 0)
        self.assertTrue(contents.find('<CellData>')> 0)
        self.assertTrue(contents.find('<PointData>')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="3">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="mass">')> 0)
    
class VtkUnstructuredGridTests(amusetest.TestCase):
    
    def test1(self):
        grid = datamodel. Grid.create([2,3,4], [1,1,1] | generic_unit_system.length)
        grid.rho = grid.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        data_file = StringIO()
        instance = vtk.VtkUnstructuredGrid("test.vtu", data_file, grid)
        instance.store()
        contents = data_file.getvalue()
        
        self.assertTrue(contents.find('<Piece NumberOfPoints="60" NumberOfCells="24"')> 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="connectivity">')> 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="offsets">')> 0)
        self.assertTrue(contents.find('11 11 11')> 0)
        self.assertTrue(contents.find('8 16 24')> 0)
        self.assertTrue(contents.find('DataArray type="Float64" NumberOfComponents="3" Name="points">')> 0)
        
    
    def test2(self):
        grid = datamodel. Grid.create([2,2,2], [1,1,1] | generic_unit_system.length)
        grid.rho = grid.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        data_file = StringIO()
        instance = vtk.VtkUnstructuredGrid("test.vtu", data_file, grid)
        instance.store()
        contents = data_file.getvalue()
        
        self.assertTrue(contents.find('<Piece NumberOfPoints="27" NumberOfCells="8"')> 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="connectivity">')> 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="offsets">')> 0)
        self.assertTrue(contents.find('11 11 11')> 0)
        self.assertTrue(contents.find('56 64')> 0)
        self.assertTrue(contents.find('DataArray type="Float64" NumberOfComponents="3" Name="points">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="rho">')>0)
        
    
    def test3(self):
        grid1 = datamodel. Grid.create([2,2,2], [1,1,1] | generic_unit_system.length)
        grid1.rho = grid1.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        grid2 = datamodel. Grid.create([2,2,2], [1,1,1] | generic_unit_system.length)
        grid2.position += [1.0,0.0,0.0] | generic_unit_system.length
        grid2.rho = grid2.x * (0.1 | generic_unit_system.mass / generic_unit_system.length ** 4)
        
        data_file = StringIO()
        instance = vtk.VtkUnstructuredGrid("test.vtu", data_file, [grid1, grid2])
        instance.store()
        contents = data_file.getvalue()
        index_of_piece_1 = contents.find('<Piece NumberOfPoints="27" NumberOfCells="8"')
        index_of_piece_2 = contents.find('<Piece NumberOfPoints="27" NumberOfCells="8"', index_of_piece_1 + 1)
        self.assertTrue(index_of_piece_1 > 0)
        self.assertTrue(index_of_piece_2 > 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="connectivity">')> 0)
        self.assertTrue(contents.find('DataArray type="Int32" NumberOfComponents="1" Name="offsets">')> 0)
        self.assertTrue(contents.find('11 11 11')> 0)
        self.assertTrue(contents.find('56 64')> 0)
        self.assertTrue(contents.find('DataArray type="Float64" NumberOfComponents="3" Name="points">')> 0)
        self.assertTrue(contents.find('<DataArray type="Float64" NumberOfComponents="1" Name="rho">')>0)
        
