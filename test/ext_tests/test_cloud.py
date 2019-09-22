import sys
import os
import numpy

from amuse.test import amusetest


from amuse.ext import cloud

from amuse.units import units
from amuse.units import generic_unit_system
from amuse import datamodel
class CloudTests(amusetest.TestCase):

    def test1(self):
        grid = datamodel.Grid.create((10,10,10), [10.0, 10.0, 10.0] | units.m)
        
        grid.rho = 0.0 | units.kg / units.m**3
        grid.rhovx = 0.0 | units.kg / (units.s * units.m**2)
        grid.rhovy = 0.0 | units.kg / (units.s * units.m**2)
        grid.rhovz = 0.0 | units.kg / (units.s * units.m**2)
        grid.energy = 1.0 | units.kg / ((units.s**2) * units.m)
        
        datamodel.Grid.add_global_vector_attribute("position", ["x","y","z"])
        
        cloud.fill_grid_with_spherical_cloud(
            grid, 
            center = [5.0, 5.0, 5.0] | units.m,
            radius = 2.0 | units.m,
            rho = 1.0 | units.kg / units.m**3,
            rhovx = 0.0 | units.kg / (units.s * units.m**2),
            rhovy = 0.1 | units.kg / (units.s * units.m**2),
            rhovz = 0.0 | units.kg / (units.s * units.m**2),
            energy = 1.0 | units.kg / ((units.s**2) * units.m)
        )
        
        self.assertEqual(grid.shape, (10,10,10))
        self.assertEqual(grid.rho[5][5][5], 1.0 | units.kg / units.m**3)
        
        #from mpl_toolkits.axes_grid1 import ImageGrid
        #from matplotlib import pyplot
        #figure = pyplot.figure()
        #grids = ImageGrid(figure, 111, nrows_ncols = (2, 2), axes_pad=0.1)
        #for i in range(4):
        #    z = grid.rho[i+2].value_in( units.kg / units.m**3,)
        #    grids[i].imshow(z)
        #figure.savefig('ax.png')
        
        self.assertEqual(grid.rho[5][6][5], 0.828125| units.kg / units.m**3)
        self.assertEqual(grid.rho[5][3][5], 0.828125| units.kg / units.m**3)
        self.assertEqual(grid.rho[5][5][6], 0.828125| units.kg / units.m**3)
        self.assertEqual(grid.rho[5][5][3], 0.828125| units.kg / units.m**3)


    def test2(self):
        n = 10
        grid = datamodel.Grid.create((n,n,n), [10.0, 10.0, 10.0] | generic_unit_system.length)
        datamodel.Grid.add_global_vector_attribute("position", ["x","y","z"])
        
        cloud.fill_grid_with_cloud_shock(
            grid, 
            center = [5.0, 5.0, 5.0] | generic_unit_system.length,
            radius = 1.0 | generic_unit_system.length,
        )
        
        self.assertAlmostRelativeEquals(grid.rho[5][4][4], 0.5921875 | generic_unit_system.mass / generic_unit_system.length**3)
    
    
