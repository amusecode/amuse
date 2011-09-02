import os.path
import numpy
from amuse.test.amusetest import get_path_to_results, TestWithMPI
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, Particle, ParticlesSuperset, Grid


from amuse.support.exceptions import AmuseException
from amuse.ext.grid_to_sph import Grid2SPH, convert_grid_to_SPH
from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import nbody_system
from amuse.units import constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
class TestGrid2SPH(TestWithMPI):
    
    def setup_simple_grid(self):
        test_grid = Grid.create((4,3,2), [1.0, 1.0, 1.0] | units.m)
        test_grid.rho = numpy.linspace(1.0, 2.0, num=24).reshape(test_grid.shape) | units.kg/units.m**3
        test_grid.rhovx = test_grid.rho * (3.0 | units.m/units.s)
        test_grid.rhovy = test_grid.rho * (4.0 | units.m/units.s)
        test_grid.rhovz = test_grid.rho * (0.0 | units.m/units.s)
        test_grid.energy = test_grid.rho * (1.0 | (units.m/units.s)**2)
        return test_grid
    
    def test0(self):
        print "Testing the simple example grid"
        test_grid = self.setup_simple_grid()
        self.assertEqual(test_grid.position[ 0][ 0][ 0],  [1.0/8.0, 1.0/6.0, 1.0/4.0] | units.m)
        self.assertEqual(test_grid.position[-1][-1][-1],  [7.0/8.0, 5.0/6.0, 3.0/4.0] | units.m)
        self.assertEqual(test_grid.momentum[ 0][ 0][ 0],  [3.0, 4.0, 0.0] | (units.kg/units.m**3) * (units.m/units.s))
        self.assertEqual(test_grid.momentum[-1][-1][-1],  [6.0, 8.0, 0.0] | (units.kg/units.m**3) * (units.m/units.s))
        self.assertEqual(test_grid.energy[ 0][ 0][ 0],  1.0 | (units.J/units.m**3))
        self.assertEqual(test_grid.energy[-1][-1][-1],  2.0 | (units.J/units.m**3))
        
    def test1(self):
        print "Testing the converter"
        number_of_particles = 10000
        test_grid = self.setup_simple_grid()
        converter = Grid2SPH(test_grid, number_of_particles)
        self.assertTrue(converter.grid is test_grid)
        self.assertEqual(converter.shape, (4,3,2))
        self.assertEqual(converter.number_of_sph_particles, number_of_particles)
        self.assertEqual(converter.base_distribution_type, "uniform")
        
        converter.setup_lookup_tables()
        converter.setup_variates()
        self.assertEqual(converter.cumulative_weight[0],  1.0/(1.5*4*3*2))
        self.assertEqual(converter.cumulative_weight[-1], 1.0)
        self.assertEqual(converter.position_lookup_table[0],  [1.0/8.0, 1.0/6.0, 1.0/4.0] | units.m)
        self.assertEqual(converter.position_lookup_table[-1], [7.0/8.0, 5.0/6.0, 3.0/4.0] | units.m)
        self.assertEqual(converter.position_lookup_table[9],  [3.0/8.0, 3.0/6.0, 3.0/4.0] | units.m)
        self.assertAlmostEqual(converter.velocity_lookup_table,  [3.0, 4.0, 0.0] | units.m/units.s)
        self.assertEqual(converter.specific_internal_energy_lookup_table,  1.0 | units.J/units.kg)
        self.assertEqual(converter.cellsize_unit,  units.m)
        self.assertTrue(converter.cellsize_unit is units.m)
        self.assertAlmostEqual(converter.cellsize_number, [0.25, 1/3.0, 0.5])
        self.assertAlmostEqual(converter.mass, 1.5 | units.kg)
        # The number of particles in a cell should scale with the amount of mass in the cell:
        self.assertAlmostRelativeEqual(
            converter.mass * numpy.histogram(converter.indices, bins=4*3*2)[0] * 1.0/number_of_particles, 
            test_grid.rho.flatten()*test_grid.cellsize().prod(), 
            places = 2
        )
    
    def test2(self):
        print "Testing the user interface"
        number_of_particles = 10000
        test_grid = self.setup_simple_grid()
        sph_particles = convert_grid_to_SPH(test_grid, number_of_particles)
        self.assertEqual(len(sph_particles), number_of_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(),  1.5 | units.kg)
        self.assertAlmostEqual(sph_particles.velocity,  [3.0, 4.0, 0.0] | units.m/units.s)
        self.assertAlmostEqual(sph_particles.u,  1.0 | (units.m/units.s)**2)
        # The number of particles in a cell should scale with the amount of mass in the cell:
        self.assertAlmostRelativeEqual(
            (1.5 | units.kg)/number_of_particles * numpy.histogramdd(
                sph_particles.position.value_in(units.m), bins=(4,3,2))[0], 
            test_grid.rho*test_grid.cellsize().prod(), 
            places = 2
        )
        self.assertAlmostEqual(sph_particles.h_smooth,  (50.0/number_of_particles)**(1.0/3) | units.m)
    
    def test3(self):
        print "Testing the user interface, random base_distribution_type"
        number_of_particles = 10000
        test_grid = self.setup_simple_grid()
        sph_particles = convert_grid_to_SPH(test_grid, number_of_particles, 
            base_distribution_type = "random", seed = 12345)
        self.assertEqual(len(sph_particles), number_of_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(),  1.5 | units.kg)
        self.assertAlmostEqual(sph_particles.velocity,  [3.0, 4.0, 0.0] | units.m/units.s)
        self.assertAlmostEqual(sph_particles.u,  1.0 | (units.m/units.s)**2)
        # For 'random', the number of particles in a cell should scale only on average 
        # with the amount of mass in the cell:
        self.assertAlmostRelativeEqual(
            ((1.5 | units.kg)/number_of_particles * numpy.histogramdd(
                sph_particles.position.value_in(units.m), bins=(4,3,2))[0]).sum(), 
            (test_grid.rho*test_grid.cellsize().prod()).sum(), 
            places = 2
        )
        self.assertRaises(AssertionError, 
            self.assertAlmostRelativeEqual,
                (1.5 | units.kg)/number_of_particles * numpy.histogramdd(sph_particles.position.value_in(units.m), bins=(4,3,2))[0], 
                test_grid.rho*test_grid.cellsize().prod(), 
                places = 2,
        )
        self.assertAlmostEqual(sph_particles.h_smooth,  (50.0/number_of_particles)**(1.0/3) | units.m)
    
    def test4(self):
        print "Testing exceptions"
        number_of_particles = 10000
        test_grid = self.setup_simple_grid()
        self.assertEqual(test_grid[0].number_of_dimensions(), 2)
        self.assertRaises(AmuseException, convert_grid_to_SPH, test_grid[0], number_of_particles,
            expected_message = "Grid must be 3D")
        self.assertRaises(AmuseException, convert_grid_to_SPH, test_grid, 
            number_of_particles, base_distribution_type = "bogus",
            expected_message = "Unknown base_distribution_type: bogus. Possible "
                "options are: 'random' or 'uniform'.")
    
