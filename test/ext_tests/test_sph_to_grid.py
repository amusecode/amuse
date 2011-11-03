import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.support.exceptions import AmuseException
from amuse.ext.sph_to_grid import convert_SPH_to_grid
from amuse.units import units, generic_unit_system, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particles
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi

class TestSPH2Grid(TestWithMPI):
    
    def setup_sph_code(self, sph_code, number_of_particles, L, rho, u):
        converter = ConvertBetweenGenericAndSiUnits(L, rho, constants.G)
        sph_code = sph_code(converter, mode = 'periodic')#, redirection = 'none')
        sph_code.parameters.periodic_box_size = 10.0 | units.parsec
        gas = Particles(number_of_particles)
        gas.mass = (rho * L**3) / number_of_particles
        numpy.random.seed(12345)
        gas.x = L * numpy.random.uniform(0.0, 1.0, number_of_particles)
        gas.y = L * numpy.random.uniform(0.0, 1.0, number_of_particles)
        gas.z = L * numpy.random.uniform(0.0, 1.0, number_of_particles)
        gas.vx = numpy.zeros(number_of_particles) | units.cm / units.s
        gas.vy = numpy.zeros(number_of_particles) | units.cm / units.s
        gas.vz = numpy.zeros(number_of_particles) | units.cm / units.s
        gas.u = u
        if isinstance(sph_code, Fi):
            sph_code.parameters.self_gravity_flag = False
            sph_code.parameters.timestep = 0.1 | generic_unit_system.time
            gas.h_smooth = L / number_of_particles**(1/3.0)
            gas.position -= 0.5 * L
        
        sph_code.gas_particles.add_particles(gas)
        sph_code.commit_particles()
#~        sph_code.evolve_model(0.01 | units.Myr)
        return sph_code
    
    def test1(self):
        print "Testing convert_SPH_to_grid with Gadget2"
        number_of_particles = 10000
        L = 10.0 | units.parsec
        rho = 1.14 | units.amu/units.cm**3
        u = 5.e11 | units.cm**2 / units.s**2
        
        sph_code = self.setup_sph_code(Gadget2, number_of_particles, L, rho, u)
        grid = convert_SPH_to_grid(sph_code, (10,10,10))
        sph_code.stop()
        
        self.assertEqual(grid.shape, (10,10,10))
        self.assertAlmostEqual(grid.cellsize(), (L/10).as_vector_with_length(3))
        self.assertAlmostEqual(grid[0,0,0].position, (L/20).as_vector_with_length(3))
        self.assertEqual(
            grid.contains([[x,x,x] for x in numpy.linspace(-1.0, 11.0, 4)] | units.parsec),
            [False, True, True, False])
        self.assertAlmostEqual(grid.get_volume(), 1000.0 | units.parsec**3)
        
        mean_rho = grid.rho.mean()
        print "Mean density:", mean_rho.as_quantity_in(units.amu/units.cm**3), " -  original:", rho
        self.assertAlmostRelativeEqual(mean_rho, rho, 1)
        self.assertAlmostRelativeEqual(grid.energy.mean(), mean_rho * u, 1)
        self.assertEqual(grid.momentum, 0 * mean_rho * u.sqrt())
    
    def test2(self):
        print "Testing convert_SPH_to_grid with Fi"
        number_of_particles = 10000
        L = 10.0 | units.parsec
        rho = 1.14 | units.amu/units.cm**3
        u = 5.e11 | units.cm**2 / units.s**2
        
        sph_code = self.setup_sph_code(Fi, number_of_particles, L, rho, u)
        grid = convert_SPH_to_grid(sph_code, (10,10,10))
        sph_code.stop()
        
        self.assertEqual(grid.shape, (10,10,10))
        self.assertAlmostEqual(grid.cellsize(), (L/10).as_vector_with_length(3))
        self.assertAlmostEqual(grid[0,0,0].position, (L/20).as_vector_with_length(3))
        self.assertEqual(
            grid.contains([[x,x,x] for x in numpy.linspace(-1.0, 11.0, 4)] | units.parsec),
            [False, True, True, False])
        self.assertAlmostEqual(grid.get_volume(), 1000.0 | units.parsec**3)
        
        mean_rho = grid.rho.mean()
        print "Mean density:", mean_rho.as_quantity_in(units.amu/units.cm**3), " -  original:", rho
        self.assertAlmostRelativeEqual(mean_rho, rho, 1)
        self.assertAlmostRelativeEqual(grid.energy.mean(), mean_rho * u, 1)
        self.assertEqual(grid.momentum, 0 * mean_rho * u.sqrt())
    
    def test3(self):
        print "Testing convert_SPH_to_grid with Gadget2 and do_scale"
        number_of_particles = 10000
        L = 10.0 | units.parsec
        rho = 1.14 | units.amu/units.cm**3
        u = 5.e11 | units.cm**2 / units.s**2
        
        sph_code = self.setup_sph_code(Gadget2, number_of_particles, L, rho, u)
        grid = convert_SPH_to_grid(sph_code, (10,10,10), do_scale = True)
        sph_code.stop()
        
        self.assertEqual(grid.shape, (10,10,10))
        self.assertAlmostEqual(grid.cellsize(), (L/10).as_vector_with_length(3))
        self.assertAlmostEqual(grid[0,0,0].position, (L/20).as_vector_with_length(3))
        self.assertEqual(
            grid.contains([[x,x,x] for x in numpy.linspace(-1.0, 11.0, 4)] | units.parsec),
            [False, True, True, False])
        self.assertAlmostEqual(grid.get_volume(), 1000.0 | units.parsec**3)
        
        mean_rho = grid.rho.mean()
        print "Mean density:", mean_rho.as_quantity_in(units.amu/units.cm**3), " -  original:", rho
        self.assertAlmostRelativeEqual(mean_rho, rho, 7)
        self.assertAlmostRelativeEqual(grid.energy.mean(), mean_rho * u, 7)
        self.assertEqual(grid.momentum, 0 * mean_rho * u.sqrt())
    
    def test4(self):
        print "Testing convert_SPH_to_grid with Fi and do_scale"
        number_of_particles = 10000
        L = 10.0 | units.parsec
        rho = 1.14 | units.amu/units.cm**3
        u = 5.e11 | units.cm**2 / units.s**2
        
        sph_code = self.setup_sph_code(Fi, number_of_particles, L, rho, u)
        grid = convert_SPH_to_grid(sph_code, (10,10,10), do_scale = True)
        sph_code.stop()
        
        self.assertEqual(grid.shape, (10,10,10))
        self.assertAlmostEqual(grid.cellsize(), (L/10).as_vector_with_length(3))
        self.assertAlmostEqual(grid[0,0,0].position, (L/20).as_vector_with_length(3))
        self.assertEqual(
            grid.contains([[x,x,x] for x in numpy.linspace(-1.0, 11.0, 4)] | units.parsec),
            [False, True, True, False])
        self.assertAlmostEqual(grid.get_volume(), 1000.0 | units.parsec**3)
        
        mean_rho = grid.rho.mean()
        print "Mean density:", mean_rho.as_quantity_in(units.amu/units.cm**3), " -  original:", rho
        self.assertAlmostRelativeEqual(mean_rho, rho, 7)
        self.assertAlmostRelativeEqual(grid.energy.mean(), mean_rho * u, 7)
        self.assertEqual(grid.momentum, 0 * mean_rho * u.sqrt())
    
    def test5(self):
        print "Testing exceptions"
        class BogusHydroCode(object):
            MODE_PERIODIC_BOUNDARIES = "periodic"
            def __init__(self, mode = "normal"):
                self.mode = mode
        
        self.assertRaises(AmuseException, convert_SPH_to_grid, BogusHydroCode(), (10,10,10),
            expected_message = "Only periodic boundary conditions supported")
        
        self.assertRaises(AmuseException, convert_SPH_to_grid, BogusHydroCode(mode = 'periodic'), (10,10),
            expected_message = "Argument dimensions must contain exactly three numbers")
        
        self.assertRaises(AmuseException, convert_SPH_to_grid, BogusHydroCode(mode = 'periodic'), (10,10,10),
            expected_message = "Unknown hydrodynamics code: BogusHydroCode - don't know whether the "
            "box runs from 0 to L or from -0.5 L to 0.5 L.")
    
