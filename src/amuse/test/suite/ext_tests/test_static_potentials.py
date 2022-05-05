import numpy
from amuse.test import amusetest
from amuse.units import units, quantities

from amuse.ext import static_potentials

from amuse.support.console import set_printing_strategy

class TestGalacticCenterPotential(amusetest.TestCase):
    def setUp(self):
        set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.Myr, units.kms, units.kms**2, units.m/units.s**2])

    def tearDown(self):
        set_printing_strategy("default")
        
    def test1(self):
        """ test enclosed mass profile of Kruijssen potential """

        static_potential = static_potentials.Galactic_Center_Potential_Kruijssen()

        m_zero = static_potential.enclosed_mass(0|units.parsec)
        self.assertAlmostEqual(m_zero, 0|units.MSun)

        factor = 8.3/8.5
        m_on_bin = static_potential.enclosed_mass(7.489|units.parsec * factor)
        self.assertAlmostEqual(m_on_bin, factor**2 * 17170000.|units.MSun)

        m_in_bin = static_potential.enclosed_mass(200|units.parsec)
        self.assertAlmostEqual(m_in_bin, 1.644965507|1E9*units.MSun)

        # test array
        r = numpy.logspace(1, 3, 500) | units.parsec
        m = static_potential.enclosed_mass(r)
        self.assertAlmostEqual(m[10], 2.0680055041|1E7*units.MSun)


        # test 2D array
        x = numpy.logspace(1, 3, 500) | units.parsec
        matrix, m = quantities.meshgrid(x, x)
        mm = static_potential.enclosed_mass(matrix)
        self.assertAlmostEqual(mm[20,20], 2.223803269|1E7*units.MSun)

        # quick plot for testing
        # from matplotlib import pyplot
        # from amuse import plot as aplot
        # aplot.loglog(r, m, marker="*")
        # aplot.scatter(static_potential.radius, static_potential.enclosed_mass_profile)
        # pyplot.show()

    def test2(self):
        """ Test get_potential_at_point """
        static_potential = static_potentials.Galactic_Center_Potential_Kruijssen()

        points = [[1,0,0],
                  [0,10,0],
                  [1e3,0,1e2],
                  ] | units.parsec
        potentials = static_potential.get_potential_at_point(0, *points.transpose())
        # print potentials.in_(units.kms**2)
        self.assertAlmostEqual(potentials[0], -17280.958078|units.kms**2, 6)
        self.assertAlmostEqual(potentials[1], -8296.5300805|units.kms**2, 6)
        self.assertAlmostEqual(potentials[2], -14965.2205261|units.kms**2, 6)

    def test3(self):
        """ Test get_gravity_at_point """
        static_potential = static_potentials.Galactic_Center_Potential_Kruijssen()

        points = [[1,0,0],
                  [0,10,0],
                  [1e3,0,1e2],
                  ] | units.parsec
        grav = static_potential.get_gravity_at_point(0, *points.transpose()).transpose()
        print(grav)
        self.assertAlmostEqual(grav[0], [-5.6003771, 0, 0]|1e-7*units.m/units.s**2, 6)
        self.assertAlmostEqual(grav[1], [0, -2.6887222, 0]|1e-8*units.m/units.s**2, 6)
        self.assertAlmostEqual(grav[2], [-4.7661597, 0, -1.200846]|1e-10*units.m/units.s**2, 6)


class TestGalacticPotential(amusetest.TestCase):
    def test1(self):
        """ Test the Galactic potential """
        static_potential = static_potentials.Disc_Bulge_Halo_Potential()

        points = [[1,0,0],
                  [10,0,0],
                  [1e3,0,1e4],
                  ] | units.kpc
        potentials = static_potential.get_potential_at_point(0, *points.transpose()).transpose()

        self.assertAlmostEqual(potentials[0], -96648.29247|units.kms**2, 6)
        self.assertAlmostEqual(potentials[1], 8827.40917774|units.kms**2, 6)
        self.assertAlmostEqual(potentials[2], 266129.25389|units.kms**2, 6)

