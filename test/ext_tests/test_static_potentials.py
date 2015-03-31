
from amuse.test import amusetest
from amuse.units import units, quantities

from amuse.ext import static_potentials


class TestStaticPotentials(amusetest.TestCase):

    def test1(self):
        """ Test the Galactic potential """
        static_potential = static_potentials.Disc_Bulge_Halo_Potential()

        points = [[1,0,0],
                  [10,0,0],
                  [1e3,0,1e4],
                  ] | units.kpc
        potentials = static_potential.get_potential_at_point(0, *points.transpose()).transpose()

        self.assertAlmostEqual(potentials[0], 96648.29247|units.kms**2, 6)
        self.assertAlmostEqual(potentials[1], -8827.40917774|units.kms**2, 6)
        self.assertAlmostEqual(potentials[2], -266129.25389|units.kms**2, 6)

    def test2(self):
        """ Test the Galactic center potential """
        static_potential = static_potentials.Galactic_Center_Potential_Kruijssen()

        points = [[1,0,0],
                  [0,10,0],
                  [1e3,0,1e2],
                  ] | units.parsec
        potentials = static_potential.get_potential_at_point(0, *points.transpose()).transpose()
        print potentials.in_(units.kms**2)
        self.assertAlmostEqual(potentials[0], 1265.5236737|units.kms**2, 6)
        self.assertAlmostEqual(potentials[1], -683.13285637|units.kms**2, 6)
        self.assertAlmostEqual(potentials[2], -528906.185054|units.kms**2, 6)

