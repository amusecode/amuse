import numpy

from amuse.test import amusetest
from amuse.units import units

from amuse.ext.roche_radius import Roche_Orbit, sepinsky_formula

class TestRocheRadius(amusetest.TestCase):

    def test1(self):
        """ Test that the basic numbers from Sepinsky are not unreasonable. """
        self.assertAlmostEqual(sepinsky_formula(), 1.0, places=1)

        q_values = 10.0**numpy.arange(-8, 9)
        answers = [1.0] * len(q_values)
        self.assertAlmostEqual(sepinsky_formula(q=q_values), answers, places=1)

        a_values = 10.0**numpy.arange(-4, 2)
        answers = [1.0]*len(a_values)
        self.assertIsOfOrder(sepinsky_formula(q=1.0, A=a_values), answers)

        self.assertIsOfOrder(sepinsky_formula(q=1e-4, A=a_values), answers)

        self.assertIsOfOrder(sepinsky_formula(q=1e4, A=a_values), answers)


    def test2(self):
        """ Test the Roche_Orbit class """

        roche_orbit = Roche_Orbit()

        self.assertAlmostEquals(roche_orbit.eggleton_roche_radius(), 81.5 | units.RSun, places=2)
        self.assertAlmostEquals(roche_orbit.sepinsky_roche_radius(), 81.32 | units.RSun, places=2)

        roche_orbit.eccentricity = 0.5

        self.assertAlmostEquals(roche_orbit.eggleton_roche_radius(), 40.75 | units.RSun, places=2)
        self.assertAlmostEquals(roche_orbit.sepinsky_roche_radius(), 39.62 | units.RSun, places=2)

        roche_orbit.semimajor_axis = 20.0 | units.AU
        roche_orbit.mass_2 = 10000.0 | units.MSun
        roche_orbit.eccentricity = 0.95
        roche_orbit.angular_velocity_ratio = 0.1

        self.assertAlmostEquals(roche_orbit.eggleton_roche_radius(), 4.87 | units.RSun, places=2)
        self.assertAlmostEquals(roche_orbit.sepinsky_roche_radius(), 5.39 | units.RSun, places=2)

        roche_orbit.eccentricity = numpy.array([0.0, 0.5, 0.95])
        self.assertAlmostEquals(roche_orbit.sepinsky_roche_radius(), [107.92, 53.93, 5.39] | units.RSun, places=2)
