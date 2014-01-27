
from amuse.test import amusetest
from amuse.units import units, quantities
from amuse.datamodel.particles import Particles

from amuse.ext.stellar_wind import StarsWithSimpleWind, StarsWithAcceleratingWind


class TestStellarWind(amusetest.TestCase):

    def assertGreaterEqual(self, value, expected):
        self.assertTrue(value >= expected, "Expected {0} >= {1}".format(value, expected))
    
    def assertLessEqual(self, value, expected):
        self.assertTrue(value <= expected, "Expected {0} <= {1}".format(value, expected))
        
    def test1(self):
        """ Test the particles created for a simple wind """
        star = self.create_star()
        star_wind = StarsWithSimpleWind(star, 1e-8|units.MSun)

        star_wind.evolve_model(1|units.yr)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 100)
        min_dist = (wind.position - star.position).lengths().min()
        max_dist = (wind.position - star.position).lengths().max()
        self.assertGreaterEqual(min_dist, 2|units.RSun)
        self.assertLessEqual(max_dist, 24.7|units.RSun)
        plus_x_wind = wind[wind.x > star.x]
        minus_x_wind = wind[wind.x < star.x]
        plus_y_wind = wind[wind.y > star.y]
        minus_y_wind = wind[wind.y < star.y]
        plus_z_wind = wind[wind.z > star.z]
        minus_z_wind = wind[wind.z < star.z]
        self.assertGreaterEqual(minus_x_wind.vx.min(), -1500|units.ms)
        self.assertLessEqual(minus_x_wind.vx.max(), -1000|units.ms)
        self.assertGreaterEqual(plus_x_wind.vx.min(), -1000|units.ms)
        self.assertLessEqual(plus_x_wind.vx.max(), -500|units.ms)
        self.assertGreaterEqual(minus_y_wind.vy.min(), -500|units.ms)
        self.assertLessEqual(minus_y_wind.vy.max(), 0|units.ms)
        self.assertGreaterEqual(plus_y_wind.vy.min(), 0|units.ms)
        self.assertLessEqual(plus_y_wind.vy.max(), 500|units.ms)
        self.assertGreaterEqual(minus_z_wind.vz.min(), 500|units.ms)
        self.assertLessEqual(minus_z_wind.vz.max(), 1000|units.ms)
        self.assertGreaterEqual(plus_z_wind.vz.min(), 1000|units.ms)
        self.assertLessEqual(plus_z_wind.vz.max(), 1500|units.ms)

    def test2(self):
        """ Test the accelerating wind """
        star = self.create_star()
        star_wind = StarsWithAcceleratingWind(star, 1e-8|units.MSun, init_v_wind_ratio=0.1)

        star_wind.evolve_model(1|units.yr)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 100)
        min_dist = (wind.position - star.position).lengths().min()
        max_dist = (wind.position - star.position).lengths().max()
        self.assertGreaterEqual(min_dist, 2|units.RSun)
        self.assertLessEqual(max_dist, 24.7|units.RSun)
        plus_x_wind = wind[wind.x > star.x]
        minus_x_wind = wind[wind.x < star.x]
        plus_y_wind = wind[wind.y > star.y]
        minus_y_wind = wind[wind.y < star.y]
        plus_z_wind = wind[wind.z > star.z]
        minus_z_wind = wind[wind.z < star.z]
        self.assertGreaterEqual(minus_x_wind.vx.min(), -1050|units.ms)
        self.assertLessEqual(minus_x_wind.vx.max(), -1000|units.ms)
        self.assertGreaterEqual(plus_x_wind.vx.min(), -1000|units.ms)
        self.assertLessEqual(plus_x_wind.vx.max(), -950|units.ms)
        self.assertGreaterEqual(minus_y_wind.vy.min(), -50|units.ms)
        self.assertLessEqual(minus_y_wind.vy.max(), 0|units.ms)
        self.assertGreaterEqual(plus_y_wind.vy.min(), 0|units.ms)
        self.assertLessEqual(plus_y_wind.vy.max(), 50|units.ms)
        self.assertGreaterEqual(minus_z_wind.vz.min(), 950|units.ms)
        self.assertLessEqual(minus_z_wind.vz.max(), 1000|units.ms)
        self.assertGreaterEqual(plus_z_wind.vz.min(), 1000|units.ms)
        self.assertLessEqual(plus_z_wind.vz.max(), 1050|units.ms)

    def test3(self):
        """ Test the wind acceleration """
        star = self.create_star()
        star.position = [1, 1, 1] | units.RSun
        star_wind = StarsWithAcceleratingWind(star, 1e-8|units.MSun)

        # unaffected points
        points = [[1,1,1], [1,1,0], [1,4,1], [12,1,1]]|units.RSun
        x, y, z = points.transpose()
        ax, ay, az = star_wind.get_gravity_at_point(1, x, y, z)
        zeros = [0,0,0,0]|units.m/units.s**2
        self.assertEqual(ax, zeros)
        self.assertEqual(ay, zeros)
        self.assertEqual(az, zeros)

        # affected points
        points = [[5.1,1,1], [1,1,5.1], [1,7,1], [1,7,8], [-8,1,1]]|units.RSun
        x, y, z = points.transpose()
        ax, ay, az = star_wind.get_gravity_at_point(1, x, y, z)
        a = quantities.as_vector_quantity([ax, ay, az]).transpose()
        self.assertAlmostEquals(a[0], [5.987e-5, 0, 0]|units.m/units.s**2, places=8)
        self.assertAlmostEquals(a[1], [0, 0, 5.987e-5]|units.m/units.s**2, places=8)
        self.assertAlmostEquals(a[2], [0, 2.796e-5, 0]|units.m/units.s**2, places=8)
        self.assertAlmostEquals(a[3], [0, 7.706e-6, 8.990e-6]|units.m/units.s**2, places=8)
        self.assertAlmostEquals(a[4], [-1.243e-5, 0, 0]|units.m/units.s**2, places=8)

    def create_star(self):
        star = Particles(1)
        star.mass = 2|units.MSun
        star.radius = 2|units.RSun
        star.temperature = 5000|units.K
        star.position = [1, 1, 1] | units.parsec
        star.velocity = [-1000, 0, 1000] | units.ms
        star.wind_mass_loss = 1e-6 | units.MSun / units.yr
        star.terminal_wind_velocity = 500 | units.ms
        return star
