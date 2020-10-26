import numpy

from amuse.test import amusetest
from amuse.units import units, quantities
from amuse.datamodel.particles import Particles, Particle
# from amuse.community.sse.interface import SSE
from amuse.community.seba.interface import SeBa

from amuse.ext import stellar_wind

from amuse.support.console import set_preferred_units


from amuse.support.console import set_printing_strategy


class TestStellarWind(amusetest.TestCase):
    def setUp(self):
        set_preferred_units(units.MSun, units.RSun, units.yr, units.kms,
                            units.kms**2, units.m/units.s**2, units.W, units.erg)

    def skip_no_scipy(self):
        try:
            import scipy
        except ImportError:
            self.skip("Failed to import scipy")

    def assertAllAlmostEqual(self, values, expected, *args, **kwargs):
        for value in values:
            self.assertAlmostEqual(value, expected, *args, **kwargs)

    def assertDecreasing(self, array):
        delta = array[1:] - array[:-1]
        for d in delta:
            self.assertLess(d, quantities.zero)

    def create_star(self, N=1):
        star = Particles(N)
        star.mass = 2 | units.MSun
        star.radius = 2 | units.RSun
        star.temperature = 5000 | units.K
        star.position = [1, 1, 1] | units.parsec
        star.velocity = [-1000, 0, 1000] | units.ms
        star.wind_mass_loss_rate = 1e-6 | units.MSun / units.yr
        star.initial_wind_velocity = 50 | units.ms
        star.terminal_wind_velocity = 500 | units.ms
        return star

    def create_stars_without_wind_attributes(self):
        stars = Particles(2)
        stars.mass = 2 | units.MSun
        stars.radius = 2 | units.RSun
        stars.temperature = 5000 | units.K
        stars.luminosity = 2 | units.LSun
        stars.age = 0 | units.Myr
        stars[0].position = [1, 1, 1] | units.parsec
        stars[1].position = [-1, -1, -1] | units.parsec
        stars[0].velocity = [-1000, 0, 1000] | units.ms
        stars[1].velocity = [0, 0, 0] | units.ms

        return stars


    def tearDown(self):
        set_printing_strategy('default')


class TestStarsWithMassLoss(TestStellarWind):
    def test_add_particles(self):
        star_particles = self.create_star(2)
        print("created parts")
        stars = stellar_wind.StarsWithMassLoss()
        print("created stars")
        stars.add_particles(star_particles)
        print("added", stars)

        self.assertEqual(stars.lost_mass, [0, 0] | units.MSun)
        self.assertEqual(stars.wind_release_time, [0, 0] | units.yr)

        self.assertAlmostEqual(stars.mu[0], 1.145934546 | units.amu)

    def test_add_particle(self):
        star_particles = self.create_star(1)
        stars = stellar_wind.StarsWithMassLoss()
        star = stars.add_particle(star_particles[0])

        self.assertEqual(star.lost_mass, 0 | units.MSun)
        self.assertEqual(star.wind_release_time, 0 | units.yr)

        self.assertAlmostEqual(star.mu, 1.145934546 | units.amu)

        attrs = stars.get_attribute_names_defined_in_store()
        self.assertFalse("mechanical_energy" in attrs)
        self.assertFalse("previous_mechanical_luminosity" in attrs)

    def test_default_lmech(self):
        star_particles = self.create_star(1)
        stars = stellar_wind.StarsWithMassLoss()
        stars.track_mechanical_energy()
        star = stars.add_particle(star_particles[0])

        self.assertEqual(star.mechanical_energy, 0 | units.J)

    def test_add_particle_later(self):
        star_particles = self.create_star(2)
        stars = stellar_wind.StarsWithMassLoss()
        stars.add_particle(star_particles[0])

        stars.evolve_mass_loss(1 | units.yr)

        new_star = stars.add_particle(star_particles[1])
        self.assertAlmostEqual(new_star.wind_release_time, 1 | units.yr)

    def test_superfluous_attributes(self):
        star_particles = self.create_star(2)
        star_particles.other_attr = 10 | units.yr

        stars = stellar_wind.StarsWithMassLoss()
        stars.add_particles(star_particles)

        self.assertRaises(AttributeError, stars.__getattr__, "other_attr")

    def test_add_wrong_attribute(self):
        star_particles = self.create_star(2)
        star_particles

        stars = stellar_wind.StarsWithMassLoss()
        stars.add_particles(star_particles)

        self.assertRaises(AttributeError, stars.__setattr__, "other_attr", 10)


class FakeStevCode(object):
    def __init__(self, *args, **kwargs):
        self.particles = Particles()
        self.model_time = 0 | units.Myr
        self.mass_loss_rate = 1e-11 | units.MSun/units.yr

    def evolve_model(self, time):
        dt = time - self.model_time
        self.model_time = time
        self.particles.mass -= self.mass_loss_rate * dt
        self.particles.age = time


class TestPositionGenerator(TestStellarWind):
    def test_random_cube(self):
        numpy.random.seed(1232367)
        gen = stellar_wind.PositionGenerator('random')
        cube = gen.cube_generator(1000)
        self.assertEqual(cube.shape, (1000, 3))
        self.assertAlmostEqual(cube.min(), -0.99970440813)
        self.assertAlmostEqual(cube.max(), 0.997978975085)

    def test_regular_cube(self):
        numpy.random.seed(1232367)
        gen = stellar_wind.PositionGenerator(grid_type="regular")
        cube = gen.cube_generator(25)
        print(cube)
        self.assertEqual(cube.shape, (25, 3))
        self.assertAlmostEqual(cube.min(), -2./3.)
        self.assertAlmostEqual(cube.max(), 2./3.)

    def test_cutout_sphere(self):
        points = [[0., 0., 0.],
                  [1., 1., 1.],
                  [-1., -1., -1.],
                  [-0., 0., -.99],
                  [0., .5, 0.],
                  ]
        points = numpy.array(points)

        gen = stellar_wind.PositionGenerator()
        remaining = gen.cutout_sphere(points, 0.1)

        self.assertEqual(len(remaining), 2)
        self.assertEqual(remaining, points[-2:])

    def test_uniform_hollow_sphere(self):
        numpy.random.seed(1232367)
        gen = stellar_wind.PositionGenerator()
        points = gen.uniform_hollow_sphere(2, 0.6)
        self.assertEqual(len(points), 2)

        points = gen.uniform_hollow_sphere(1000, 0.1)
        self.assertEqual(len(points), 1000)

    def test_regular_hollow_sphere(self):
        numpy.random.seed(1232367)
        gen = stellar_wind.PositionGenerator(grid_type="regular")
        points = gen.uniform_hollow_sphere(2, 0.6)
        self.assertEqual(len(points), 2)

        points = gen.uniform_hollow_sphere(1000, 0.1)
        self.assertEqual(len(points), 1000)

    def test_density_distribution(self):
        numpy.random.seed(1234567)
        N = 100000
        rmin = 1 | units.RSun
        rmax = 10 | units.RSun
        gen = stellar_wind.PositionGenerator()
        p, _ = gen.generate_positions(N, rmin, rmax)
        r = p.lengths()

        print('rmin', r.min())
        print('rmax', r.max())
        print(r.mean())
        self.assertEqual(len(r), N)
        self.assertGreaterEqual(r.min(), rmin)
        self.assertLessEqual(r.max(), rmax)
        self.assertAlmostEqual(r.mean(), 5.49955427602 | units.RSun)

        return

        r = r.value_in(units.RSun)
        n_bins = 50
        n, bins = numpy.histogram(r, n_bins, range=(1, 10))
        bin_volume = 4./3. * numpy.pi * (bins[1:]**3 - bins[:-1]**3)
        dens = n / bin_volume
        bin_means = (bins[1:] + bins[:-1])/2.

        s = p[abs(p[:, 2]) < (0.1 | units.RSun)]
        from matplotlib import pyplot
        from amuse import plot as aplot
        aplot.plot(s[:, 0], s[:, 1], '.')
        pyplot.axis('equal')
        pyplot.savefig("scatter.pdf")
        pyplot.clf()

        x = numpy.linspace(1, 10, num=200)
        y = 300. / x**2
        # y = 0. * x + dens.mean()

        from matplotlib import pyplot
        # pyplot.plot(x, y)
        pyplot.loglog(x, y)
        pyplot.plot(bin_means, dens, '*')
        pyplot.savefig("dens.pdf")

    def test_alternate_radius_function(self):
        numpy.random.seed(123)
        gen = stellar_wind.PositionGenerator()
        # func = stellar_wind.LogisticVelocityAcceleration().radius_from_number
        func = stellar_wind.ConstantVelocityAcceleration().radius_from_number

        N = 100000
        rmin = 2 | units.RSun
        rmax = 6 | units.RSun
        star = Particle()
        star.radius = rmin
        star.acc_cutoff = 10 | units.RSun
        star.initial_wind_velocity = 4 | units.kms
        star.terminal_wind_velocity = 12 | units.kms

        p, _ = gen.generate_positions(N, rmin, rmax, func, star)
        r = p.lengths()

        self.assertEqual(len(r), N)
        self.assertGreaterEqual(r.min(), rmin)
        self.assertLessEqual(r.max(), rmax)
        print(r.mean())
        self.assertAlmostEqual(r.mean(), 4.00196447056 | units.RSun)

        return

        r = r.value_in(units.RSun)
        n_bins = 50
        n, bins = numpy.histogram(r, n_bins)
        bin_volume = 4./3. * numpy.pi * (bins[1:]**3 - bins[:-1]**3)
        dens = n / bin_volume
        bin_means = (bins[1:] + bins[:-1])/2.

        from matplotlib import pyplot
        pyplot.plot(bin_means, dens, '*')
        pyplot.savefig("dens.pdf")

    def test_random_rotation(self):
        numpy.random.seed(1232367)
        gen = stellar_wind.PositionGenerator()
        axis, angle = gen.random_rotation()
        print(axis)
        print(angle)
        self.assertEqual(len(axis), 3)
        self.assertAlmostEqual((axis**2).sum(), 1.)
        self.assertGreaterEqual(angle, 0.)
        self.assertLessEqual(angle, 2. * numpy.pi)

    def test_rotation_matrix(self):
        gen = stellar_wind.PositionGenerator()
        axis = [0, 0, 1]
        angle = numpy.pi
        matrix = gen.rotation_matrix(axis, angle)
        print(numpy.rint(matrix))

        expected = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
        self.assertAlmostEqual(matrix, expected)

    def test_rotate_positions(self):
        gen = stellar_wind.PositionGenerator()
        axis = [0, 1, 0]
        angle = numpy.pi
        positions = [[0.5, 0.5, 0.5], [0.2, -0.5, 1.3]]
        rotated = gen.rotate_positions(positions, axis, angle)
        print(rotated)
        expected = [[-0.5, 0.5, -0.5], [-0.2, -0.5, -1.3]]
        self.assertAlmostEqual(rotated, expected)


class TestSimpleWind(TestStellarWind):
    def test_wind_particles(self):
        star = self.create_star()
        star_wind = stellar_wind.new_stellar_wind(1e-8 | units.MSun)
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(1 | units.yr)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 100)
        min_dist = (wind.position - star.position).lengths().min()
        max_dist = (wind.position - star.position).lengths().max()
        self.assertGreaterEqual(min_dist, 2 | units.RSun)
        self.assertLessEqual(max_dist, 24.7 | units.RSun)
        plus_x_wind = wind[wind.x > star.x]
        minus_x_wind = wind[wind.x < star.x]
        plus_y_wind = wind[wind.y > star.y]
        minus_y_wind = wind[wind.y < star.y]
        plus_z_wind = wind[wind.z > star.z]
        minus_z_wind = wind[wind.z < star.z]
        self.assertGreaterEqual(minus_x_wind.vx.min(), -1500 | units.ms)
        self.assertLessEqual(minus_x_wind.vx.max(), -1000 | units.ms)
        self.assertGreaterEqual(plus_x_wind.vx.min(), -1000 | units.ms)
        self.assertLessEqual(plus_x_wind.vx.max(), -500 | units.ms)
        self.assertGreaterEqual(minus_y_wind.vy.min(), -500 | units.ms)
        self.assertLessEqual(minus_y_wind.vy.max(), 0 | units.ms)
        self.assertGreaterEqual(plus_y_wind.vy.min(), 0 | units.ms)
        self.assertLessEqual(plus_y_wind.vy.max(), 500 | units.ms)
        self.assertGreaterEqual(minus_z_wind.vz.min(), 500 | units.ms)
        self.assertLessEqual(minus_z_wind.vz.max(), 1000 | units.ms)
        self.assertGreaterEqual(plus_z_wind.vz.min(), 1000 | units.ms)
        self.assertLessEqual(plus_z_wind.vz.max(), 1500 | units.ms)

    def test_create_wind_short_time(self):
        numpy.random.seed(123457)
        star = self.create_star()
        star_wind = stellar_wind.new_stellar_wind(
            1e-9 | units.MSun)
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(0.01 | units.yr)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 10)
        radii = (wind.position - star.position).lengths()
        self.assertGreaterEqual(radii.min(), star.radius)

    def test_target_gas(self):
        star = self.create_star()
        gas = Particles()
        star_wind = stellar_wind.new_stellar_wind(
            1e-8 | units.MSun, target_gas=gas, timestep=1 | units.day)
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(1 | units.yr)
        self.assertEqual(len(gas), 99)

        star_wind.evolve_model(1.5 | units.yr)
        self.assertEqual(len(gas), 149)

    def test_compensate_gravity(self):
        numpy.random.seed(123457)
        star = self.create_star()
        star_wind = stellar_wind.new_stellar_wind(
            1e-8 | units.MSun, compensate_gravity=True,
            grid_type="random", rotate=False)
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(1 | units.yr)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 100)
        radii = (wind.position - star.position).lengths()
        wind_velocities = (wind.velocity - star.velocity).lengths()
        r_sort = radii.argsort()
        v_sorted = wind_velocities[r_sort]

        print(v_sorted)
        self.assertDecreasing(v_sorted)

        self.assertAlmostEqual(v_sorted[0], 615.945250201 | units.kms)
        self.assertAlmostEqual(v_sorted[-1], 176.050310315 | units.kms)

        return

        from matplotlib import pyplot
        from amuse import plot as aplot
        aplot.plot(radii[r_sort], v_sorted)
        pyplot.show()

    def test_later_star_add(self):
        star1 = self.create_star()
        gas = Particles()
        star_wind = stellar_wind.new_stellar_wind(
            1e-8 | units.MSun, target_gas=gas, timestep=1 | units.day)
        star_wind.particles.add_particles(star1)

        star_wind.evolve_model(1 | units.yr)
        self.assertEqual(len(gas), 99)

        star2 = self.create_star()[0]
        star2.x = 50 | units.parsec
        star_wind.particles.add_particle(star2)

        star_wind.evolve_model(1.5 | units.yr)
        self.assertEqual(len(gas), 198)
        star1_gas = gas[gas.x < 25 | units.parsec]
        star2_gas = gas - star1_gas
        self.assertEqual(len(star1_gas), 149)
        self.assertEqual(len(star2_gas), 49)

    def test_derive_from_evolution(self):
        stars = self.create_stars_without_wind_attributes()
        stars[1].temperature = 16000 | units.K
        stars[1].luminosity = 20 | units.LSun

        star_wind = stellar_wind.new_stellar_wind(
            1e-8 | units.MSun, derive_from_evolution=True, tag_gas_source=True)

        stars = star_wind.particles.add_particles(stars)
        stev = FakeStevCode()
        stev.mass_loss_rate = [1e-7, 2e-7] | units.MSun/units.yr
        stev.particles.add_particles(stars)
        chan = stev.particles.new_channel_to(
            star_wind.particles,
            attributes=["age", "radius", "mass", "luminosity", "temperature"])

        time = 0 | units.yr
        while time < 2.1 | units.yr:
            stev.evolve_model(time)
            chan.copy()
            star_wind.evolve_model(time)
            time += 0.5 | units.yr

        self.assertTrue(star_wind.has_new_wind_particles(),
                        "No new wind particles")
        wind = star_wind.create_wind_particles()
        self.assertEqual(len(wind), 59)

        wind_0 = wind[wind.source == stars[0].key]
        wind_1 = wind[wind.source == stars[1].key]
        self.assertEqual(len(wind_0), 19)
        self.assertEqual(len(wind_1), 40)

        vel_0 = (wind_0.velocity - stars[0].velocity).lengths()
        vel_1 = (wind_1.velocity - stars[1].velocity).lengths()
        self.assertAllAlmostEqual(vel_0, 617.83472287 | units.kms)
        self.assertAllAlmostEqual(vel_1, 864.87956297 | units.kms)

    def test_static_wind_from_evolution(self):
        stars = self.create_stars_without_wind_attributes()
        star_wind = stellar_wind.new_stellar_wind(
            1e-8 | units.MSun, derive_from_evolution=True)
        stev = FakeStevCode()
        stev.mass_loss_rate = [1e-11, 2e-11] | units.MSun/units.yr
        stev.particles.add_particles(stars)

        stellar_wind.static_wind_from_stellar_evolution(
            star_wind, stev, 1.49 | units.Gyr, 1.495 | units.Gyr)

        self.assertAlmostRelativeEquals(
            star_wind.particles.wind_mass_loss_rate,
            [1e-11, 2e-11] | units.MSun/units.yr,
            places=5
            )

        self.assertAlmostRelativeEquals(
            star_wind.particles.terminal_wind_velocity,
            [615.521, 613.198] | units.kms, places=5)


class ConstantSubclass(stellar_wind.AccelerationFunction):
    def acceleration_from_radius(self, radius, star):
        return numpy.zeros_like(radius, dtype=float) | units.m/units.s**2


class RSquaredSubclass(stellar_wind.AccelerationFunction):
    def scaling(self, star):
        return 0.5 * ((star.terminal_wind_velocity**2
                       - star.initial_wind_velocity**2)
                      / (1./star.radius - 1./star.acc_cutoff))

    def acceleration_from_radius(self, r, star):
        return self.scaling(star)/r**2


class TestAccelerationFunctions(TestStellarWind):
    def create_star(self):
        star = Particle()
        star.radius = 2 | units.RSun
        star.acc_cutoff = 10 | units.RSun
        star.initial_wind_velocity = 4 | units.kms
        star.terminal_wind_velocity = 12 | units.kms
        return star

    def test_default_acceleration_1(self):
        self.skip_no_scipy()

        func = ConstantSubclass()
        self.constant_asserts(func)

    def slowtest_default_acceleration_2(self):
        self.skip_no_scipy()

        func = RSquaredSubclass()
        self.r_squared_asserts(func)

    def test_constant_velocity(self):
        func = stellar_wind.ConstantVelocityAcceleration()
        self.constant_asserts(func)

    def test_r_squared_velocity(self):
        self.skip_no_scipy()

        func = stellar_wind.RSquaredAcceleration()
        self.r_squared_asserts(func)

    def constant_asserts(self, func):
        star = self.create_star()
        star.initial_wind_velocity = star.terminal_wind_velocity
        radii = [2., 4.] | units.RSun

        accelerations = func.acceleration_from_radius(radii, star)
        self.assertAlmostEqual(accelerations, [0., 0.] | units.m/units.s**2)

        velocities = func.velocity_from_radius(radii, star)
        self.assertAlmostEqual(velocities, [12., 12.] | units.kms)

        times = (numpy.array([0, 1, 2, 3]) * star.radius
                 / star.terminal_wind_velocity)
        new_radii = func.radius_from_time(times, star)
        self.assertAlmostEqual(new_radii, [2., 4., 6., 8.] | units.RSun)

        random_numbers = numpy.linspace(0., 1., 4)

        radii = func.radius_from_number(random_numbers, 5 | units.RSun, star)
        self.assertAlmostEqual(radii, [2, 3., 4., 5] | units.RSun)

    def r_squared_asserts(self, func):
        star = self.create_star()
        radii = [2., 4., 8., 10.] | units.RSun

        accelerations = func.acceleration_from_radius(radii, star)
        expected = (numpy.array([40., 10., 2.5, 1.6]) * 0.0014378145
                    | units.m/units.s**2)
        self.assertAlmostEqual(accelerations, expected)

        velocities = func.velocity_from_radius(radii, star)
        expected = [4., 9.79795897113, 11.6619037897, 12.] | units.kms
        self.assertAlmostEqual(velocities, expected)

        times = (numpy.array([0., 1., 2., 3.]) * star.radius
                 / star.initial_wind_velocity)
        new_radii = func.radius_from_time(times, star)
        exp = [2.0, 6.50863038599, 12.4140692531, 18.4140691706] | units.RSun
        self.assertAlmostEqual(new_radii, exp)

        x = numpy.linspace(0., 1., 4)
        new_radii = func.radius_from_number(x, 5 | units.RSun, star)
        expected = [2,  2.72665604616,  3.7783430146, 5] | units.RSun
        self.assertAlmostEqual(new_radii, expected)

    def test_radius_from_time(self):
        self.skip_no_scipy()

        func = stellar_wind.RSquaredAcceleration()
        star = Particle()
        star.radius = 312. | units.RSun
        star.acc_cutoff = 312*5. | units.RSun
        star.initial_wind_velocity = 2 | units.kms
        star.terminal_wind_velocity = 20 | units.kms

        times = [0., 25., 100.] | units.yr
        print("times", times)
        radii = func.radius_from_time(times, star)
        print("radii", radii)

        self.assertAlmostEqual(radii[0], 312. | units.RSun)
        self.assertAlmostEqual(radii[1], 22644.6086263 | units.RSun)
        self.assertAlmostEqual(radii[2], 90704.1183516 | units.RSun)

        return

        times = numpy.linspace(0., 1, 100) | units.yr
        radii = func.radius_from_time(times, star)

        print(radii)
        from matplotlib import pyplot
        from amuse import plot as aplot
        aplot.plot(times, radii/star.radius)
        pyplot.show()

    def test_radius_from_number(self):
        self.skip_no_scipy()

        star = Particle(
            initial_wind_velocity=200. | units.kms, mass=20. | units.MSun,
            radius=10. | units.RSun, terminal_wind_velocity=1000. | units.kms)
        star.acc_cutoff = 5. * star.radius
        r_max = 24.2642526438 | units.RSun
        x = 0.1

        func = stellar_wind.RSquaredAcceleration()
        radius = func.radius_from_number(x, r_max, star)
        self.assertAlmostEqual(radius,  10.6029030808 | units.RSun)


class TestAcceleratingWind(TestStellarWind):
    def test_wind_creation_constant(self):
        numpy.random.seed(123457)
        star = self.create_star()
        star.initial_wind_velocity = star.terminal_wind_velocity
        star_wind = stellar_wind.new_stellar_wind(
            3e-11 | units.MSun,
            mode="accelerate",
            acceleration_function="constant_velocity",
            )
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(10 | units.day)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 912)
        radii = (wind.position - star.position).lengths()
        velocities = (wind.velocity - star.velocity).lengths()

        self.assertGreaterEqual(radii.min(), star.radius)
        self.assertLessEqual(radii.max(), 8.212 | units.RSun)
        self.assertAlmostEqual(radii.mean(), 2.32897955194 | units.RSun)

        self.assertAlmostEqual(velocities.min(), star.initial_wind_velocity)
        self.assertAlmostEqual(velocities.max(), star.initial_wind_velocity)

        return

        from matplotlib import pyplot
        from amuse import plot as aplot
        import pylab
        n, bins, patches = pylab.hist(radii.value_in(units.RSun), 10)
        pylab.show()
        pyplot.clf()
        aplot.scatter(radii, velocities)
        pyplot.show()

    def test_wind_creation_rsquared(self):
        self.skip_no_scipy()

        numpy.random.seed(123457)
        star = self.create_star()
        star_wind = stellar_wind.new_stellar_wind(
            3e-9 | units.MSun,
            mode="accelerate",
            acceleration_function="rsquared",
            )
        star_wind.particles.add_particles(star)

        star_wind.evolve_model(10 | units.day)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 9)
        radii = (wind.position - star.position).lengths()
        velocities = (wind.velocity - star.velocity).lengths()

        self.assertGreaterEqual(radii.min(), star.radius)
        self.assertLessEqual(radii.max(), 12.023 | units.RSun)
        self.assertAlmostEqual(radii.mean(), 2.06872960025 | units.RSun)

        self.assertGreaterEqual(velocities.min(), star.initial_wind_velocity)
        self.assertLessEqual(velocities.max(), star.terminal_wind_velocity)
        self.assertAlmostEqual(velocities.mean(), 0.103550548126 | units.kms)

    def test_acceleration(self):
        star = self.create_star()
        star.position = [1, 1, 1] | units.RSun
        star_wind = stellar_wind.new_stellar_wind(
            3e-11 | units.MSun,
            mode="accelerate",
            acceleration_function="rsquared",
            compensate_gravity=False,
            r_out_ratio=5,
            )
        star = star_wind.particles.add_particles(star)

        distances = [2, 4, 5, 10] | units.RSun
        accelerations = star_wind.acceleration(star, distances)
        expected = ([0.00011120596693, 2.78014917326e-5, 1.77929547088e-5, 0.]
                    | units.m/units.s**2)

        self.assertAlmostEqual(accelerations, expected)

    def test_gravity_compensation(self):
        star = self.create_star()
        star.position = [1, 1, 1] | units.RSun
        star_wind = stellar_wind.new_stellar_wind(
            3e-11 | units.MSun,
            mode="accelerate",
            acceleration_function="constant_velocity",
            compensate_gravity=True,
            )
        star = star_wind.particles.add_particles(star)

        distances = [2, 4, 5, 10, 20] | units.RSun
        accelerations = star_wind.acceleration(star, distances)
        expected = ([137.213699216, 34.3034248039, 21.9541918745,
                     5.48854796862, 1.37213699216] | units.m/units.s**2)

        self.assertAlmostEqual(accelerations, expected)

    def test_pressure_compensation(self):
        star = self.create_star()
        star.position = [1, 1, 1] | units.RSun
        star_wind = stellar_wind.new_stellar_wind(
            3e-11 | units.MSun,
            mode="accelerate",
            acceleration_function="constant_velocity",
            compensate_gravity=False,
            compensate_pressure=True
            )
        star = star_wind.particles.add_particles(star)

        distances = [2, 4, 5, 10] | units.RSun
        accelerations = star_wind.acceleration(star, distances)
        print(accelerations)
        expected = ([-0.0187296577097, -0.00371643479392, -0.00220802076676,
                     -0.000438126811] | units.m/units.s**2)

        self.assertAlmostEqual(accelerations, expected)

    def test_get_gravity_at_point(self):
        star = self.create_star()
        star.position = [1, 1, 1] | units.RSun
        star_wind = stellar_wind.new_stellar_wind(
            3e-11 | units.MSun,
            mode="accelerate",
            acceleration_function="rsquared",
            compensate_gravity=False,
            compensate_pressure=False,
            r_out_ratio=5,
            )
        star_wind.particles.add_particles(star)

        # unaffected points
        points = [[1, 1, -10], [1, 12, 0], [1, -10, 1],
                  [12, 0, 1]] | units.RSun
        x, y, z = points.transpose()
        ax, ay, az = star_wind.get_gravity_at_point(1, x, y, z)
        zeros = [0, 0, 0, 0] | units.m/units.s**2
        self.assertEqual(ax, zeros)
        self.assertEqual(ay, zeros)
        self.assertEqual(az, zeros)

        # affected points
        points = [[5.1, 1, 1], [1, 1, 5.1], [1, 7, 1], [1, 7, 8],
                  [-8, 1, 1]] | units.RSun
        x, y, z = points.transpose()
        ax, ay, az = star_wind.get_gravity_at_point(1, x, y, z)

        result = quantities.as_vector_quantity([ax, ay, az]).transpose()
        expected = [[2.64618600667e-05, 0.0, 0.0],
                    [0.0, 0.0, 2.64618600667e-05],
                    [0.0, 1.23562185478e-05, 0.0],
                    [0.0, 3.40573571553e-06, 3.97335833479e-06],
                    [-5.49165268791e-06, 0.0, 0.0]] | units.m / units.s**2
        self.assertAlmostEqual(result, expected)


class TestHeatingWind(TestStellarWind):
    def test_wind_creation(self):
        numpy.random.seed(123456789)
        star = self.create_star()
        feedback_efficiency = 0.01
        r_max = 10 | units.RSun
        star_wind = stellar_wind.new_stellar_wind(
            3e-10 | units.MSun,
            mode="heating",
            r_max=r_max,
            feedback_efficiency=feedback_efficiency,
            )
        star_wind.particles.add_particles(star)

        evolve_time = 1 | units.day
        expected_u = 0.5*feedback_efficiency*star.terminal_wind_velocity**2

        star_wind.evolve_model(evolve_time)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 9)
        max_dist = (wind.position - star.position).lengths().max()
        self.assertLessEqual(max_dist, r_max + star.radius)
        self.assertAlmostEqual(max_dist, 8.69397162579 | units.RSun)

        self.assertAllAlmostEqual(wind.u, expected_u)

        for t in numpy.arange(2, 10) * evolve_time:
            star_wind.evolve_model(t)
        wind = star_wind.create_wind_particles()

        self.assertEqual(len(wind), 73)
        max_dist = (wind.position - star.position).lengths().max()
        self.assertLessEqual(max_dist, r_max + star.radius)

        self.assertAllAlmostEqual(wind.u, expected_u)

    def setup_supernova(self):
        numpy.random.seed(123456789)

        stars = Particles(2)
        stars.mass = [9, 10] | units.MSun
        stars[0].position = [1, 1, 1] | units.parsec
        stars[1].position = [-1, -1, -1] | units.parsec
        stars[0].velocity = [-1000, 0, 1000] | units.ms
        stars[1].velocity = [0, 0, 0] | units.ms

        stev = SeBa()
        stars = stev.particles.add_particles(stars)

        r_max = .1 | units.parsec
        star_wind = stellar_wind.new_stellar_wind(
            3e-5 | units.MSun,
            mode="heating",
            r_max=r_max,
            derive_from_evolution=True,
            tag_gas_source=True
            )
        star_wind.particles.add_particles(stars)

        return stev, star_wind, stars

    def test_supernova(self):
        stev, star_wind, stars = self.setup_supernova()
        chan = stev.particles.new_channel_to(
            star_wind.particles,
            attributes=["age", "radius", "mass", "luminosity", "temperature", "stellar_type"])


        dt = 5 | units.Myr
        t = 0 | units.Myr
        t_end = 31 | units.Myr

        wind_N = []
        while t < t_end:
            stev.evolve_model(t)
            chan.copy()
            star_wind.evolve_model(t)
            if star_wind.has_new_wind_particles():
                wind = star_wind.create_wind_particles()
                wind_1 = wind[wind.source == stars[0].key]
                wind_2 = wind - wind_1
                wind_N.append([len(wind_1), len(wind_2)])
                if len(wind_2) > 0:
                    print("time", t, "wind energy", (wind_2.u * wind_2.mass).sum())
            else:
                wind_N.append([0, 0])

            t += dt


        self.assertEqual(wind_N, [[0, 0], [32, 45], [57, 59], [114, 130], [302, 635], [1231, 5810], [2981, 285777]])
        supernova = wind_2
        sn_energy = (supernova.u * supernova.mass).sum()
        self.assertAlmostRelativeEqual(sn_energy, 1e49 | units.erg, 2)


    def test_supernova_manual(self):
        stev, star_wind, stars = self.setup_supernova()
        stev.stopping_conditions.supernova_detection.enable()

        chan = stev.particles.new_channel_to(
            star_wind.particles,
            attributes=["age", "radius", "mass", "luminosity", "temperature"])


        dt = 5 | units.Myr
        t = 0 | units.Myr
        t_end = 31 | units.Myr

        wind_N = []
        wind_E = [] | units.erg
        while stev.model_time < t_end:
            stev.evolve_model(t)
            chan.copy()

            if stev.stopping_conditions.supernova_detection.is_set():
                key = stev.stopping_conditions.supernova_detection.particles(0)[0].key
                star = star_wind.particles._get_particle(key)
                star.mass_loss_type = "supernova"
                print("supernova detected at time", t, stev.model_time)

            star_wind.evolve_model(t)

            if star_wind.has_new_wind_particles():
                wind = star_wind.create_wind_particles()
                wind_2 = wind[wind.source == stars[1].key]
                wind_N.append(len(wind_2))
                wind_E.append((wind_2.u * wind_2.mass).sum())
                print(len(wind), len(wind_2))

            if not stev.stopping_conditions.supernova_detection.is_set():
                t += dt

        print(wind_N)
        print(wind_E)
        self.assertEqual(wind_N[-2:], [606194, 0])
        self.assertAlmostRelativeEqual(wind_E[-2:], [1.e49, 0] | units.erg,7)
