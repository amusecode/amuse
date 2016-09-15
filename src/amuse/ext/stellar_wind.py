import numpy

from amuse.support.exceptions import AmuseException
from amuse.datamodel import Particles
from amuse.units import units, quantities, constants


def kudritzki_wind_velocity(mass, radius, luminosity, temperature,
                            Y=0.25, I_He=2):
    """
      This routine calculates the escape and terminal wind velocity. The
      Equations are taken from Kudritzki & Puls, Annual Reviews of Astronomy
      and Astrophysics, 2000, Vol. 38, p.613-666 Equation (8) and (9) and
      Kudritzki et al., 1989, A&A 219, 205 Equation (64) and (65).

      I_He:    Number of electrons per He nucleus (= 2 in O-Stars)
      sigma_e:   Thomson absorption coefficient
      Gamma:     Ratio of radiative Thomson to gravitational acceleration
    """
    sigma_e = 0.398 * (1 + I_He*Y)/(1 + 4*Y)
    Gamma = 7.66E-5 * sigma_e * (luminosity.value_in(units.LSun)
                                 / mass.value_in(units.MSun))
    v_esc = (2*constants.G * mass / radius*(1 - Gamma))**0.5

    condlist = [temperature >= 21000. | units.K,
                (10000. | units.K < temperature) &
                (temperature < 21000. | units.K),
                temperature <= 10000. | units.K]
    choicelist = [2.65, 1.4, 1.0]

    return v_esc * numpy.select(condlist, choicelist)


class PositionGenerator(object):
    def __init__(self, grid_type="regular", rotate=True):
        self.cube_generator = {
            "random": self.random_cube,
            "regular": self.regular_grid_unit_cube,
            }[grid_type]

        self.rotate = rotate

    def as_three_vector(self, array):
        number = array
        if quantities.is_quantity(array):
            number = array.number
        three_vector = numpy.transpose([number]*3)
        if quantities.is_quantity(array):
            three_vector = three_vector | array.unit
        return three_vector

    def regular_grid_unit_cube(self, N):
        n = int(numpy.ceil(N**(1./3.)))
        start = -1. + 1./n
        stop = 1. - 1./n
        # complex step number tells mgrid to work like linspace
        step = n*1j

        grid = numpy.mgrid[start: stop: step,
                           start: stop: step,
                           start: stop: step]

        grid = grid.reshape(3, n**3)
        grid = grid.transpose()
        grid = grid[numpy.random.choice(n**3, size=N, replace=False)]
        return grid

    def random_cube(self, N):
        numbers = numpy.random.uniform(-1., 1., 3 * N)
        return numpy.reshape(numbers, (N, 3))

    def cutout_sphere(self, positions, rmin):
        r = numpy.sqrt((positions**2).sum(1))
        return positions[(r >= rmin) & (r < 1)]

    def random_rotation(self):
        u, v, r = numpy.random.rand(3)
        theta = 2. * numpy.pi * u
        phi = numpy.arccos(2. * v - 1.)
        axis = (numpy.sin(theta) * numpy.cos(phi),
                numpy.sin(theta) * numpy.sin(phi),
                numpy.cos(theta))
        axis = numpy.asarray(axis)
        angle = 2. * numpy.pi * r
        return axis, angle

    def rotation_matrix(self, axis, angle):
        """ Using the Euler-Rodrigues formula """
        axis = numpy.asarray(axis)
        theta = numpy.asarray(angle)
        axis = axis/numpy.sqrt(numpy.dot(axis, axis))
        a = numpy.cos(theta/2.)
        b, c, d = -axis * numpy.sin(theta/2.)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        m = [[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
             [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
             [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]]
        return numpy.asarray(m)

    def rotate_positions(self, positions, axis, angle):
        matrix = self.rotation_matrix(axis, angle)
        rotated = [matrix.dot(p) for p in positions]
        return numpy.asarray(rotated)

    def uniform_hollow_sphere(self, N, rmin):
        cube_sphere_ratio = 4/3. * numpy.pi * 0.5**3 * (1 - rmin**3)
        estimatedN = N / cube_sphere_ratio

        while True:
            estimatedN = estimatedN * 1.1 + 1
            cube = self.cube_generator(int(estimatedN))
            hollow_sphere = self.cutout_sphere(cube, rmin)
            if len(hollow_sphere) >= N:
                break

        if self.rotate:
            axis, angle = self.random_rotation()
            hollow_sphere = self.rotate_positions(hollow_sphere, axis, angle)

        return hollow_sphere[:N]

    def generate_positions(self, N, rmin, rmax, radius_function=None,
                           star=None):
        """
            The particles start out in a (random) position between
            the surface of the star and the distance that the
            previously released particles have reached.
            This assumes that the current wind velocity is
            comparable to the previous wind velocity.

            Note that the stellar position is not added yet here.
        """
        positions = self.uniform_hollow_sphere(N, 1. * rmin / rmax)
        vector_lengths = numpy.sqrt((positions**2).sum(1))

        unit_vectors = positions/self.as_three_vector(vector_lengths)

        int_v_over_total = (((vector_lengths * rmax)**3 - rmin**3)
                            / (rmax**3 - rmin**3))

        if radius_function is not None:
            distance = radius_function(int_v_over_total, rmax, star)
        else:
            distance = int_v_over_total * (rmax - rmin) + rmin

        position = unit_vectors * self.as_three_vector(distance)

        return position, unit_vectors


class StarsWithMassLoss(Particles):
    def __init__(self, *args, **kwargs):
        super(StarsWithMassLoss, self).__init__(*args, **kwargs)
        self._private.timestamp = 0. | units.yr
        self._private.previous_time = 0. | units.yr
        self._private.track_mechanical_energy = False
        self._private.new_unset_lmech_particles = False

        self._private.attribute_names = set(["lost_mass",
                                             "wind_release_time",
                                             "mu",
                                             "mass",
                                             "radius",
                                             "age",
                                             "temperature",
                                             "luminosity",
                                             "stellar_type",
                                             "x",
                                             "y",
                                             "z",
                                             "vx",
                                             "vy",
                                             "vz",
                                             "wind_mass_loss_rate",
                                             "initial_wind_velocity",
                                             "terminal_wind_velocity",
                                             "mass_loss_type",
                                             ])
        self._private.defaults = dict(lost_mass=0 | units.MSun,
                                      mass=0 | units.MSun,
                                      radius=0 | units.RSun,
                                      age=0 | units.Myr,
                                      temperature=0 | units.K ,
                                      luminosity = 0 | units.LSun,
                                      stellar_type = 1 | units.stellar_type,
                                      x=0 | units.m,
                                      y=0 | units.m,
                                      z=0 | units.m,
                                      vx=0 | units.ms,
                                      vy=0 | units.ms,
                                      vz=0 | units.ms,
                                      wind_mass_loss_rate=0 | units.MSun/units.yr,
                                      initial_wind_velocity=0 | units.ms,
                                      terminal_wind_velocity=0 | units.ms,
                                      mechanical_energy=0 | units.J,
                                      mass_loss_type="wind",
                                      )

        self.set_global_mu()


    def add_particles(self, particles, *args, **kwargs):
        new_particles = super(StarsWithMassLoss, self).add_particles(
            particles, *args, **kwargs)

        return new_particles

    def can_extend_attributes(self):
        return False

    def get_attribute_names_defined_in_store(self):
        return list(self._private.attribute_names) if len(self) > 0 else []

    def add_particles_to_store(self, keys, attributes=[], values=[]):
        good_attributes = self._private.attribute_names

        # add default values if missing
        for attr in good_attributes:
            if attr not in attributes:
                attributes.append(attr)

                if attr == "wind_release_time":
                    value = self.collection_attributes.timestamp or 0 | units.yr
                elif attr == "previous_age":
                    if "age" in attributes:
                        value = values[attributes.index("age")]
                    else:
                        value = 0 | units.yr
                elif attr == "previous_mass":
                    if "mass" in attributes:
                        value = values[attributes.index("mass")]
                    else:
                        value = 0 | units.MSun
                elif attr == 'previous_mechanical_luminosity':
                    value = -1 | units.W
                    self._private.new_unset_lmech_particles = True
                else:
                    value = self._private.defaults[attr]

                values.append(value)

        # remove unsupported attributes
        if len(attributes) > len(good_attributes):
            values = [v for i, v in enumerate(values)
                      if attributes[i] in good_attributes]
            attributes = [a for a in attributes if a in good_attributes]

        super(StarsWithMassLoss, self).add_particles_to_store(keys, attributes, values)

    def set_values_in_store(self, indices, attributes, list_of_values_to_set):
        for attr in attributes:
            if attr not in self._private.attribute_names:
                raise AttributeError("You tried to set attribute '{0}'"
                    " but this attribute is not accepted for this set."
                    .format(attr))

        # TODO
        super(StarsWithMassLoss, self).set_values_in_store(
            indices, attributes, list_of_values_to_set)

    def add_calculated_attribute(self, name_of_the_attribute, *args, **kwargs):
        if name_of_the_attribute in self._private.attribute_names:
            self._private.attribute_names.remove(name_of_the_attribute)
            del self._private.defaults[name_of_the_attribute]

        super(StarsWithMassLoss, self).add_calculated_attribute(
            name_of_the_attribute, *args, **kwargs)

    def evolve_mass_loss(self, time):
        if self._private.previous_time >= time:
            return

        elapsed_time = time - self._private.previous_time
        self.lost_mass += elapsed_time * self.wind_mass_loss_rate

        if self._private.track_mechanical_energy:
            new_mechanical_luminosity = (0.5 * self.wind_mass_loss_rate
                                         * self.terminal_wind_velocity**2)

            if self._private.new_unset_lmech_particles:
                i_new = self.previous_mechanical_luminosity < quantities.zero
                self[i_new].previous_mechanical_luminosity =\
                    new_mechanical_luminosity[i_new]
                self._private.new_unset_lmech_particles = False

            average_mechanical_luminosity = 0.5 * (
                self.previous_mechanical_luminosity
                + new_mechanical_luminosity)
            self.mechanical_energy += (elapsed_time
                                       * average_mechanical_luminosity)

            self.previous_mechanical_luminosity = new_mechanical_luminosity

        self.collection_attributes.timestamp = time
        self._private.previous_time = time

    def track_mechanical_energy(self, track=True):
        self._private.track_mechanical_energy = track
        mech_attrs = set(["mechanical_energy", "previous_mechanical_luminosity"])
        if track:
            self._private.attribute_names |= mech_attrs
        else:
            self._private.attribute_names -= mech_attrs

    def set_global_mu(self, mu=None, Y=0.25, Z=0.02, x_ion=0.1):
        """
            Set the global value of mu used to create stellar wind.
            if mu is added directly, Y (Helium fraction), Z (metal fraction)
            and x_ion (percentage of ionized atoms) are ignored.
            An alternative way is to set mu for each star separately.
        """
        if mu is None:
            X = 1.0 - Y - Z
            ion_num = X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0
            mu = constants.proton_mass / ion_num

        self.mu = mu
        self._private.defaults['mu'] = mu

    def reset(self):
        self.lost_mass = 0.0 | units.MSun
        self.set_begin_time(0. | units.yr)

    def set_begin_time(self, time):
        self.wind_release_time = time
        self.collection_attributes.timestamp = time
        self._private.previous_time = time


class EvolvingStarsWithMassLoss(StarsWithMassLoss):
    """
        Derive the stellar wind from stellar evolution.
        You have to copy the relevant attributes from the stellar evolution.
        This can be done using a channel like:

        chan = stellar_evolution.particles.new_channel_to(
            stellar_wind.particles,
            attributes=["age", "radius", "mass", "luminosity", "temperature"])

        while <every timestep>:
            chan.copy()
    """

    def __init__(self, *args, **kwargs):
        super(EvolvingStarsWithMassLoss, self).__init__(*args, **kwargs)
        attrs = set(["previous_age",
                     "previous_mass",
                     ])
        self._private.attribute_names |= attrs

    def add_particles(self, particles, *args, **kwargs):
        new_particles = super(EvolvingStarsWithMassLoss, self).add_particles(
            particles, *args, **kwargs)

        return new_particles

    def evolve_mass_loss(self, time):
        if self._private.previous_time <= time:
            self.update_from_evolution()
            StarsWithMassLoss.evolve_mass_loss(self, time)

    def update_from_evolution(self):
        if (self.age != self.previous_age).any():
            mass_loss = self.previous_mass - self.mass
            timestep = self.age - self.previous_age
            self.wind_mass_loss_rate = mass_loss / timestep

            self.previous_age = self.age
            self.previous_mass = self.mass


class SimpleWind(PositionGenerator):
    """
        The simple wind model creates SPH particles moving away
        from the star at the terminal velocity.
        This is a safe assumption if the distance to other objects
        is (far) larger then the stellar radius.
    """

    def __init__(self, sph_particle_mass, derive_from_evolution=False,
                 tag_gas_source=False, compensate_gravity=False, **kwargs):
        self.r_max = kwargs.pop("r_max", None)
        super(SimpleWind, self).__init__(**kwargs)
        self.sph_particle_mass = sph_particle_mass
        self.model_time = 0.0 | units.yr

        if derive_from_evolution:
            self.particles = EvolvingStarsWithMassLoss()
            self.particles.add_calculated_attribute(
                "terminal_wind_velocity", kudritzki_wind_velocity,
                attributes_names=['mass', 'radius',
                                  'luminosity', 'temperature'])
        else:
            self.particles = StarsWithMassLoss()

        self.target_gas = self.timestep = None
        self.tag_gas_source = tag_gas_source
        self.compensate_gravity = compensate_gravity

        self.internal_energy_formula = self.internal_energy_from_temperature
        self.set_initial_wind_velocity()

    def set_initial_wind_velocity(self):
        self.particles.add_calculated_attribute(
            "initial_wind_velocity", lambda v: v,
            attributes_names=['terminal_wind_velocity'])

    def evolve_particles(self):
        self.particles.evolve_mass_loss(self.model_time)

    def evolve_model(self, time):
        if self.has_target():
            while self.model_time <= time:
                self.evolve_particles()
                if self.has_new_wind_particles():
                    wind_gas = self.create_wind_particles()
                    self.target_gas.add_particles(wind_gas)
                self.model_time += self.timestep
        else:
            self.model_time = time
            self.evolve_particles()

    def set_target_gas(self, target_gas, timestep):
        self.target_gas = target_gas
        self.timestep = timestep

    def has_target(self):
        return self.target_gas is not None

    def internal_energy_from_temperature(self, star, wind=None):
        """
            set the internal energy from the stellar surface temperature.
        """
        return (3./2. * constants.kB * star.temperature / star.mu)

    def internal_energy_from_velocity(self, star, wind=None):
        """
            set the internal energy from the terminal wind velocity.
        """

        return 0.5 * star.terminal_wind_velocity**2

    def wind_sphere(self, star, Ngas):
        wind = Particles(Ngas)

        wind_velocity = star.initial_wind_velocity

        outer_wind_distance = star.radius + wind_velocity * (
            self.model_time - star.wind_release_time)

        if self.r_max is not None and outer_wind_distance < self.r_max:
            outer_wind_distance = self.r_max

        wind.position, direction = self.generate_positions(
            Ngas, star.radius, outer_wind_distance)

        if self.compensate_gravity:
            r = wind.position.lengths()
            escape_velocity_squared = 2. * constants.G * star.mass / r
            speed = (wind_velocity**2 + escape_velocity_squared).sqrt()
            wind.velocity = self.as_three_vector(speed) * direction
        else:
            wind.velocity = direction * wind_velocity

        return wind

    def create_wind_particles_for_one_star(self, star):
        Ngas = int(star.lost_mass/self.sph_particle_mass)
        star.lost_mass -= Ngas * self.sph_particle_mass

        wind = self.wind_sphere(star, Ngas)

        wind.mass = self.sph_particle_mass
        wind.u = self.internal_energy_formula(star, wind)
        wind.position += star.position
        wind.velocity += star.velocity

        if self.tag_gas_source:
            wind.source = star.key

        return wind

    def create_wind_particles(self):
        wind = Particles(0)

        for star in self.particles:
            if star.lost_mass > self.sph_particle_mass:
                new_particles = self.create_wind_particles_for_one_star(star)
                wind.add_particles(new_particles)
                star.wind_release_time = self.model_time

        return wind

    def has_new_wind_particles(self):
        return self.particles.lost_mass.max() > self.sph_particle_mass

    def create_initial_wind_for_time(self, time, check_length=True):
        """
            Particles are created as if the wind has already been blowing for
            'time'.  Note that this does not work if the mass loss is derived
            from stellar evolution.
        """
        self.model_time = time
        self.particles.evolve_mass_loss(self.model_time)

        if self.has_new_wind_particles():
            wind_gas = self.create_wind_particles()
            if self.has_target():
                self.target_gas.add_particles(wind_gas)
        elif check_length:
            raise AmuseException("create_initial_wind time was too small to"
                                 "create any particles.")
        else:
            wind_gas = Particles()

        self.reset()

        return wind_gas

    def create_initial_wind(self, number):
        """
            This is a convenience method that creates some initial particles.
        """
        required_mass = number * self.sph_particle_mass
        total_mass_loss = self.particles.wind_mass_loss_rate.sum()
        time = 1.0 * required_mass/total_mass_loss
        wind_gas = Particles()
        while len(wind_gas) < number:
            time = 1.1 * time
            wind_gas = self.create_initial_wind_for_time(time, False)

        return wind_gas[:number]

    def reset(self):
        self.particles.reset()
        self.model_time = 0.0 | units.yr

    def set_begin_time(self, time):
        self.model_time = time
        self.particles.set_begin_time(time)

    def get_gravity_at_point(self, eps, x, y, z):
        return [0, 0, 0] | units.m/units.s**2

    def get_potential_at_point(self, radius, x, y, z):
        return [0, 0, 0] | units.J


class AccelerationFunction(object):
    """
    Abstact superclass of all acceleration functions.
    It numerically derives everything using acceleration_from_radius
    Overwrite as many of these functions with analitic solutions as possible.
    """

    def __init__(self):
        try:
            from scipy import integrate, optimize
            self.quad = integrate.quad
            self.brentq = optimize.brentq
        except ImportError:
            self.quad = self.unsupported
            self.brentq = self.unsupported

    def unsupported(self, *args, **kwargs):
        raise AmuseException("Importing SciPy has failed")

    def acceleration_from_radius(self, radius, star):
        """
            to be overridden
        """
        pass

    def velocity_from_radius(self, radius, star):
        def stripped_acceleration(r1):
            acc = self.acceleration_from_radius(r1 | units.RSun, star)
            return acc.value_in(units.RSun/units.yr**2)

        def acc_integral(r):
            start = star.radius.value_in(units.RSun)

            result = self.quad(stripped_acceleration, start, r)
            return result[0]

        integral = numpy.vectorize(acc_integral)(radius.value_in(units.RSun))
        integral = integral | units.RSun**2/units.yr**2

        return (2. * integral + star.initial_wind_velocity**2).sqrt()

    def radius_from_time(self, time, star):
        """
            following http://math.stackexchange.com/questions/54586/
            converting-a-function-for-velocity-vs-position-vx-to-position-vs-time
        """

        def inverse_velocity(r1):
            velocity = self.velocity_from_radius(r1 | units.RSun, star)
            return 1. / velocity.value_in(units.RSun/units.yr)

        def one_radius(t):
            def residual(r2):
                start = star.radius.value_in(units.RSun)
                result = self.quad(inverse_velocity, start, r2)
                return result[0] - t.value_in(units.yr)

            start = star.radius.value_in(units.RSun)
            end = 1e5 * start
            result = self.brentq(residual, start, end)

            return result

        radius = numpy.vectorize(one_radius)(time)
        return radius | units.RSun

    def radius_from_number(self, numbers, max_radius, star):
        """
            See http://www.av8n.com/physics/arbitrary-probability.htm
            for some good info on this.
        """

        rmin = star.radius.value_in(units.RSun)
        rmax = max_radius.value_in(units.RSun)

        def inverse_velocity(r1):
            velocity = self.velocity_from_radius(r1 | units.RSun, star)
            velocity = velocity.value_in(units.RSun/units.s)
            return 1. / velocity

        def cumulative_inverse_velocity(q):
            res = self.quad(inverse_velocity, rmin, q)
            return res[0]

        d_max = cumulative_inverse_velocity(rmax)

        def one_radius(x):
            def residual(r2):
                return cumulative_inverse_velocity(r2) / d_max - x

            return self.brentq(residual, rmin, rmax)

        radius = numpy.vectorize(one_radius)(numbers)

        return radius | units.RSun

    def fix_cutoffs(self, test, value, star, default):
        if hasattr(value, "__len__"):
            value[test] = default
        elif test:
            value = default
        return value

    def fix_acc_cutoff(self, r, acc, star):
        if star.acc_cutoff is None:
            return acc

        test = r > star.acc_cutoff
        return self.fix_cutoffs(test, acc, star, quantities.zero)

    def fix_v_cutoff(self, r, v, star):
        if star.acc_cutoff is None:
            return v
        test = r > star.acc_cutoff
        return self.fix_cutoffs(test, v, star, star.terminal_wind_velocity)


class ConstantVelocityAcceleration(AccelerationFunction):
    """
        A very basic "acceleration" function that ensures a constant velocity,
    """
    def __init__(self, use_initial=False):
        super(ConstantVelocityAcceleration, self).__init__()
        self.use_initial = use_initial

    def velocity(self, star):
        if self.use_initial:
            return star.initial_wind_velocity
        else:
            return star.terminal_wind_velocity

    def acceleration_from_radius(self, radius, star):
        return numpy.zeros_like(radius, dtype=float) | units.m/units.s**2

    def velocity_from_radius(self, radius, star):
        return numpy.ones_like(radius, dtype=float) * self.velocity(star)

    def radius_from_time(self, t, star):
        return star.radius + t * self.velocity(star)

    def radius_from_number(self, x, r_max, star):
        r_star = star.radius
        return x * (r_max - r_star) + r_star


class RSquaredAcceleration(AccelerationFunction):
    def scaling(self, star):
        denominator = 1./star.radius

        if star.acc_cutoff is not None:
            denominator = denominator - 1./star.acc_cutoff

        numerator = star.terminal_wind_velocity**2 - star.initial_wind_velocity**2
        return 0.5 * numerator / denominator

    def acceleration_from_radius(self, r, star):
        acc = self.scaling(star)/r**2
        return self.fix_acc_cutoff(r, acc, star)

    def velocity_from_radius(self, r, star):
        v = (2 * self.scaling(star) * (1./star.radius - 1./r)
             + star.initial_wind_velocity**2).sqrt()
        return self.fix_v_cutoff(r, v, star)


class DelayedRSquaredAcceleration(AccelerationFunction):
    def scaling(self, star):
        denominator = 1./star.acc_start

        if star.acc_cutoff is not None:
            denominator = denominator - 1./star.acc_cutoff

        numerator = star.terminal_wind_velocity**2 - star.initial_wind_velocity**2
        return 0.5 * numerator / denominator

    def fix_acc_start_cutoff(self, r, acc, star):
        return self.fix_cutoffs(r < star.acc_start, acc, star, quantities.zero)

    def fix_v_start_cutoff(self, r, v, star):
        test = r < star.acc_start
        return self.fix_cutoffs(test, v, star, star.initial_wind_velocity)

    def acceleration_from_radius(self, r, star):
        acc = self.scaling(star)/r**2
        acc = self.fix_acc_start_cutoff(r, acc, star)
        return self.fix_acc_cutoff(r, acc, star)

    def velocity_from_radius(self, r, star):
        v = (2 * self.scaling(star) * (1./star.acc_start - 1./r)
             + star.initial_wind_velocity**2).sqrt()
        v = self.fix_v_start_cutoff(r, v, star)
        return self.fix_v_cutoff(r, v, star)


class BetaLawAcceleration(AccelerationFunction):
    """ Following Lamers 1999 and Maciel 2005 """

    def __init__(self, beta=.8):
        super(BetaLawAcceleration, self).__init__()
        self.beta = beta

    def acceleration_from_radius(self, r, star):
        v_diff = star.terminal_wind_velocity - star.initial_wind_velocity
        dvdr = (v_diff * star.radius / r**2
                * self.beta * (1 - star.radius / r)**(self.beta - 1))
        return dvdr * self.velocity_from_radius(r, star)

    def velocity_from_radius(self, r, star):
        v_start = star.initial_wind_velocity
        v_end = star.terminal_wind_velocity
        return v_start + (v_end - v_start) * (1 - star.radius/r)**self.beta


class LogisticVelocityAcceleration(AccelerationFunction):
    """ The velocity follows the Logistic (Sigmoid) Function """

    def __init__(self, steepness=10, r_mid=3.):
        super(LogisticVelocityAcceleration, self).__init__()
        self.steepness = steepness
        self.r_mid = r_mid

    def short(self, r, star):
        v_init = star.initial_wind_velocity
        v_end = star.terminal_wind_velocity

        r_mid = self.r_mid * star.radius
        exp = numpy.exp(-self.steepness * (r - r_mid) / (r_mid))

        return v_init, v_end, r_mid, exp

    def acceleration_from_radius(self, r, star):
        v_init, v_end, r_mid, exp = self.short(r, star)
        dvdr = (self.steepness * (v_end - v_init) * exp
                / (r_mid * (1. + exp)**2))
        v = v_init + (v_end - v_init) / (1. + exp)
        acc = v * dvdr
        return self.fix_acc_cutoff(r, acc, star)

    def velocity_from_radius(self, r, star):
        v_init, v_end, r_mid, exp = self.short(r, star)
        v = v_init + (v_end - v_init) / (1. + exp)
        return self.fix_v_cutoff(r, v, star)


class AGBAcceleration(AccelerationFunction):
    """ fit by Onno to the profiles by Nowotny 2005 """

    def __init__(self, alpha=10, r_mid=3.):
        super(AGBAcceleration, self).__init__()
        self.alpha = alpha
        self.r_mid = r_mid

    def short(self, r, star):
        v_init = star.initial_wind_velocity
        v_end = star.terminal_wind_velocity

        scaling = self.r_mid**self.alpha

        return v_init, v_end, scaling, r/star.radius

    def acceleration_from_radius(self, r, star):
        v_init, v_end, scaling, r_over_R = self.short(r, star)
        top = self.alpha * scaling * r_over_R**self.alpha
        bottom = r * (r_over_R**(-self.alpha) + scaling)**2
        dvdr = (v_end - v_init) * top / bottom
        acc = self.velocity_from_radius(r, star) * dvdr
        return self.fix_acc_cutoff(r, acc, star)

    def velocity_from_radius(self, r, star):
        v_init, v_end, scaling, r_over_R = self.short(r, star)

        denominator = 1 + scaling * r_over_R**(-self.alpha)
        v = v_init + (v_end - v_init)/denominator
        return self.fix_v_cutoff(r, v, star)


class AcceleratingWind(SimpleWind):
    """
       This wind model returns SPH particles moving away from the star at sub
       terminal velocity. It also adds a potential around the star that
       represents the radiation pressure. This potential can accelerate all
       particles away from the star using bridge. This is good for simulating
       processes within a few stellar radii.
    """

    acc_functions = {"constant_velocity": ConstantVelocityAcceleration,
                     "rsquared": RSquaredAcceleration,
                     "delayed_rsquared": DelayedRSquaredAcceleration,
                     "logistic": LogisticVelocityAcceleration,
                     "beta_law": BetaLawAcceleration,
                     "agb": AGBAcceleration,
                     }

    def __init__(self, *args, **kwargs):
        r_out_ratio = kwargs.pop("r_out_ratio", None)
        acc_start_ratio = kwargs.pop("acc_start_ratio", 2)
        grav_r_out_ratio = kwargs.pop("grav_r_out_ratio", None)
        acc_func = kwargs.pop("acceleration_function", "constant_velocity")
        acc_func_args = kwargs.pop("acceleration_function_args", {})
        self.critical_timestep = kwargs.pop("critical_timestep", None)
        self.v_init_ratio = kwargs.pop("v_init_ratio", None)
        self.compensate_pressure = kwargs.pop("compensate_pressure", False)

        super(AcceleratingWind, self).__init__(*args, **kwargs)

        if isinstance(acc_func, basestring):
            acc_func = self.acc_functions[acc_func]

        self.acc_function = acc_func(**acc_func_args)

        def r_out_function(r):
            if r_out_ratio is None:
                return None
            else:
                return r_out_ratio * r

        self.particles.add_calculated_attribute(
            "acc_cutoff", r_out_function,
            attributes_names=['radius'])

        def grav_r_out_function(r):
            if r_out_ratio is None:
                return None
            else:
                return grav_r_out_ratio * r

        self.particles.add_calculated_attribute(
            "grav_acc_cutoff", grav_r_out_function,
            attributes_names=['radius'])

        self.particles.add_calculated_attribute(
            "acc_start", lambda r: acc_start_ratio * r,
            attributes_names=['radius'])

        self.internal_energy_formula = self.scaled_u_from_T
        self.gamma = 5./3.

        self.staging_radius = None

    def set_initial_wind_velocity(self):
        if self.v_init_ratio is not None:
            self.particles.add_calculated_attribute(
                "initial_wind_velocity", lambda v: self.v_init_ratio * v,
                attributes_names=['terminal_wind_velocity'])

    def scaled_u_from_T(self, star, wind=None):
        """
            set the internal energy from the stellar surface temperature.
        """
        u_0 = (3./2. * constants.kB * star.temperature / star.mu)
        if wind is None:
            return u_0

        m_dot = star.wind_mass_loss_rate
        v_0 = star.initial_wind_velocity
        rho_0 = m_dot / (4. * numpy.pi * star.radius**2 * v_0)

        r = wind.position.lengths()
        v = wind.velocity.lengths()
        rho = m_dot / (4. * numpy.pi * r**2 * v)

        u = rho_0**(1 - self.gamma) * rho**(self.gamma - 1) * u_0
        return u

    def wind_sphere(self, star, Ngas):
        wind = Particles(Ngas)

        dt = (self.model_time - star.wind_release_time)
        if self.critical_timestep is None or dt > self.critical_timestep:
            acc_function = self.acc_function
        else:
            acc_function = ConstantVelocityAcceleration(use_initial=True)

        outer_wind_distance = acc_function.radius_from_time(dt, star)

        wind.position, direction = self.generate_positions(
            Ngas, star.radius, outer_wind_distance,
            acc_function.radius_from_number, star=star)

        velocities = acc_function.velocity_from_radius(
            wind.position.lengths(), star)
        wind.velocity = direction * self.as_three_vector(velocities)

        return wind

    def pressure_accelerations(self, indices, radii, star):
        v = self.acc_function.velocity_from_radius(radii, star)
        a = self.acc_function.acceleration_from_radius(radii, star)
        u = self.internal_energy_formula(star)
        m_dot = star.wind_mass_loss_rate
        v_init = star.initial_wind_velocity

        rho = m_dot / (4 * numpy.pi * v * radii**2)
        rho_init = m_dot / (4. * numpy.pi * v_init * star.radius**2)

        k = (self.gamma-1) * rho_init**(1-self.gamma) * u

        dvdr = a/v

        acceleration = (self.gamma * k * rho**(self.gamma-1)
                        * (2./radii + dvdr/v))

        return acceleration

    def radial_velocities(self, gas, star):
        rad_velocity = [] | units.ms
        pos_vel = zip(gas.position-star.position, gas.velocity-star.velocity)
        for pos, vel in pos_vel:
            rad_direction = pos/pos.length()
            scalar_projection = vel.dot(rad_direction)

            rad_velocity.append(scalar_projection)

        return rad_velocity

    def set_staging_radius(self, staging_radius, gas, timestep):
        self.staging_radius = staging_radius
        self.st_target_gas = gas
        self.st_timestep = timestep

    def staging_accelerations(self, indices, radii, star):
        particles = self.st_target_gas[indices]
        v_now = self.radial_velocities(particles, star)
        v_target = self.acc_function.velocity_from_radius(radii, star)
        dt = self.st_timestep
        acc = (v_target - v_now) / dt
        return acc

    def acceleration(self, star, radii):
        accelerations = numpy.zeros(radii.shape) | units.m/units.s**2

        i_acc = radii >= star.radius
        if star.acc_cutoff is not None:
            i_acc = i_acc & (radii < star.acc_cutoff)


        accelerations[i_acc] += self.acc_function.acceleration_from_radius(radii[i_acc], star)

        if self.compensate_pressure:
            if self.staging_radius is not None:
                i_pres = radii > star.radius * self.staging_radius
            else:
                i_pres = i_acc
            accelerations[i_pres] -= self.pressure_accelerations(
                i_pres, radii[i_pres], star)

        if self.compensate_gravity:
            if star.grav_acc_cutoff is not None:
                indices = radii < star.grav_acc_cutoff
                accelerations[indices] += constants.G * star.mass / radii[indices]**2
            else:
                accelerations += constants.G * star.mass / radii**2

        if self.staging_radius is not None:
            i_stag = radii < star.radius * self.staging_radius
            if i_stag.any():
                accelerations[i_stag] += self.staging_accelerations(i_stag, radii[i_stag], star)

        return accelerations

    def get_gravity_at_point(self, eps, x, y, z):
        total_acceleration = (
            numpy.zeros(shape=(len(x), 3)) | units.m/units.s**2)

        positions = quantities.as_vector_quantity([x, y, z]).transpose()
        for star in self.particles:
            relative_position = positions - star.position
            distance = relative_position.lengths()
            acceleration = self.acceleration(star, distance)

            direction = relative_position / self.as_three_vector(distance)
            direction = numpy.nan_to_num(direction)

            acc_vector = self.as_three_vector(acceleration)
            total_acceleration += direction * acc_vector

        return total_acceleration.transpose()


class HeatingWind(SimpleWind):
    """
        This wind model returns SPH particles that have no initial velocity
        with respect to the star. The energy of the integrated mechanical
        luminosity is added as internal energy. This is a numerical
        integration, so the timescale with which evolve_model is called should
        be small enough for convergence.

        This method good for simulating processes far from the star, and when
        the SPH particle mass is larger then the stellar mass loss per
        timestep.  It can make a big difference when the wind is derived from
        evolution.
    """

    def __init__(self, *args, **kwargs):
        self.feedback_efficiency = kwargs.pop("feedback_efficiency", 0.01)
        self.r_max_ratio = kwargs.pop("r_max_ratio", 5)
        super(HeatingWind, self).__init__(*args, **kwargs)

        self.internal_energy_formula = self.mechanical_internal_energy

        self.previous_time = 0 | units.Myr
        self.particles.track_mechanical_energy()

        self.supernova_energy = 1e51 | units.erg

    def evolve_particles(self):
        self.particles.evolve_mass_loss(self.model_time)

    def went_supernova(self, star, mass_lost):
        manual = star.mass_loss_type == "supernova"
        post_SN = star.stellar_type in [13, 14, 15] | units.stellar_type
        enough_lost = mass_lost > (.1 | units.MSun)

        return manual or (post_SN and enough_lost)

    def mechanical_internal_energy(self, star, wind):
        mass_lost = wind.mass.sum()
        lmech = star.mechanical_energy

        lmech_wind = lmech / (star.lost_mass/mass_lost + 1)
        star.mechanical_energy -= lmech_wind

        if self.went_supernova(star, mass_lost):
            lmech_wind = self.supernova_energy
            star.mass_loss_type = "wind"

        return self.feedback_efficiency * lmech_wind / mass_lost

    def wind_sphere(self, star, Ngas):
        wind = Particles(Ngas)

        r_max = self.r_max or self.r_max_ratio * star.radius
        wind.position, direction = self.generate_positions(Ngas, star.radius,
                                                           r_max)
        wind.velocity = [0, 0, 0] | units.kms

        return wind

    def reset(self):
        super(HeatingWind, self).reset()
        self.previous_time = 0 | units.Myr


def new_stellar_wind(sph_particle_mass, target_gas=None, timestep=None,
                     derive_from_evolution=False, mode="simple", **kwargs):
    """
        Create a new stellar wind code.
        target_gas: a gas particle set into which the wind particles should be
            put (requires timestep)
        timestep: the timestep at which the wind particles should be generated.
        derive_from_evolution: derive the wind parameters from stellar
            evolution (you still need to update the stellar parameters)
        mode: set to 'simple', 'accelerate' or 'mechanical'
    """
    if (target_gas is None) ^ (timestep is None):
        raise AmuseException("Must specify both target_gas and timestep"
                             "(or neither)")

    wind_modes = {"simple": SimpleWind,
                  "accelerate": AcceleratingWind,
                  "heating": HeatingWind,
                  }

    stellar_wind = wind_modes[mode](sph_particle_mass, derive_from_evolution,
                                    **kwargs)

    if target_gas is not None:
        stellar_wind.set_target_gas(target_gas, timestep)

    return stellar_wind


def static_wind_from_stellar_evolution(stellar_wind, stellar_evolution,
                                       start_time, end_time):
    """
        Convenience method that sets up the stellar wind parameters using a
        stellar evolution code. The change in the stars between start_time and
        end_time determines the stellar wind. Do not add the star particles to
        the stellar_wind code before calling this function.
    """
    stellar_evolution.evolve_model(start_time)
    stellar_wind.particles.add_particles(stellar_evolution.particles)

    stellar_evolution.evolve_model(end_time)
    chan = stellar_evolution.particles.new_channel_to(stellar_wind.particles)
    chan.copy_attributes(["age", "radius", "mass", "luminosity",
                          "temperature"])
    stellar_wind.evolve_model(0 | units.yr)
    stellar_wind.reset()

    return stellar_wind
