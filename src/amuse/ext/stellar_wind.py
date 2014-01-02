import numpy

from amuse.support.exceptions import AmuseException
from amuse.datamodel.particles import Particles, ParticlesOverlay
from amuse.units import units, quantities, constants

from amuse.ext.evrard_test import uniform_unit_sphere

def as_three_vector(array):
    number = array
    if quantities.is_quantity(array):
        number = array.number
    three_vector = numpy.transpose([number]*3)
    if quantities.is_quantity(array):
        three_vector = three_vector | array.unit
    return three_vector

class StarsWithSimpleWind(object):
    """
        The simple wind model creates SPH particles moving away
        from the star at the terminal velocity.
        This is a safe assumption if the distance to other object
        is far larger then the stellar radius.
    """
    def __init__(self, stars, sph_particle_mass, evolving_stars=None):
        self.sph_particle_mass = sph_particle_mass
        self.derive_from_evolution = evolving_stars is not None
        self.system_time = 0.0|units.Myr

        self.setup_stars(stars, evolving_stars)
        self.set_mu()

    def initial_wind_velocity(self, stars):
        return stars.terminal_wind_velocity

    def update_stars(self):
        if not self.derive_from_evolution:
            return
        self.stev_channel.copy()
        self.stars.lost_mass = self.stars.init_mass - self.stars.mass - self.stars.released_mass

        """
            TODO: get the right formula for the terminal wind velocity.
            This one taken from Simon's example
        """
        t4 = numpy.log10(self.stars.temperature.value_in(units.K))-4
        t4 = t4.clip(0.,1.)
        self.stars.terminal_wind_velocity = (30. | units.kms) + ((4000. | units.kms)*t4)

    def setup_stars(self, stars, evolving_stars):
        if stars.can_extend_attributes():
            self.stars = stars
        else:
            self.stars = ParticlesOverlay(stars)

        if self.derive_from_evolution:
            self.stars.init_mass = self.stars.mass
            self.stars.released_mass = 0. | units.MSun
            try:
                evolving_particles = evolving_stars.particles
            except AttributeError:
                evolving_particles = evolving_stars

            self.stev_channel = evolving_particles.new_channel_to(self.stars,
                attributes=["mass", "radius", "temperature", "luminosity"])

            self.update_stars()
        else:
            self.stars.lost_mass = 0. | units.MSun

        self.stars.wind_release_time = 0. | units.Myr
        self.stars.wind_velocity = self.initial_wind_velocity(self.stars)

    def evolve_model(self, time):
        if not self.derive_from_evolution:
            elapsed_time = time - self.system_time
            self.stars.lost_mass += elapsed_time * self.stars.wind_mass_loss
        self.system_time = time

    def random_positions(self, N, rmin, rmax):
        """
            The particles start out in a random position between
            the surface of the star and the distance that the
            previously released particles have reached.
            This assumes that the current wind velocity is
            comparable to the previous wind velocity.

            Note that the stellar position is not added yet here.
        """
        random_positions = numpy.transpose(uniform_unit_sphere(N).make_xyz())
        vector_lengths = numpy.sqrt((random_positions**2).sum(1))

        unit_vectors = random_positions/as_three_vector(vector_lengths)

        distance = rmin + (vector_lengths * rmax)
        position = unit_vectors * as_three_vector(distance)

        return position, unit_vectors

    def set_mu(self, Y=0.25, Z=0.02, x_ion=0.1):
        X = 1.0 - Y - Z
        self.mu = constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

    def calculate_internal_energy(self, star, wind):
        return (3./2. * constants.kB * star.temperature / self.mu ) * 0.8

    def create_wind_particles_for_one_star(self, star):
        Ngas = int(star.lost_mass/self.sph_particle_mass)

        if Ngas == 0:
            return Particles(0)

        star.lost_mass -= Ngas * self.sph_particle_mass

        wind=Particles(Ngas)

        outer_wind_distance = star.wind_velocity * (self.system_time - star.wind_release_time)

        wind.position, direction = self.random_positions(Ngas, star.radius, outer_wind_distance)
        wind.velocity = direction * star.wind_velocity

        wind.u = self.calculate_internal_energy(star, wind)

        wind.mass = self.sph_particle_mass

        wind.position += star.position
        wind.velocity += star.velocity
        return wind

    def create_wind_particles(self):
        wind=Particles(0)

        self.update_stars()

        for star in self.stars:
            new_particles = self.create_wind_particles_for_one_star(star)
            if not new_particles.is_empty():
                wind.add_particles(new_particles)
                star.wind_release_time = self.system_time
                if self.derive_from_evolution:
                    star.released_mass += new_particles.mass.sum()

        return wind

    def has_new_wind_particles(self):
        self.update_stars()
        return self.stars.lost_mass.max() > self.sph_particle_mass

    def create_initial_wind(self, time=None, number=None, check_length=True):
        """
            This is a convenience method that creates some initial particles.
            They are created as if the wind has already been blowing for 'time'.

            If 'number' is given, the required time to get that number of particles
            is calculated. This assumes that the number of expected particles
            is far larger then the number of stars
        """
        if not number is None:
            required_mass = number * self.sph_particle_mass
            total_mass_loss = self.stars.wind_mass_loss.sum()
            time = 1.1 * required_mass/total_mass_loss

        self.evolve_model(time)

        wind = self.create_wind_particles()
        if check_length and len(wind) < 1:
            raise AmuseException("create_initial_wind time was too small to create any particles.")

        self.stars.lost_mass = 0.0|units.MSun
        self.stars.wind_release_time = 0.0|units.Myr
        self.system_time = 0.0|units.Myr

        return wind

    def get_gravity_at_point(self, eps, x, y, z):
        return [0, 0, 0]|units.m/units.s**2

    def get_potential_at_point(self, radius, x, y, z):
        return [0, 0, 0]|units.J

class StarsWithAcceleratingWind(StarsWithSimpleWind):
    """
       This wind model returns SPH particles moving away from the star at sub terminal velocity.
       It also adds a potential around the star that represents the radiation pressure.
       This potential can accelerate all particles away from the star using bridge.
    """

    def __init__(self, *args, **kwargs):
        self.init_v_wind_ratio = kwargs.pop("init_v_wind_ratio", 0.4)
        super(StarsWithAcceleratingWind, self).__init__(*args, **kwargs)
        self.r_min_ratio = 2
        self.r_max_ratio = 5

    def initial_wind_velocity(self, stars):
        return self.init_v_wind_ratio * stars.terminal_wind_velocity

    def wind_accelation_formula(self, star, distance):
        """
            A simple formula that approximates acceleration due to
            radiation pressure at a given distance from the star.
        """

        r_min = self.r_min_ratio * star.radius
        r_max = self.r_max_ratio * star.radius
        integrated_acceleration = 0.5 * ( star.terminal_wind_velocity**2 - star.wind_velocity**2 )
        scaling_constant = integrated_acceleration / (r_min**-1 - r_max**-1)

        radii_in_range = numpy.logical_and(distance > r_min, distance < r_max)

        acceleration = numpy.zeros(distance.shape) | units.m/units.s**2
        acceleration[radii_in_range] = scaling_constant * distance[radii_in_range]**-2

        return acceleration

    def get_gravity_at_point(self, eps, x, y, z):
        total_acceleration = numpy.zeros(shape=(len(x), 3))|units.m/units.s**2

        positions = quantities.as_vector_quantity(numpy.transpose([x, y, z]))
        for star in self.stars:
            relative_position = positions - star.position
            distance = relative_position.lengths()
            acceleration = self.wind_accelation_formula(star, distance)
            direction = relative_position / as_three_vector(distance)
            # Correct for directionless vectors with length 0
            direction[numpy.isnan(direction)] = 0
            total_acceleration += direction * as_three_vector(acceleration)

        return total_acceleration.transpose()

