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

class StarsWithMassLoss(Particles):
    def __init__(self, *args, **kwargs):
        super(StarsWithMassLoss, self).__init__(*args, **kwargs)
        self.collection_attributes.timestamp = 0. | units.yr
        self.collection_attributes.previous_time = 0. | units.yr

    def add_particles(self, *args, **kwargs):
        new_particles = super(StarsWithMassLoss, self).add_particles(*args, **kwargs)
        new_particles.lost_mass = 0. | units.MSun
        new_particles.wind_release_time = self.collection_attributes.timestamp
        new_particles.mu = self.collection_attributes.global_mu

        return new_particles

    def evolve_mass_loss(self, time):
        if self.collection_attributes.previous_time < time:
            elapsed_time = time - self.collection_attributes.previous_time
            self.lost_mass += elapsed_time * self.wind_mass_loss

            self.collection_attributes.timestamp = time

            self.collection_attributes.previous_time = time

    def set_global_mu(self, mu):
        self.mu = mu
        self.collection_attributes.global_mu = mu

    def reset(self):
        self.lost_mass = 0.0|units.MSun
        self.wind_release_time = 0.0|units.yr
        self.collection_attributes.timestamp = 0. | units.yr
        self.collection_attributes.previous_time = 0. | units.yr

class EvolvingStarsWithMassLoss(StarsWithMassLoss):
    """
        Derive the stellar wind from stellar evolution.
        You have to copy the relevant attributes from the stellar evolution.
        This can be done using a channel like:

        chan = stellar_evolution.particles.new_channel_to(stellar_wind.particles,
            attributes=["age", "radius", "mass", "luminosity", "temperature"])

        while <every timestep>:
            chan.copy()
    """
    def add_particles(self, *args, **kwargs):
        new_particles = super(EvolvingStarsWithMassLoss, self).add_particles(*args, **kwargs)
        new_particles.wind_mass_loss = 0. | units.MSun/units.yr
        new_particles.terminal_wind_velocity = 0. | units.ms
        new_particles.previous_age = new_particles.age
        new_particles.previous_mass = new_particles.mass
        return new_particles

    def evolve_mass_loss(self, time):
        if self.collection_attributes.previous_time < time:
            self.update_from_evolution()
            StarsWithMassLoss.evolve_mass_loss(self, time)

    def calculate_terminal_wind_velocity(self, Y=0.25, I_He = 2):
        """
          This routine calculates the escape and terminal wind velocity. The Equations are taken
          from Kudritzki & Puls, Annual Reviews of Astronomy and Astrophysics, 2000, Vol. 38,
          p.613-666 Equation (8) and (9) and Kudritzki et al., 1989, A&A 219, 205 Equation (64)
          and (65).

          I_He:    Number of electrons per He nucleus (= 2 in O-Stars)
          sigma_e:   Thomson absorption coefficient
          Gamma:     Ratio of radiative Thomson to gravitational acceleration
        """
        T_eff = self.temperature

        sigma_e = 0.398 * (1 + I_He*Y)/(1 + 4*Y)
        Gamma = 7.66E-5 * sigma_e * self.luminosity.value_in(units.LSun)/ self.mass.value_in(units.MSun)
        g_star = constants.G * self.mass / self.radius**2
        v_esc = (2*g_star*self.radius*(1 - Gamma))**0.5

        condlist = [T_eff >= 21000. | units.K, (10000. | units.K  < T_eff) &  (T_eff < 21000. | units.K), T_eff <= 10000. | units.K]
        choicelist = [2.65, 1.4, 1.0]

        return v_esc * numpy.select(condlist, choicelist)

    def update_from_evolution(self):
        if (self.age != self.previous_age).any():
            mass_loss = self.previous_mass - self.mass
            timestep = self.age - self.previous_age
            self.wind_mass_loss = mass_loss / timestep

            self.terminal_wind_velocity = self.calculate_terminal_wind_velocity()

    def reset(self):
        super(EvolvingStarsWithMassLoss, self).reset()
        self.previous_age = 0|units.yr

class SimpleWind(object):
    """
        The simple wind model creates SPH particles moving away
        from the star at the terminal velocity.
        This is a safe assumption if the distance to other objects
        is far larger then the stellar radius.
    """
    def __init__(self, sph_particle_mass, derive_from_evolution=False):
        self.sph_particle_mass = sph_particle_mass
        self.system_time = 0.0|units.yr

        if derive_from_evolution:
            self.particles = EvolvingStarsWithMassLoss()
        else:
            self.particles = StarsWithMassLoss()
        self.target_gas = self.timestep = None

        self.set_global_mu()
        self.internal_energy_formula = self.calculate_internal_energy_from_temperature

    def initial_wind_velocity(self, stars):
        return stars.terminal_wind_velocity

    def evolve_model(self, time):
        if self.has_target():
            while self.system_time < time:
                self.particles.evolve_mass_loss(self.system_time)
                if self.has_new_wind_particles():
                    wind_gas = self.create_wind_particles()
                    self.target_gas.add_particles(wind_gas)
                self.system_time += self.timestep
        else:
            self.system_time = time
            self.particles.evolve_mass_loss(self.system_time)

    def set_target_gas(self, target_gas, timestep):
        self.target_gas = target_gas
        self.timestep = timestep

    def has_target(self):
        return self.target_gas is not None

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

    def set_global_mu(self, Y=0.25, Z=0.02, x_ion=0.1):
        """
            Set the global value of mu used to create stellar wind.
            If the value of mu is known directly,
            use <stellar_wind>.particles.set_global_mu().
            An alternative way is to set mu for each star separately.
        """
        X = 1.0 - Y - Z
        mu = constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)
        self.particles.set_global_mu(mu)

    def calculate_internal_energy_from_temperature(self, star, wind):
        return (3./2. * constants.kB * star.temperature / star.mu ) * 0.8

    def create_wind_particles_for_one_star(self, star):
        Ngas = int(star.lost_mass/self.sph_particle_mass)

        if Ngas == 0:
            return Particles(0)

        star.lost_mass -= Ngas * self.sph_particle_mass

        wind=Particles(Ngas)

        wind_velocity = self.initial_wind_velocity(star)
        outer_wind_distance = wind_velocity * (self.system_time - star.wind_release_time)

        wind.position, direction = self.random_positions(Ngas, star.radius, outer_wind_distance)
        wind.velocity = direction * wind_velocity

        wind.u = self.internal_energy_formula(star, wind)

        wind.mass = self.sph_particle_mass

        wind.position += star.position
        wind.velocity += star.velocity
        return wind

    def create_wind_particles(self):
        wind=Particles(0)

        for star in self.particles:
            new_particles = self.create_wind_particles_for_one_star(star)
            if not new_particles.is_empty():
                wind.add_particles(new_particles)
                star.wind_release_time = self.system_time

        return wind

    def has_new_wind_particles(self):
        return self.particles.lost_mass.max() > self.sph_particle_mass

    def create_initial_wind(self, time=None, number=None, check_length=True):
        """
            This is a convenience method that creates some initial particles.
            They are created as if the wind has already been blowing for 'time'.
            Note that this does not work if the mass loss is derived from stellar evolution.

            If 'number' is given, the required time to get that number of particles
            is calculated. This assumes that the number of expected particles
            is far larger then the number of stars
        """
        if not number is None:
            required_mass = number * self.sph_particle_mass
            total_mass_loss = self.particles.wind_mass_loss.sum()
            time = 1.1 * required_mass/total_mass_loss

        if self.has_target():
            start_target_gas = self.target_gas.copy()

        self.evolve_model(time)

        if self.has_target():
            wind = self.target_gas - start_target_gas
        else:
            wind = self.create_wind_particles()

        if check_length and len(wind) < 1:
            raise AmuseException("create_initial_wind time was too small to create any particles.")

        self.particles.reset()
        self.system_time = 0.0|units.yr

        return wind

    def get_gravity_at_point(self, eps, x, y, z):
        return [0, 0, 0]|units.m/units.s**2

    def get_potential_at_point(self, radius, x, y, z):
        return [0, 0, 0]|units.J

class AcceleratingWind(SimpleWind):
    """
       This wind model returns SPH particles moving away from the star at sub terminal velocity.
       It also adds a potential around the star that represents the radiation pressure.
       This potential can accelerate all particles away from the star using bridge.
    """

    def __init__(self, *args, **kwargs):
        super(AcceleratingWind, self).__init__(*args, **kwargs)
        self.init_v_wind_ratio = 0.4
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
        init_wind_velocity = self.initial_wind_velocity(star)
        integrated_acceleration = 0.5 * ( star.terminal_wind_velocity**2 - init_wind_velocity**2 )
        scaling_constant = integrated_acceleration / (r_min**-1 - r_max**-1)

        radii_in_range = numpy.logical_and(distance > r_min, distance < r_max)

        acceleration = numpy.zeros(distance.shape) | units.m/units.s**2
        acceleration[radii_in_range] = scaling_constant * distance[radii_in_range]**-2

        return acceleration

    def get_gravity_at_point(self, eps, x, y, z):
        total_acceleration = numpy.zeros(shape=(len(x), 3))|units.m/units.s**2

        positions = quantities.as_vector_quantity(numpy.transpose([x, y, z]))
        for star in self.particles:
            relative_position = positions - star.position
            distance = relative_position.lengths()
            acceleration = self.wind_accelation_formula(star, distance)
            direction = relative_position / as_three_vector(distance)
            # Correct for directionless vectors with length 0
            direction[numpy.isnan(direction)] = 0
            total_acceleration += direction * as_three_vector(acceleration)

        return total_acceleration.transpose()

def new_stellar_wind(sph_particle_mass, target_gas=None, timestep=None, derive_from_evolution=False, accelerate=False):
    """
        Create a new stellar wind code.
        target_gas: a gas particle set into which the wind particles should be put (requires timestep)
        timestep: the timestep at which the wind particles should be generated.
        derive_from_evolution: derive the wind parameters from stellar evolution (you still need to update the stellar parameters)
        accelerate: start at subterminal velocity and accelerate the gas near the star (requires a bridge coupling)

        TODO: add option to setup non evolving wind from stellar evolution.
    """
    if (target_gas is None) ^ (timestep is None):
        raise AmuseException("Must specify both target_gas and timestep (or neither)")

    if accelerate:
        stellar_wind = AcceleratingWind(sph_particle_mass, derive_from_evolution)
    else:
        stellar_wind = SimpleWind(sph_particle_mass, derive_from_evolution)

    if target_gas is not None:
        stellar_wind.set_target_gas(target_gas, timestep)

    return stellar_wind
