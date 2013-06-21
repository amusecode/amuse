import numpy
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles
from amuse.plot import scatter, native_plot
#~from amuse.support.exceptions import AmuseException
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model
from amuse.community.kepler.interface import Kepler
from amuse.community.hop.interface import Hop


class StellarEncounterInHydrodynamics(object):
    """
    Resolves collisions between stars by converting them to SPH models, let them 
    collide in an SPH code, and converting the resulting SPH particle distribution 
    back to a 1D stellar evolution model.
    
    Requires a stellar evolution code to supply the internal structure of the 
    stars for the convert_stellar_model_to_SPH routine.
    Requires a gravity code to set up the initial configuration. The stars in the 
    gravity code have typically already collided, so they are first "evolved" back 
    in time up to a certain separation, assuming Keplerian motion.
    
    :argument number_of_particles: Total number of gas particles in the SPH simulation
    :argument hydrodynamics: SPH code class for the simulation
    :argument initial_separation: can be a physical length or a factor relative to 
        the sum of the radii (1 means in contact, default: 5)
    """
    
    stellar_evolution_code_required = True
    gravity_code_required = True
    
    def __init__(
            self, 
            number_of_particles,
            hydrodynamics, 
            initial_separation = 5, 
            verbose = False, 
            debug = False, 
            hydrodynamics_arguments = dict(),
            hydrodynamics_parameters = dict(),
            star_to_sph_arguments = dict(),
            sph_to_star_arguments = dict(),
        ):
        
        self.number_of_particles = number_of_particles
        self.hydrodynamics = hydrodynamics
        self.initial_separation = initial_separation
        self.verbose = verbose
        self.debug = debug
        self.hydrodynamics_arguments = hydrodynamics_arguments
        self.hydrodynamics_parameters = hydrodynamics_parameters
        self.star_to_sph_arguments = star_to_sph_arguments
        self.sph_to_star_arguments = sph_to_star_arguments
    
    def handle_collision(self, primary, secondary, stellar_evolution_code=None, gravity_code=None):
        particles = primary + secondary
        self.collect_required_attributes(particles, gravity_code, stellar_evolution_code)
        self.backtrack_particles(particles)
        gas_particles = self.convert_stars(particles, stellar_evolution_code)
        self.simulate_collision(gas_particles)
        
        self.models = []
        for group in self.groups_after_encounter:
            self.models.append(convert_SPH_to_stellar_model(group, **self.sph_to_star_arguments))
        
        def internal_structure(set, particle=None):
            return self.models[(set.key == particle.key).nonzero()[0]]
        
        result = Particles(len(self.models))
        result.add_function_attribute("internal_structure", None, internal_structure)
        result.mass = [model["dmass"].sum().as_quantity_in(self.mass_unit) for model in self.models]
        result.radius = [model["radius"][-1].as_quantity_in(self.radius_unit) for model in self.models]
        result.position = (self.original_center_of_mass + self.stars_after_encounter.position).as_quantity_in(self.position_unit)
        result.velocity = (self.original_center_of_mass_velocity + self.stars_after_encounter.velocity).as_quantity_in(self.velocity_unit)
        return result
    
    def collect_required_attributes(self, particles, gravity_code, stellar_evolution_code):
        # Collect the required attributes and copy to the particles in memory
        required_attributes = set(["mass", "x","y","z", "vx","vy","vz", "radius"])
        required_attributes -= set(particles.get_attribute_names_defined_in_store())
        for code in [stellar_evolution_code, gravity_code]:
            attrs_in_code = required_attributes & set(code.particles.get_attribute_names_defined_in_store())
            if len(attrs_in_code) > 0:
                code.particles.copy_values_of_attributes_to(list(attrs_in_code), particles)
                required_attributes -= attrs_in_code
        
        self.mass_unit = particles.mass.unit
        self.radius_unit = particles.radius.unit
        self.position_unit = particles.position.unit
        self.velocity_unit = particles.velocity.unit
        self.dynamical_timescale = (particles.radius.sum()**3 / (2 * constants.G * particles.total_mass())).sqrt()
    
    def backtrack_particles(self, particles):
        self.original_center_of_mass = particles.center_of_mass()
        self.original_center_of_mass_velocity = particles.center_of_mass_velocity()
        total_mass = particles.total_mass()
        relative_position = particles[1].position - particles[0].position
        relative_velocity = particles[1].velocity - particles[0].velocity
        
        if not hasattr(self.initial_separation, "unit") or self.initial_separation.unit is units.none:
            self.initial_separation *= particles.radius.sum()
        
        if self.verbose:
            print "Particles at collision:"
            print particles
            print "Backtrack particles to initial separation", self.initial_separation.as_string_in(units.RSun)
        
        unit_converter = nbody_system.nbody_to_si(total_mass, self.initial_separation)
        kepler = Kepler(unit_converter, redirection = "none" if self.debug else "null")
        kepler.initialize_code()
        kepler.initialize_from_dyn(
            total_mass, 
            relative_position[0], relative_position[1], relative_position[2],
            relative_velocity[0], relative_velocity[1], relative_velocity[2]
        )
        kepler.return_to_radius(self.initial_separation)
        self.begin_time = kepler.get_time()
        
        particles[1].position = kepler.get_separation_vector()
        particles[1].velocity = kepler.get_velocity_vector()
        kepler.stop()
        particles[0].position = [0, 0, 0] | units.m
        particles[0].velocity = [0, 0, 0] | units.m / units.s
        particles.move_to_center()
        if self.verbose:
            print "Backtracking particles done. Initial conditions:"
            print particles
    
    def convert_stars(self, particles, stellar_evolution_code):
        n_particles = self.divide_number_of_particles(particles)
        se_colliders = particles.get_intersecting_subset_in(stellar_evolution_code.particles)
        if self.verbose:
            print "Converting stars of {0} to SPH models of {1} particles, respectively.".format(particles.mass, n_particles)
        sph_models = (
            convert_stellar_model_to_SPH(se_colliders[0], n_particles[0], **self.star_to_sph_arguments),
            convert_stellar_model_to_SPH(se_colliders[1], n_particles[1], **self.star_to_sph_arguments)
        )
        gas_particles = Particles()
        for particle, sph_model in zip(particles, sph_models):
            sph_model.gas_particles.position += particle.position
            sph_model.gas_particles.velocity += particle.velocity
            gas_particles.add_particles(sph_model.gas_particles)
        if self.verbose:
            print "Converting stars to SPH particles done"
        if self.debug:
            print gas_particles
        return gas_particles
    
    def divide_number_of_particles(self, particles):
        n1 = int(0.5 + self.number_of_particles * particles[0].mass / particles.total_mass())
        return (n1, self.number_of_particles - n1)
    
    def simulate_collision(self, gas_particles):
        unit_converter = nbody_system.nbody_to_si(gas_particles.total_mass(), self.initial_separation)
        hydro = self.hydrodynamics(unit_converter, **self.hydrodynamics_arguments)
        hydro.initialize_code()
        for par, value in self.hydrodynamics_parameters.iteritems():
            setattr(hydro.parameters, par, value)
        hydro.commit_parameters()
        hydro.gas_particles.add_particles(gas_particles)
        hydro.commit_particles()
        channel = hydro.gas_particles.new_channel_to(gas_particles)
        
        if self.verbose:
            print "Simulating collision with {0} from {1} to {2}.".format(
                self.hydrodynamics.__name__, 
                self.begin_time.as_string_in(units.day), 
                (3*self.dynamical_timescale).as_string_in(units.day))
        
        hydro.evolve_model(3*self.dynamical_timescale - self.begin_time)
        channel.copy_attributes(["x","y","z","vx","vy","vz","pressure","density","u"])
        while not self.encounter_is_over(gas_particles):
            hydro.evolve_model(hydro.model_time + self.time_to_collision + 3*self.dynamical_timescale)
            channel.copy_attributes(["x","y","z","vx","vy","vz","pressure","density","u"])
        
        hydro.stop()
    
    def encounter_is_over(self, gas_particles):
        groups = self.group_bound_particles(gas_particles)
        stars = Particles(len(groups))
        for star, group in zip(stars, groups):
            star.mass = group.total_mass()
            star.radius = group.total_radius()
            star.position = group.center_of_mass()
            star.velocity = group.center_of_mass_velocity()
        
        if len(stars) > 1:
            # Should do full check for stable binaries, triple, multiples, two escapers,
            # escaping star + binary, etc.
            # For now we only check whether the two most massive groups will (re)collide
            a, b = stars.sorted_by_attribute("mass")[-2:]
            if self.binary_will_collide(a, b):
                return False
        
        self.groups_after_encounter = groups
        self.stars_after_encounter = stars
        if self.verbose:
            print "Encounter is over, {0} stars after encounter.".format(len(groups))
        return True
    
    def group_bound_particles(self, gas_particles):
        groups, lost = self.analyze_particle_distribution(gas_particles)
        while True:
            if self.debug:
                self.group_plot(groups, lost)
            previous_number_of_lost_particles = len(lost)
            groups, lost = self.select_bound_particles(groups, lost)
            if len(lost) == previous_number_of_lost_particles:
                break
        return groups
    
    def analyze_particle_distribution(self, gas_particles):
        if self.verbose:
            print "Analyzing particle distribution using Hop"
        converter = nbody_system.nbody_to_si(gas_particles.total_mass(), 1.0 | units.RSun)
        hop = Hop(unit_converter=converter, redirection = "none" if self.debug else "null")
        hop.parameters.density_method = 0
        hop.parameters.number_of_neighbors_for_hop = 100
        hop.parameters.number_of_neighbors_for_local_density = min(64, len(gas_particles) / 10)
        hop.particles.add_particles(gas_particles)
        hop.calculate_densities()
        hop.parameters.outer_density_threshold = 0.5 * hop.particles.density.mean()
        hop.parameters.saddle_density_threshold_factor = 0.8
        hop.parameters.relative_saddle_density_threshold = True
        
        hop.do_hop()
        result = []
        for group in hop.groups():
            result.append(group.get_intersecting_subset_in(gas_particles))
        lost = hop.no_group().get_intersecting_subset_in(gas_particles)
        hop.stop()
        return result, lost
    
    def select_bound_particles(self, groups, lost):
        specific_total_energy_relative_to_group = [] | (units.m / units.s)**2
        for group in groups:
            group_mass = group.total_mass()
            group_com = group.center_of_mass()
            group_com_velocity = group.center_of_mass_velocity()
            specific_total_energy_relative_to_group.append(
                (lost.velocity - group_com_velocity).lengths_squared() + lost.u - 
                constants.G * group_mass / (lost.position - group_com).lengths())
        index_minimum = specific_total_energy_relative_to_group.argmin(axis=0)
        bound=lost[:0]
        for i, group in enumerate(groups):
            bound_to_group = lost[numpy.logical_and(
                index_minimum == i, 
                specific_total_energy_relative_to_group[i] < 0 | (units.m / units.s)**2
            )]
            bound += bound_to_group
            groups[i] = group + bound_to_group
        return groups, lost - bound
    
    def group_plot(self, groups, no_group):
        colors = ["r", "g", "b", "y", "k", "w"]*100
        for group, color in zip(groups, colors):
            scatter(group.x, group.y, c=color)
        
        if len(no_group):
            scatter(no_group.x, no_group.y, c="m", marker="s")
        
        native_plot.gca().set_aspect("equal", adjustable = "datalim")
        native_plot.show()
    
    def binary_will_collide(self, a, b):
        if self.verbose:
            print "Using Kepler to check whether the two stars will (re)collide."
        total_mass = a.mass + b.mass
        relative_position = a.position - b.position
        relative_velocity = a.velocity - b.velocity
        unit_converter = nbody_system.nbody_to_si(total_mass, relative_position.length())
        kepler = Kepler(unit_converter, redirection = "none" if self.debug else "null")
        kepler.initialize_code()
        kepler.initialize_from_dyn(
            total_mass, 
            relative_position[0], relative_position[1], relative_position[2],
            relative_velocity[0], relative_velocity[1], relative_velocity[2]
        )
        true_anomaly = kepler.get_angles()[1]
        eccentricity = kepler.get_elements()[1]
        if true_anomaly > 0.0 and eccentricity >= 1.0:
            if self.verbose:
                print "Stars are on hyperbolic/parabolic orbits and moving away from each other, interaction is over."
            return False
        
        periastron = kepler.get_periastron()
        will_collide = periastron < a.radius + b.radius
        if will_collide:
            kepler.advance_to_periastron()
            self.time_to_collision = kepler.get_time()
        kepler.stop()
        if self.verbose:
            print "Stars {0} collide. Distance at periastron: {1}, sum of radii: {2}".format(
                "will" if will_collide else "won't",
                periastron.as_string_in(units.RSun), (a.radius + b.radius).as_string_in(units.RSun))
        return will_collide
    

