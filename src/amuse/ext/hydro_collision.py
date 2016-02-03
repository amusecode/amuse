import numpy
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles
from amuse.plot import plot, scatter, xlabel, ylabel, native_plot
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
    :argument initial_separation: a factor relative to the sum of the radii (1 means in contact, default: 5)
    """
    
    stellar_evolution_code_required = True
    gravity_code_required = True
    
    def __init__(
            self, 
            number_of_particles,
            hydrodynamics, 
            initial_separation = 5, 
            relax_sph_models = True,
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
        if not relax_sph_models:
            self.relax = self.no_relax
        self.verbose = verbose
        self.debug = debug
        self.hydrodynamics_arguments = hydrodynamics_arguments
        self.hydrodynamics_parameters = hydrodynamics_parameters
        self.star_to_sph_arguments = star_to_sph_arguments
        self.sph_to_star_arguments = sph_to_star_arguments
        
        self.dynamical_timescales_per_step = 1.0 # encounter_is_over check is performed at this interval
        self.extra_steps_when_encounter_is_over = 3
        
        self.continue_with_kepler = False
    
    def handle_collision(self, primary, secondary, stellar_evolution_code=None, gravity_code=None):
        particles = self.local_copy_of_particles(primary, secondary)
        self.collect_required_attributes(particles, gravity_code, stellar_evolution_code)
        self.backtrack_particles(particles)
        gas_particles = self.convert_stars(particles, stellar_evolution_code)
        self.simulate_collision(gas_particles)
        self.models = [convert_SPH_to_stellar_model(group, **self.sph_to_star_arguments) for group in self.groups_after_encounter]
        return self.new_particles_with_internal_structure_from_models()
    
    def new_particles_with_internal_structure_from_models(self):
        def get_internal_structure(set, particle=None):
            return self.models[(set.key == particle.key).nonzero()[0]]
        
        result = Particles(len(self.models))
        result.add_function_attribute("get_internal_structure", None, get_internal_structure)
        result.mass = [model.dmass.sum().as_quantity_in(self.mass_unit) for model in self.models]
        result.radius = [model.radius[-1].as_quantity_in(self.radius_unit) for model in self.models]
        result.position = (self.original_center_of_mass + self.stars_after_encounter.position).as_quantity_in(self.position_unit)
        result.velocity = (self.original_center_of_mass_velocity + self.stars_after_encounter.velocity).as_quantity_in(self.velocity_unit)
        return result
    
    def local_copy_of_particles(self, primary, secondary):
        particles = Particles(0)
        particles.add_particle(primary)
        particles.add_particle(secondary)
        return particles
    
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
        self.dynamical_timescale = numpy.pi * (particles.radius.sum()**3 / (8 * constants.G * particles.total_mass())).sqrt()
    
    def start_kepler(self, mass_unit, length_unit):
        unit_converter = nbody_system.nbody_to_si(mass_unit, length_unit)
        self.kepler = Kepler(unit_converter, redirection = "none" if self.debug else "null")
        self.kepler.initialize_code()
    
    def initialize_binary_in_kepler(self, star_a, star_b):
        self.kepler.initialize_from_dyn(
            star_a.mass + star_b.mass, 
            star_a.x - star_b.x, star_a.y - star_b.y, star_a.z - star_b.z,
            star_a.vx-star_b.vx, star_a.vy-star_b.vy, star_a.vz-star_b.vz
        )
        return self.kepler
    
    def backtrack_particles(self, particles):
        self.original_center_of_mass = particles.center_of_mass()
        self.original_center_of_mass_velocity = particles.center_of_mass_velocity()
        
        initial_separation = self.initial_separation * particles.radius.sum()
        if self.verbose:
            print "Particles at collision:"
            print particles
            print "Backtrack particles to initial separation", initial_separation.as_string_in(units.RSun)
        
        self.start_kepler(particles.total_mass(), initial_separation)
        kepler = self.initialize_binary_in_kepler(particles[0], particles[1])
        kepler.return_to_radius(initial_separation)
        self.begin_time = kepler.get_time()
        
        particles[1].position = kepler.get_separation_vector()
        particles[1].velocity = kepler.get_velocity_vector()
        kepler.advance_to_periastron()
        self.begin_time -= kepler.get_time()
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
            self.relax(convert_stellar_model_to_SPH(se_colliders[0], n_particles[0], **self.star_to_sph_arguments)),
            self.relax(convert_stellar_model_to_SPH(se_colliders[1], n_particles[1], **self.star_to_sph_arguments))
        )
        gas_particles = Particles()
        for particle, sph_model in zip(particles, sph_models):
            sph_model.position += particle.position
            sph_model.velocity += particle.velocity
            gas_particles.add_particles(sph_model)
        if self.verbose:
            print "Converting stars to SPH particles done"
        if self.debug:
            print gas_particles
        return gas_particles
    
    def divide_number_of_particles(self, particles):
        n1 = int(0.5 + self.number_of_particles * particles[0].mass / particles.total_mass())
        return (n1, self.number_of_particles - n1)
    
    def relax(self, sph_model):
        if self.debug:
            monitor = dict(time=[]|units.day, kinetic=[]|units.J, potential=[]|units.J, thermal=[]|units.J)
        gas_particles = sph_model.gas_particles
        hydro = self.new_hydrodynamics(gas_particles)
        hydro.parameters.artificial_viscosity_alpha = 0.0 # Viscous damping doesn't seem to be very important, but turned off just in case...
        channel_from_hydro = hydro.gas_particles.new_channel_to(gas_particles)
        channel_to_hydro = gas_particles.new_channel_to(hydro.gas_particles)
        
        dynamical_timescale = numpy.pi * (gas_particles.total_radius()**3 / (8 * constants.G * gas_particles.total_mass())).sqrt()
        t_end_in_t_dyn = 2.5 # Relax for this many dynamical timescales
        n_steps = 100
        velocity_damp_factor = 1.0 - (2.0*numpy.pi*t_end_in_t_dyn)/n_steps # Critical damping
        if self.verbose:
            print "Relaxing SPH model with {0} for {1} ({2} dynamical timescales).".format(
                self.hydrodynamics.__name__, 
                (t_end_in_t_dyn*dynamical_timescale).as_string_in(units.day),
                t_end_in_t_dyn)
        for i_step, time in enumerate(t_end_in_t_dyn*dynamical_timescale * numpy.linspace(1.0/n_steps, 1.0, n_steps)):
            hydro.evolve_model(time)
            channel_from_hydro.copy_attributes(["mass","x","y","z","vx","vy","vz","u"])
            gas_particles.position -= gas_particles.center_of_mass()
            gas_particles.velocity = velocity_damp_factor * (gas_particles.velocity - gas_particles.center_of_mass_velocity())
            channel_to_hydro.copy_attributes(["x","y","z","vx","vy","vz"])
            if self.debug:
                K, U, Q = hydro.kinetic_energy, hydro.potential_energy, hydro.thermal_energy
                print "t, K, U, Q:", time, K, U, Q
                monitor["time"].append(time)
                monitor["kinetic"].append(K)
                monitor["potential"].append(U)
                monitor["thermal"].append(Q)
                
        hydro.stop()
        if self.debug:
            energy_evolution_plot(monitor["time"], monitor["kinetic"], monitor["potential"], monitor["thermal"])
        return gas_particles
    
    def no_relax(self, sph_model):
        return sph_model.gas_particles
    
    def new_hop(self, particles):
        converter = nbody_system.nbody_to_si(particles.total_mass(), 1.0 | units.RSun)
        if self.debug:
            print "Output of Hop is redirected to hop_out.log"
            options = dict(redirection="file", redirect_file="hop_out.log")
        else:
            options = dict()
        hop = Hop(unit_converter=converter, **options)
        hop.parameters.number_of_neighbors_for_hop = 100
        hop.parameters.saddle_density_threshold_factor = 0.8
        hop.parameters.relative_saddle_density_threshold = True
        return hop
    
    def new_hydrodynamics(self, gas_particles):
        unit_converter = nbody_system.nbody_to_si(gas_particles.total_mass(), self.dynamical_timescale)
        hydro = self.hydrodynamics(unit_converter, **self.hydrodynamics_arguments)
        hydro.initialize_code()
        for par, value in self.hydrodynamics_parameters.iteritems():
            setattr(hydro.parameters, par, value)
        hydro.commit_parameters()
        hydro.gas_particles.add_particles(gas_particles)
        hydro.commit_particles()
        return hydro
    
    def simulate_collision(self, gas_particles):
        self.hop = self.new_hop(gas_particles)
        hydro = self.new_hydrodynamics(gas_particles)
        channel = hydro.gas_particles.new_channel_to(gas_particles)
        
        if self.verbose:
            print "Simulating collision with {0} from {1} to {2}.".format(
                self.hydrodynamics.__name__, 
                self.begin_time.as_string_in(units.day), 
                (self.dynamical_timescales_per_step * self.dynamical_timescale).as_string_in(units.day))
        
        hydro.evolve_model(self.dynamical_timescales_per_step * self.dynamical_timescale - self.begin_time)
        channel.copy_attributes(["x","y","z","vx","vy","vz","pressure","density","u"])
        extra_steps_counter = 0
        while True:
            if self.encounter_is_over(gas_particles):
                extra_steps_counter += 1
                if extra_steps_counter > self.extra_steps_when_encounter_is_over:
                    print "Encounter is over and finished extra steps."
                    break
                else:
                    print "Encounter is over. Now performing step {0} out of {1} extra steps".format(
                        extra_steps_counter, self.extra_steps_when_encounter_is_over)
            else:
                extra_steps_counter = 0
            print "Continuing to {0}.".format((hydro.model_time + self.next_dt + self.begin_time).as_string_in(units.day))
            if self.continue_with_kepler:
                self.evolve_with_kepler(hydro)
            hydro.evolve_model(hydro.model_time + self.next_dt)
            channel.copy_attributes(["x","y","z","vx","vy","vz","pressure","density","u"])
        
        hydro.stop()
        self.hop.stop()
        self.kepler.stop()
    
    def encounter_is_over(self, gas_particles):
        self.next_dt = self.dynamical_timescales_per_step * self.dynamical_timescale
        groups = self.group_bound_particles(gas_particles)
        stars = self.convert_groups_to_stars(groups)
        self.groups_after_encounter = groups
        self.stars_after_encounter = stars
        if len(stars) > 1:
            # Should do full check for stable binaries, triple, multiples, two escapers,
            # escaping star + binary, etc.
            # For now we only check whether the two most massive groups will (re)collide
            a, b = stars.sorted_by_attribute("mass")[-2:]
            if self.debug: print "System consists of {0} groups. The two most massive are: {1} and {2}.".format(len(stars), a.mass.as_string_in(units.MSun), b.mass.as_string_in(units.MSun))
            if self.binary_will_collide(a, b):
                return False
        
        if self.verbose:
            print "Encounter is over, {0} stars after encounter.".format(len(groups))
        return True
    
    def group_bound_particles(self, gas_particles):
        groups, lost = self.analyze_particle_distribution(gas_particles)
        while len(lost) > 0:
            if self.debug:
                group_plot(groups, lost)
            previous_number_of_lost_particles = len(lost)
            groups, lost = self.select_bound_particles(groups, lost)
            if len(lost) == previous_number_of_lost_particles:
                break
        return groups
    
    def convert_groups_to_stars(self, groups):
        stars = Particles(len(groups))
        for star, group in zip(stars, groups):
            star.mass = group.total_mass()
            star.position = group.center_of_mass()
            star.velocity = group.center_of_mass_velocity()
            star.radius = group.LagrangianRadii(mf=[0.9], cm=star.position)[0][0]
        return stars
    
    def analyze_particle_distribution(self, gas_particles):
        if self.verbose:
            print "Analyzing particle distribution using Hop"
        if "density" in gas_particles.get_attribute_names_defined_in_store():
            if self.debug: print "Using the original particles' density"
            self.hop.parameters.outer_density_threshold = 0.5 * gas_particles.density.mean()
            self.hop.particles.add_particles(gas_particles)
            gas_particles.copy_values_of_attribute_to("density", self.hop.particles)
        else:
            if self.debug: print "Using Hop to calculate the density"
            self.hop.particles.add_particles(gas_particles)
            self.hop.calculate_densities()
            self.hop.parameters.outer_density_threshold = 0.5 * self.hop.particles.density.mean()
        self.hop.do_hop()
        result = []
        for group in self.hop.groups():
            result.append(group.get_intersecting_subset_in(gas_particles))
        lost = self.hop.no_group().get_intersecting_subset_in(gas_particles)
        self.hop.particles.remove_particles(self.hop.particles)
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
    
    def binary_will_collide(self, a, b):
        self.continue_with_kepler = False
        if self.verbose:
            print "Using Kepler to check whether the two stars will (re)collide."
        kepler = self.initialize_binary_in_kepler(a, b)
        
        true_anomaly = kepler.get_angles()[1]
        eccentricity = kepler.get_elements()[1]
        if true_anomaly > 0.0 and eccentricity >= 1.0:
            if self.verbose:
                print "Stars are on hyperbolic/parabolic orbits and moving away from each other, interaction is over."
            return False
        
        periastron = kepler.get_periastron()
        will_collide = periastron < a.radius + b.radius
        if self.verbose:
            print "Stars {0} collide. Distance at periastron: {1}, sum of radii: {2}".format(
                "will" if will_collide else "won't",
                periastron.as_string_in(units.RSun), (a.radius + b.radius).as_string_in(units.RSun))
        
        if will_collide:
            # 1) check whether the stars are still relaxing: less than ~3 t_dyn passed since last moment of contact --> relax
            # 2) check whether the stars are already within 'initial_separation', else skip (dtmax?)
            
            kepler.advance_to_periastron()
            self.next_dt = kepler.get_time() + self.dynamical_timescales_per_step * self.dynamical_timescale
            if self.debug:
                print "Time to collision: {0}, next_dt: {1}".format(
                    kepler.get_time().as_string_in(units.day), self.next_dt.as_string_in(units.day))
            if kepler.get_time() > 3 * self.dynamical_timescale and kepler.get_apastron() > 2.0 * self.initial_separation * (a.radius + b.radius):
                # evolve for 3 * self.dynamical_timescale and skip the rest until ~initial_separation
                kepler.return_to_apastron()
                kepler.return_to_radius(a.radius + b.radius)
                if -kepler.get_time() > 2.9 * self.dynamical_timescale: # If ~3 t_dyn have passed since the end of the collision
                    if self.verbose: print "~3 t_dyn have passed since the end of the collision -> skip to next collision"
                    self.continue_with_kepler = True
                    kepler.advance_to_apastron()
                    kepler.advance_to_radius(2.0 * self.initial_separation * (a.radius + b.radius))
                    self.skip_to_relative_position_velocity = (kepler.get_separation_vector(), kepler.get_velocity_vector())
                    self.begin_time = kepler.get_time()
                    kepler.advance_to_periastron()
                    self.next_dt = self.dynamical_timescales_per_step * self.dynamical_timescale + kepler.get_time() - self.begin_time
                else:
                    self.next_dt = 3 * self.dynamical_timescale + kepler.get_time()
        return will_collide
    
    def evolve_with_kepler(self, hydro):
        if self.verbose: print "evolve_with_kepler"
        indices_two_most_massive = self.stars_after_encounter.mass.argsort()[-2:]
        groups = [self.groups_after_encounter[i] for i in indices_two_most_massive]
        old_particles = self.stars_after_encounter[indices_two_most_massive]
        new_particles = Particles(2)
        new_particles.mass = old_particles.mass
        new_particles[0].position, new_particles[0].velocity = self.skip_to_relative_position_velocity
        new_particles.move_to_center()
        for group, old_particle, new_particle in zip(groups, old_particles, new_particles):
            in_hydro = group.get_intersecting_subset_in(hydro.gas_particles)
            if self.verbose: print in_hydro.center_of_mass().as_quantity_in(units.RSun), old_particle.position.as_quantity_in(units.RSun), new_particle.position.as_quantity_in(units.RSun)
            in_hydro.position += new_particle.position - old_particle.position
            in_hydro.velocity += new_particle.velocity - old_particle.velocity
    
def group_plot(groups, no_group, figname="group_plot.png"):
    colors = ["r", "g", "b", "y", "k", "w"]*100
    for group, color in zip(groups, colors):
        scatter(group.x, group.y, c=color)
    
    if len(no_group):
        scatter(no_group.x, no_group.y, c="m", marker="s")
    
    native_plot.gca().set_aspect("equal", adjustable = "datalim")
    native_plot.savefig(figname)
    native_plot.clf()

def energy_evolution_plot(time, kinetic, potential, thermal, figname="energy_evolution.png"):
    native_plot.subplot(211)
    plot(time, kinetic, label='K')
    plot(time, potential, label='U')
    plot(time, thermal, label='Q')
    plot(time, kinetic + potential + thermal, label='E')
    xlabel('Time')
    ylabel('Energy')
    native_plot.legend(prop={'size':"x-small"}, loc=4)
    native_plot.subplot(212)
    plot(time, thermal, label='Q')
    native_plot.savefig(figname)
    native_plot.clf()
    

