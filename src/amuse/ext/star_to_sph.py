import numpy
from amuse.support.data.core import Particles
from amuse.support.units import units, constants
from amuse.support.exceptions import AmuseWarning, AmuseException
from amuse.ext.spherical_model import new_spherical_particle_distribution, get_enclosed_mass_from_tabulated
from amuse.legacy.gadget2.interface import Gadget2
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.data.console import set_printing_strategy


class StellarModel2SPH(object):
    """
    Requests the internal structure of the star from a Stellar Evolution 
    legacy code and converts it into an SPH model consisting of the 
    specified number of particles. Useful for merging stars.
    
    :argument particle: Star particle to be converted to an SPH model
    :argument number_of_sph_particles: Number of gas particles in the resulting model
    :argument with_core_particle: Model the core as a heavy, non-sph particle (only for "scaling method")
    :argument do_relax: Relax the SPH model - doesn't seem to work satisfactorily yet!
    """
    
    def __init__(self, particle, number_of_sph_particles, seed = None, mode = "scaling method", 
            do_relax = False, sph_legacy_code = Gadget2, compatible_converter = ConvertBetweenGenericAndSiUnits,
            with_core_particle = False):
        self.particle = particle
        self.number_of_sph_particles = number_of_sph_particles
        self.zone_index = None
        self.delta      = None
        self.with_core_particle = with_core_particle
        numpy.random.seed(seed)
        if mode in ["random sampling", "scaling method"]:
            self.mode = mode  # "random sampling" or "scaling method"
        else:
            raise AmuseException("Unknown mode: {0}. Mode can be 'scaling method' or 'random sampling'.".format(mode))
        
        self.do_relax = do_relax
        self.sph_legacy_code = sph_legacy_code # used to relax the SPH model
        self.compatible_converter = compatible_converter
    
    def retrieve_stellar_structure(self):
        self.number_of_zones   = self.particle.get_number_of_zones().number
        self.number_of_species = self.particle.get_number_of_species().number
        self.species_names     = self.particle.get_names_of_species(number_of_species = self.number_of_species)
        self.species_IDs       = self.particle.get_IDs_of_species(number_of_species = self.number_of_species)
        self.frac_mass   = self.particle.get_mass_profile(number_of_zones = self.number_of_zones)
        self.density     = self.particle.get_density_profile(number_of_zones = self.number_of_zones)
        self.radius      = self.particle.get_radius_profile(number_of_zones = self.number_of_zones)
        self.temperature = self.particle.get_temperature_profile(number_of_zones = self.number_of_zones)
        self.mu          = self.particle.get_mu_profile(number_of_zones = self.number_of_zones)
        self.composition = self.particle.get_chemical_abundance_profiles(
            number_of_zones = self.number_of_zones, number_of_species = self.number_of_species)
        self.specific_internal_energy = (1.5 * constants.kB * self.temperature / self.mu).as_quantity_in(units.m**2/units.s**2)
        # Note: self.radius is in increasing order; from center to surface
        radius_profile = [0] | units.m
        radius_profile.extend(self.radius) # outer radius of each mesh zone
        self.midpoints = -(radius_profile[1:2])/2                           # dummy element to handle boundaries correctly
        self.midpoints.extend((radius_profile[1:] + radius_profile[:-1])/2) # real midpoints of each mesh zone
        self.midpoints.append(2*self.midpoints[-1] - self.midpoints[-2])   # dummy element to handle boundaries correctly
        
        if self.with_core_particle and self.mode == "scaling method":
            mean_density = 3.0 * self.particle.mass / (4.0 * numpy.pi * self.particle.radius**3)
            max_density = 1000 * mean_density
            # We assume reverse(self.density) is in ascending order.
            index = numpy.searchsorted(self.density[::-1], max_density)
            if index < self.number_of_zones:
                print "Will model core as a separate particle."
                i_core = self.number_of_zones - index
                self.core_radius = self.radius[i_core] - ((self.radius[i_core] - self.radius[i_core-1]) *
                    (max_density - self.density[i_core]) / (self.density[i_core-1] - self.density[i_core]))
                self.core_mass = get_enclosed_mass_from_tabulated(self.core_radius, radii = self.radius, densities = self.density)
                self.mass = self.particle.mass - self.core_mass
                print "core mass:", self.core_mass.as_quantity_in(units.MSun)
                set_printing_strategy("cgs")
                print self.radius[i_core-1], "<", self.core_radius, "<", self.radius[i_core]
                print self.density[i_core-1], ">", max_density, ">", self.density[i_core]
                set_printing_strategy("default")
                self.radius = self.radius[i_core:]
                self.density = self.density[i_core:]
                self.composition = self.composition[i_core:]
                self.specific_internal_energy = self.specific_internal_energy[i_core:]
                return
        
        self.core_radius = None
        self.core_mass = None
        self.mass = self.particle.mass
    
    def coordinates_from_spherical(self, radius, theta, phi):
        result  =      radius * numpy.sin( theta ) * numpy.cos( phi )   # x
        result.extend( radius * numpy.sin( theta ) * numpy.sin( phi ) ) # y
        result.extend( radius * numpy.cos( theta ) )                    # z
        return result
    
    def set_zone_indices_and_interpolation_coeffs(self):
        self.zone_index = numpy.empty(self.number_of_sph_particles, dtype="int32")
        self.delta      = numpy.empty(self.number_of_sph_particles)
        rands     = numpy.random.uniform(size = self.number_of_sph_particles)
        self.rands = rands
        selection = numpy.arange(self.number_of_sph_particles, dtype='int32') # used for indexing
        enclosed_mass = 1.0 # fractional enclosed mass, at the surface initially
        for (i, f_mass_i) in enumerate(self.frac_mass.value_in(units.none)):
            enclosed_mass -= f_mass_i
            found = numpy.where( enclosed_mass < rands[selection])[0]
            self.zone_index[selection[found]] = i
            self.delta[selection[found]] = (rands[selection[found]] - enclosed_mass) / f_mass_i
            selection = numpy.delete(selection, found)
            if not len(selection): break
        if len(selection):
            raise AmuseException("Some zone indices could not be determined")
        self.delta = units.none.new_quantity(self.delta)
    
    def new_radial_positions(self):
        if self.zone_index is None:
            self.set_zone_indices_and_interpolation_coeffs()
        radius_profile = self.radius
        radius_profile.append(0|units.m)
        radial_positions = (   self.delta    * radius_profile[self.zone_index] + 
                            (1 - self.delta) * radius_profile[self.zone_index+1] )
        return radial_positions
    
    def new_positions_spherical_coordinates(self):
        radius = self.new_radial_positions()
        theta  = numpy.arccos(numpy.random.uniform(-1.0, 1.0, self.number_of_sph_particles))
        phi    = numpy.random.uniform(0.0, 2.0*numpy.pi, self.number_of_sph_particles)
        return (radius,theta,phi)
    
    def new_positions(self):
        return self.coordinates_from_spherical(*self.new_positions_spherical_coordinates())
    
    def get_index(self, value, sorted_vector):
        if not sorted_vector[0] <= value <= sorted_vector[-1]:
            raise AmuseException("Can't find a valid index. {0} is not in "
                "the range [{1}, {2}].".format(value, sorted_vector[0], sorted_vector[-1]))
        index = numpy.searchsorted(sorted_vector, value)
        return max(index - 1, 0)
    
    def calculate_interpolation_coefficients(self, radial_positions):
        indices = numpy.array([self.get_index(r, self.midpoints) for r in radial_positions])
        delta = (self.midpoints[indices+1] - radial_positions) / (self.midpoints[indices+1] - self.midpoints[indices])
        return indices, delta
    
    def interpolate_internal_energy(self, radial_positions, do_composition_too = True):
        indices, delta = self.calculate_interpolation_coefficients(radial_positions)
        one_minus_delta = 1 - delta
        
        extended = self.specific_internal_energy[:1]
        extended.extend(self.specific_internal_energy)
        extended.append(self.specific_internal_energy[-1])
        interpolated_energies = delta*extended[indices] + one_minus_delta*extended[indices+1]
        
        if do_composition_too:
            comp = [] | units.none
            for species in self.composition:
                extended = species[:1]
                extended.extend(species)
                extended.append(species[-1])
                comp.append(delta*extended[indices] + one_minus_delta*extended[indices+1])
            return interpolated_energies, comp.transpose()
        else:
            return interpolated_energies
    
    def convert_to_SPH(self):
        if self.mode == "scaling method":
            sph_particles = new_spherical_particle_distribution(
                self.number_of_sph_particles, 
                radii = self.radius, densities = self.density, core_radius = self.core_radius)
        else:
            sph_particles = Particles(self.number_of_sph_particles)
            sph_particles.position = self.new_positions()
        sph_particles.mass = (self.mass.number * 1.0 / 
            self.number_of_sph_particles) | self.mass.unit
        sph_particles.velocity = [0,0,0] | units.m/units.s
        return sph_particles
    
    def relax(self, particles):
        num_iterations = 20
        max_delta = 0.01  # maximum change to particle positions relative to its smoothing length
        
        result = []
        previous_acc = 0 | units.m / units.s**2
        unit_converter = self.compatible_converter(self.particle.radius, self.particle.mass, 1.0e-3 | units.s)
        hydro_legacy_code = self.sph_legacy_code(unit_converter)
        particles.u = 1.0 | (units.m / units.s)**2
        hydro_legacy_code.gas_particles.add_particles(particles)
        
        for i in range(1, num_iterations+1):
            hydro_legacy_code.gas_particles.u = self.interpolate_internal_energy(particles.position.lengths(), do_composition_too = False)
            hydro_legacy_code.evolve_model(i * (1.0e-5 | units.s))
            accelerations     = hydro_legacy_code.gas_particles.acceleration
            acc_correlated = (previous_acc * accelerations).sum() / (accelerations * accelerations).sum()
            if (acc_correlated < 0.5 | units.none and i > 2):
                break
            
            previous_acc = accelerations
            internal_energies = hydro_legacy_code.gas_particles.u
            smoothing_lengths = hydro_legacy_code.gas_particles.h_smooth
            factor = numpy.minimum((max_delta * internal_energies / (accelerations.lengths() * smoothing_lengths)).value_in(units.none), 0.5)
            result.append(str(i) + ": Accelerations correlated: " + str(acc_correlated) + ", median factor: " + str(numpy.median(factor)))
            
            particles.position += accelerations * ((smoothing_lengths * smoothing_lengths * factor) / 
                internal_energies).reshape((self.number_of_sph_particles, 1))
            hydro_legacy_code.gas_particles.position = particles.position
            hydro_legacy_code.gas_particles.velocity = particles.velocity
        
        particles.u = hydro_legacy_code.gas_particles.u
        hydro_legacy_code.stop()
        if i == num_iterations:
            print "\nUnable to converge to stable SPH model within {0} iterations.".format(num_iterations)
        else:
            print "\nSuccessfully converged to stable SPH model within {0} iterations.".format(i-1)
        return result
    
    @property
    def result(self):
        self.retrieve_stellar_structure()
        sph_particles = self.convert_to_SPH()
        if self.do_relax:
            for result_string in self.relax(sph_particles):
                print result_string
        
        sph_particles.add_vector_attribute("composition", self.species_names)
        specific_internal_energy, composition = self.interpolate_internal_energy(sph_particles.position.lengths())
        sph_particles.u = specific_internal_energy
        sph_particles.composition = composition
        if self.with_core_particle:
            if self.core_radius:
                core_particle = Particles(1)
                core_particle.mass = self.core_mass
                core_particle.position = [0.0, 0.0, 0.0] | units.m
                core_particle.velocity = [0.0, 0.0, 0.0] | units.m / units.s
                return core_particle, sph_particles, self.core_radius
            else:
                return Particles(), sph_particles, None
        return sph_particles
    

def convert_stellar_model_to_SPH(particle, number_of_sph_particles, **keyword_arguments):
    """
    Requests the internal structure of the star from a Stellar Evolution 
    legacy code and converts it into an SPH model consisting of the 
    specified number of particles. Useful for merging stars.
    
    :argument particle: Star particle to be converted to an SPH model
    :argument number_of_sph_particles: Number of gas particles in the resulting model
    """
    converter = StellarModel2SPH(particle, number_of_sph_particles, **keyword_arguments)
    return converter.result
