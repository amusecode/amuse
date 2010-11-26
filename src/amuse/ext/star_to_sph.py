import numpy
import pickle
import os.path
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
            with_core_particle = False, pickle_file = None):
        self.particle = particle
        self.number_of_sph_particles = number_of_sph_particles
        self.with_core_particle = with_core_particle
        self.pickle_file = pickle_file
        if seed:
            numpy.random.seed(seed)
        if mode in ["random sampling", "scaling method"]:
            self.mode = mode  # "random sampling" or "scaling method"
            self.zone_index = None
            self.delta      = None
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
        self.frac_mass_profile   = self.particle.get_mass_profile(number_of_zones = self.number_of_zones)
        self.density_profile     = self.particle.get_density_profile(number_of_zones = self.number_of_zones)
        self.radius_profile      = self.particle.get_radius_profile(number_of_zones = self.number_of_zones)
        self.temperature_profile = self.particle.get_temperature_profile(number_of_zones = self.number_of_zones)
        self.mu_profile          = self.particle.get_mu_profile(number_of_zones = self.number_of_zones)
        self.composition_profile = self.particle.get_chemical_abundance_profiles(
            number_of_zones = self.number_of_zones, number_of_species = self.number_of_species)
        self.specific_internal_energy_profile = (1.5 * constants.kB * self.temperature_profile / self.mu_profile).as_quantity_in(units.m**2/units.s**2)
        # Note: self.radius is in increasing order; from center to surface
        radius_profile = [0] | units.m
        radius_profile.extend(self.radius_profile) # outer radius of each mesh zone
        self.midpoints_profile = -(radius_profile[1:2])/2                           # dummy element to handle boundaries correctly
        self.midpoints_profile.extend((radius_profile[1:] + radius_profile[:-1])/2) # real midpoints of each mesh zone
        self.midpoints_profile.append(2*self.midpoints_profile[-1] - self.midpoints_profile[-2])   # dummy element to handle boundaries correctly
        
        self.mass         = self.particle.mass
        self.radius       = self.particle.radius
    
    def unpickle_stellar_structure(self):
        if os.path.isfile(self.pickle_file):
            infile = open(self.pickle_file, 'rb')
        else:
            raise AmuseException("Input pickle file '{0}' does not exist".format(self.pickle_file))
        structure = pickle.load(infile)
        self.mass   = structure['mass']
        self.radius = structure['radius']
        self.number_of_zones     = structure['number_of_zones']
        self.number_of_species   = structure['number_of_species']
        self.species_names       = structure['species_names']
        self.species_IDs         = structure['species_IDs']
        self.frac_mass_profile   = structure['frac_mass_profile']
        self.density_profile     = structure['density_profile']
        self.radius_profile      = structure['radius_profile']
        self.temperature_profile = structure['temperature_profile']
        self.mu_profile          = structure['mu_profile']
        self.composition_profile = structure['composition_profile']
        self.specific_internal_energy_profile = structure['specific_internal_energy_profile']
        self.midpoints_profile   = structure['midpoints_profile']
    
    def setup_core_parameters(self):
        self.core_radius = None
        self.core_mass = None
        self.sph_core_radius = None
        
        if self.with_core_particle and self.mode == "scaling method":
            mean_density = 3.0 * self.mass / (4.0 * numpy.pi * self.radius**3)
            max_density = 1000 * mean_density
            # We assume reverse(self.density) is in ascending order.
            index = numpy.searchsorted(self.density_profile[::-1], max_density)
            if index < self.number_of_zones:
                print "Will model core as a separate particle."
                i_core = self.number_of_zones - index
                self.core_radius = self.radius_profile[i_core] - ((self.radius_profile[i_core] - self.radius_profile[i_core-1]) *
                    (max_density - self.density_profile[i_core]) / (self.density_profile[i_core-1] - self.density_profile[i_core]))
                self.core_mass = get_enclosed_mass_from_tabulated(self.core_radius, 
                    radii = self.radius_profile, densities = self.density_profile)
                print "core mass:", self.core_mass.as_quantity_in(units.MSun)
                print "core radius:", self.core_radius.as_quantity_in(units.RSun)
                if True:
                    self.density_profile[:i_core] = max_density
                    self.core_mass -= get_enclosed_mass_from_tabulated(self.core_radius, 
                        radii = self.radius_profile, densities = self.density_profile)
                    self.sph_core_radius = 0 | units.m
                    print "core mass in DM particle:", self.core_mass.as_quantity_in(units.MSun)
                else:
                    self.radius_profile = self.radius_profile[i_core:]
                    self.density_profile = self.density_profile[i_core:]
                    self.composition_profile = self.composition_profile[i_core:]
                    self.specific_internal_energy_profile = self.specific_internal_energy_profile[i_core:]
                    self.sph_core_radius = self.core_radius
                
                self.mass = self.mass - self.core_mass
    
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
        for (i, f_mass_i) in enumerate(self.frac_mass_profile.value_in(units.none)):
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
        radius_profile = self.radius_profile
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
        indices = numpy.array([self.get_index(r, self.midpoints_profile) for r in radial_positions])
        delta = (self.midpoints_profile[indices+1] - radial_positions) / (
            self.midpoints_profile[indices+1] - self.midpoints_profile[indices])
        return indices, delta
    
    def interpolate_internal_energy(self, radial_positions, do_composition_too = True):
        indices, delta = self.calculate_interpolation_coefficients(radial_positions)
        one_minus_delta = 1 - delta
        
        extended = self.specific_internal_energy_profile[:1]
        extended.extend(self.specific_internal_energy_profile)
        extended.append(self.specific_internal_energy_profile[-1])
        interpolated_energies = delta*extended[indices] + one_minus_delta*extended[indices+1]
        
        if do_composition_too:
            comp = [] | units.none
            for species in self.composition_profile:
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
                radii = self.radius_profile, densities = self.density_profile, core_radius = self.sph_core_radius)
        else:
            sph_particles = Particles(self.number_of_sph_particles)
            sph_particles.position = self.new_positions()
        sph_particles.mass = (self.mass.number * 1.0 / 
            self.number_of_sph_particles) | self.mass.unit
        sph_particles.velocity = [0,0,0] | units.m/units.s
        # Crude estimate of the smoothing length; the SPH code will calculate the true value itself.
        sph_particles.h_smooth = self.radius * (self.number_of_sph_particles/50.0)**(-1/3.0)
        return sph_particles
    
    def relax(self, particles):
        num_iterations = 20
        max_delta = 0.01  # maximum change to particle positions relative to its smoothing length
        
        result = []
        previous_acc = 0 | units.m / units.s**2
        unit_converter = self.compatible_converter(self.radius, self.mass, 1.0e-3 | units.s)
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
        if self.pickle_file is None:
            self.retrieve_stellar_structure()
        else:
            self.unpickle_stellar_structure()
        self.setup_core_parameters()
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
                core_particle.radius = 0.0 | units.m
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

def pickle_stellar_model(particle, pickle_file_name):
    """
    Requests the internal structure of the star from a Stellar Evolution legacy 
    code and pickles it (stores it as a *.pkl file), for later use:
    convert_stellar_model_to_SPH(None, ..., pickle_file = pickle_file_name)
    Using a pickled stellar model is significantly faster for modelling giants 
    and other extremely evolved stars.
    
    :argument particle: Star particle to be converted to an SPH model later
    :argument pickle_file_name: Name of the pickle file in which to store the stellar structure
    """
    if os.path.isdir(os.path.dirname(pickle_file_name)) and not os.path.exists(pickle_file_name):
        outfile = open(pickle_file_name, 'wb')
    else:
        raise AmuseWarning("Incorrect file name '{0}'; directory must exist and "
            "file may not exist".format(pickle_file_name))
    converter = StellarModel2SPH(particle, None)
    converter.retrieve_stellar_structure()
    pickle.dump(dict(  
        mass                = converter.mass,
        radius              = converter.radius,
        number_of_zones     = converter.number_of_zones,
        number_of_species   = converter.number_of_species,
        species_names       = converter.species_names,
        species_IDs         = converter.species_IDs,
        frac_mass_profile   = converter.frac_mass_profile,
        density_profile     = converter.density_profile,
        radius_profile      = converter.radius_profile,
        temperature_profile = converter.temperature_profile,
        mu_profile          = converter.mu_profile,
        composition_profile = converter.composition_profile,
        specific_internal_energy_profile = converter.specific_internal_energy_profile,
        midpoints_profile   = converter.midpoints_profile
    ), outfile)
