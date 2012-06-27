import numpy
import os.path
import pickle

from collections import namedtuple

from amuse.community.gadget2.interface import Gadget2
from amuse.ext.spherical_model import EnclosedMassInterpolator
from amuse.ext.spherical_model import new_spherical_particle_distribution
from amuse.support.exceptions import AmuseException, AmuseWarning
from amuse.units.quantities import zero
from amuse.units import units, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits

from amuse.support.console import set_printing_strategy
from amuse.datamodel import Particles, Particle

StellarModelInSPH = namedtuple('StellarModelInSPH', ['gas_particles', 'core_particle', 'core_radius'])

class StellarModel2SPH(object):
    """
    Requests the internal structure of the star from a Stellar Evolution 
    legacy code and converts it into an SPH model consisting of the 
    specified number of particles. Useful for merging stars.
    
    :argument particle: Star particle to be converted to an SPH model
    :argument number_of_sph_particles: Number of gas particles in the resulting model
    :argument with_core_particle: Model the core as a heavy, non-sph particle
    :argument target_core_mass:   If (with_core_particle): target mass for the non-sph particle
    :argument do_relax: Relax the SPH model - doesn't seem to work satisfactorily yet!
    :argument pickle_file: If provided, read stellar structure from here instead of using 'particle'
    :argument do_store_composition:  If set, store the local chemical composition on each particle
    :argument base_grid_options: dict() with options for the initial distribution, 
        see new_uniform_spherical_particle_distribution
    """
    
    def __init__(self, particle, number_of_sph_particles, seed = None, 
            do_relax = False, sph_code = Gadget2, compatible_converter = ConvertBetweenGenericAndSiUnits,
            with_core_particle = False, target_core_mass = None, pickle_file = None, 
            do_store_composition = True, gamma=5.0/3.0, base_grid_options = dict(type = "bcc")):
        self.particle = particle
        self.number_of_sph_particles = number_of_sph_particles
        
        self.with_core_particle = with_core_particle
        self.target_core_mass = target_core_mass
        self.core_radius = None
        self.core_mass = None
        
        self.pickle_file = pickle_file
        if seed:
            numpy.random.seed(seed)
        
        self.do_store_composition = do_store_composition
        self.gamma = gamma
        self.base_grid_options = base_grid_options
        self.do_relax = do_relax
        self.sph_code = sph_code # used to relax the SPH model
        self.compatible_converter = compatible_converter
    
    def retrieve_stellar_structure(self):
        self.number_of_zones   = self.particle.get_number_of_zones()
        if self.do_store_composition:
            self.number_of_species = self.particle.get_number_of_species()
            self.species_names     = self.particle.get_names_of_species(number_of_species = self.number_of_species)
            self.composition_profile = self.particle.get_chemical_abundance_profiles(
                number_of_zones = self.number_of_zones, number_of_species = self.number_of_species)
        self.density_profile     = self.particle.get_density_profile(number_of_zones = self.number_of_zones)
        self.radius_profile      = self.particle.get_radius_profile(number_of_zones = self.number_of_zones)
        temperature_profile = self.particle.get_temperature_profile(number_of_zones = self.number_of_zones)
        self.mu_profile          = self.particle.get_mu_profile(number_of_zones = self.number_of_zones)
        self.specific_internal_energy_profile = (1.5 * constants.kB * temperature_profile / self.mu_profile).as_quantity_in(units.m**2/units.s**2)
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
        self.density_profile     = structure['density_profile']
        self.radius_profile      = structure['radius_profile']
        self.mu_profile          = structure['mu_profile']
        self.composition_profile = structure['composition_profile']
        self.specific_internal_energy_profile = structure['specific_internal_energy_profile']
        self.midpoints_profile   = structure['midpoints_profile']
    
    def setup_core_parameters(self):
        self.original_entropy = (self.gamma - 1.0) * (self.specific_internal_energy_profile * 
            self.density_profile**(1.0-self.gamma))
        
        interpolator = EnclosedMassInterpolator()
        interpolator.initialize(self.radius_profile, self.density_profile)
        
        i_edge = numpy.searchsorted(interpolator.enclosed_mass, self.target_core_mass)-1
        min_i = i_edge
        max_i = len(self.radius_profile)-3
        enclosed_mass_edge = interpolator.enclosed_mass[min_i+1]
        min_enclosed_mass_residual = self.construct_model_with_core(min_i, enclosed_mass_edge, self.gamma)
        enclosed_mass_edge = interpolator.enclosed_mass[max_i+1]
        max_enclosed_mass_residual = self.construct_model_with_core(max_i, enclosed_mass_edge, self.gamma)
        
        if (min_enclosed_mass_residual > zero) or (max_enclosed_mass_residual < zero):
            raise AmuseException("Requested target_core_mass of {0} is out of range.".format(self.target_core_mass))
        
        while max_i - min_i > 1:
            next_i = (max_i + min_i)/2
            enclosed_mass_edge = interpolator.enclosed_mass[next_i+1]
            enclosed_mass_residual = self.construct_model_with_core(next_i, enclosed_mass_edge, self.gamma)
            if enclosed_mass_residual >= zero:
                max_i = next_i
            else:
                min_i = next_i
        if enclosed_mass_residual < zero:
            enclosed_mass_edge = interpolator.enclosed_mass[max_i+1]
            self.construct_model_with_core(max_i, enclosed_mass_edge, self.gamma)
        self.density_profile = self.rho
        self.specific_internal_energy_profile = self.u
        interpolator.initialize(self.radius_profile, self.density_profile)
        self.core_mass = self.mass - interpolator.enclosed_mass[-1]
        self.core_radius = self.radius_profile[max_i] / 2.8
        print "core mass in DM particle:", self.core_mass.as_quantity_in(units.MSun)
        self.mass = self.mass - self.core_mass
    
    def construct_model_with_core(self, i_edge, m_enc_edge, gamma):
        r = self.radius_profile
        rho = self.density_profile * 1
        u = self.specific_internal_energy_profile * 1
        m_enc = m_enc_edge
        r_c = r[i_edge]
        entropy = self.original_entropy * 1
#~            entropy[:i_edge+1] = (entropy[:i_edge+1] + entropy[i_edge+1]) / 2.0
        entropy[:i_edge+1] = entropy[i_edge+1]
        d_entropy = entropy[1:] - entropy[:-1]
        
        def int_Wr2_A(x):
            return (x**3 / 3.0) - (1.2 * x**5 / r_c**2) + (x**6 / r_c**3)
        
        def int_Wr2_B(x):
            return (x**3 / 1.5) - (1.5 * x**4 / r_c) + (1.2 * x**5 / r_c**2) - (x**6 / (3.0 * r_c**3))
        
        for i in range(i_edge, 0, -1):
            r_out = r[i]
            r_in = r[i-1]
            if r_out < 0.5 * r_c:
                W_int_m_enc = int_Wr2_A(r_out) - int_Wr2_A(r_in)
            elif r_in > 0.5 * r_c:
                W_int_m_enc = int_Wr2_B(r_out) - int_Wr2_B(r_in)
            else:
                W_int_m_enc = int_Wr2_B(r_out) - int_Wr2_B(r_c/2.0) + int_Wr2_A(r_c/2.0) - int_Wr2_A(r_in)
            delta_rho = (rho[i]**(2.0-gamma) * constants.G * m_enc * (r_out -r_in) / (gamma * entropy[i] * r[i]**2) +
                rho[i] * d_entropy[i] / (gamma * entropy[i]))
            
            m_enc -= 4.0 * constants.pi * rho[i] * (r_out**3 - r_in**3) / 3.0 + 32 * self.target_core_mass * W_int_m_enc / r_c**3
            if m_enc < zero:
                break
            rho[i-1] = rho[i] + delta_rho
            u[i-1] = entropy[i-1] * rho[i-1]**(gamma-1.0) / (gamma-1.0)
        
        self.rho = rho
        self.u = u
        return m_enc
            
    
    def get_index(self, value, sorted_vector):
        if not sorted_vector[0] <= value <= sorted_vector[-1]:
            raise AmuseException("Can't find a valid index. {0} is not in "
                "the range [{1}, {2}].".format(value, sorted_vector[0], sorted_vector[-1]))
        index = numpy.searchsorted(sorted_vector, value)
        return max(index - 1, 0)
    
    def get_indices(self,values,sorted_vector):
        if values.amin() < sorted_vector[0] or values.amax()>sorted_vector[-1]:
            raise AmuseException("Can't find a valid index.  not in "
                "the range [{0}, {1}].".format(sorted_vector[0], sorted_vector[-1]))
        indices=numpy.maximum(numpy.searchsorted(sorted_vector,values)-1,0)
        return indices

    
    def calculate_interpolation_coefficients(self, radial_positions):
#        indices = numpy.array([self.get_index(r, self.midpoints_profile) for r in radial_positions])
        indices=self.get_indices(radial_positions,self.midpoints_profile)
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
            comp = [] 
            for species in self.composition_profile:
                extended = list(species[:1])
                extended.extend(species)
                extended.append(species[-1])
                extended = numpy.asarray(extended)
                comp.append(delta*extended[indices] + one_minus_delta*extended[indices+1])
            
            comp = numpy.asarray(comp)
            
            extended = self.mu_profile[:1]
            extended.extend(self.mu_profile)
            extended.append(self.mu_profile[-1])
            mu = delta*extended[indices] + one_minus_delta*extended[indices+1]
            
            return interpolated_energies, comp.transpose(), mu
        else:
            return interpolated_energies, None, None
    
    def convert_to_SPH(self):
        sph_particles = new_spherical_particle_distribution(
            self.number_of_sph_particles, 
            radii = self.radius_profile, densities = self.density_profile, 
            **self.base_grid_options
        )
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
        hydro_code = self.sph_code(unit_converter)
        particles.u = 1.0 | (units.m / units.s)**2
        hydro_code.gas_particles.add_particles(particles)
        
        for i in range(1, num_iterations+1):
            hydro_code.gas_particles.u, tmp, tmp2 = self.interpolate_internal_energy(particles.position.lengths(), do_composition_too = False)
            hydro_code.evolve_model(i * (1.0e-5 | units.s))
            accelerations     = hydro_code.gas_particles.acceleration
            acc_correlated = (previous_acc * accelerations).sum() / (accelerations * accelerations).sum()
            if (acc_correlated < 0.5 | units.none and i > 2):
                break
            
            previous_acc = accelerations
            internal_energies = hydro_code.gas_particles.u
            smoothing_lengths = hydro_code.gas_particles.h_smooth
            factor = numpy.minimum((max_delta * internal_energies / (accelerations.lengths() * smoothing_lengths)), 0.5)
            result.append(str(i) + ": Accelerations correlated: " + str(acc_correlated) + ", median factor: " + str(numpy.median(factor)))
            
            particles.position += accelerations * ((smoothing_lengths * smoothing_lengths * factor) / 
                internal_energies).reshape((self.number_of_sph_particles, 1))
            hydro_code.gas_particles.position = particles.position
            hydro_code.gas_particles.velocity = particles.velocity
        
        particles.u = hydro_code.gas_particles.u
        hydro_code.stop()
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
        
        if self.with_core_particle:
            self.setup_core_parameters()
        
        sph_particles = self.convert_to_SPH()
        if self.do_relax:
            for result_string in self.relax(sph_particles):
                print result_string
        
        specific_internal_energy, composition, mu = self.interpolate_internal_energy(
            sph_particles.position.lengths(),
            do_composition_too = self.do_store_composition
        )
        sph_particles.u = specific_internal_energy
        if self.do_store_composition:
            sph_particles.add_vector_attribute("composition", self.species_names)
            sph_particles.composition = composition
            sph_particles.mu = mu
        
        if self.with_core_particle and self.core_radius:
            core_particle = Particle()
            core_particle.mass = self.core_mass
            core_particle.position = [0.0, 0.0, 0.0] | units.m
            core_particle.velocity = [0.0, 0.0, 0.0] | units.m / units.s
            core_particle.radius = self.core_radius
            return StellarModelInSPH(gas_particles=sph_particles, core_particle=core_particle, core_radius=self.core_radius)
        return StellarModelInSPH(gas_particles=sph_particles, core_particle=Particle(), core_radius=None)

    

def convert_stellar_model_to_SPH(particle, number_of_sph_particles, **keyword_arguments):
    """
    Requests the internal structure of the star from a Stellar Evolution 
    legacy code and converts it into an SPH model consisting of the 
    specified number of particles. Useful for merging stars.
    
    :argument particle: Star particle to be converted to an SPH model
    :argument number_of_sph_particles: Number of gas particles in the resulting model
    :argument with_core_particle: Model the core as a heavy, non-sph particle
    :argument target_core_mass:   If (with_core_particle): target mass for the non-sph particle
    :argument do_relax: Relax the SPH model - doesn't seem to work satisfactorily yet!
    :argument pickle_file: If provided, read stellar structure from here instead of using 'particle'
    :argument do_store_composition:  If set, store the local chemical composition on each particle
    :argument base_grid_options: dict() with options for the initial distribution, 
        see new_uniform_spherical_particle_distribution
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
        density_profile     = converter.density_profile,
        radius_profile      = converter.radius_profile,
        mu_profile          = converter.mu_profile,
        composition_profile = converter.composition_profile,
        specific_internal_energy_profile = converter.specific_internal_energy_profile,
        midpoints_profile   = converter.midpoints_profile
    ), outfile)
