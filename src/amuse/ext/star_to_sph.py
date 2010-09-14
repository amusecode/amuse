import numpy
from amuse.support.data.core import Particles
from amuse.support.units import units, constants
from amuse.support.exceptions import AmuseWarning, AmuseException
from amuse.ext.spherical_model import new_spherical_particle_distribution
from amuse.legacy.gadget2.interface import Gadget2
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits


class StellarModel2SPH(object):
    def __init__(self, particle, number_of_sph_particles, seed = None, mode = "scaling method", 
            do_relax = False, sph_legacy_code = Gadget2, compatible_converter = ConvertBetweenGenericAndSiUnits):
        self.particle = particle
        self.number_of_sph_particles = number_of_sph_particles
        self.zone_index = None
        self.delta      = None
        numpy.random.seed(seed)
        self.mode = mode  # "random sampling" or "scaling method"
        
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
    
    def get_index(self, radius, midpoints):
        index = 1
        while radius < midpoints[index]: index += 1
        return index - 1
    
    def interpolate_hydro_quantities(self, radial_positions):
        if self.mode == "scaling method":
            # Note: self.radius is in decreasing order; from surface to center
            radius_profile = self.radius # outer radius of each mesh zone
            radius_profile.append(0|units.m)
            midpoints = (3*radius_profile[0:1] - radius_profile[1:2])/2    # dummy element to handle boundaries correctly
            midpoints.extend((radius_profile[1:] + radius_profile[:-1])/2) # real midpoints of each mesh zone
            midpoints.append(-midpoints[-1])                               # dummy element to handle boundaries correctly
            indices = numpy.array([self.get_index(r, midpoints) for r in radial_positions])
            delta2 = (midpoints[indices+1] - radial_positions) / (midpoints[indices+1] - midpoints[indices])
        elif self.mode == "random sampling":
            delta2 = self.delta - 0.5
            indices = self.zone_index + 1
            neg = numpy.where(delta2 < (0|units.none))[0]
            delta2[neg] = delta2[neg] + 1
            indices[neg] = indices[neg] + 1
        else:
            raise AmuseException("Unknown mode: {0}. Mode can be 'scaling method' or 'random sampling'.".format(self.mode))
        one_minus_delta2 = 1 - delta2
        
        hydro_quantities = [self.specific_internal_energy] #, self.temperature, self.luminosity, self.density]
        result = []
        for hydro_quantity in hydro_quantities:
            extended = hydro_quantity[:1]
            extended.extend(hydro_quantity)
            extended.append(hydro_quantity[-1])
            result.append(delta2*extended[indices] + one_minus_delta2*extended[indices+1])
        
        comp = [] | units.none
        for hydro_quantity in self.composition:
            extended = hydro_quantity[:1]
            extended.extend(hydro_quantity)
            extended.append(hydro_quantity[-1])
            comp.append(delta2*extended[indices] + one_minus_delta2*extended[indices+1])
        return result, comp.transpose()
    
    def convert_to_SPH(self):
        if self.mode == "scaling method":
            sph_particles = new_spherical_particle_distribution(
                self.number_of_sph_particles, 
                radii = self.radius, densities = self.density)
        elif self.mode == "random sampling":
            sph_particles = Particles(self.number_of_sph_particles)
            sph_particles.position = self.new_positions()
            sph_particles.rands = units.none.new_quantity(self.rands)
            sph_particles.int = self.delta + self.zone_index
        else:
            raise AmuseException("Unknown mode: {0}. Mode can be 'scaling method' or 'random sampling'.".format(self.mode))
        sph_particles.mass = (self.particle.mass.number * 1.0 / 
            self.number_of_sph_particles) | self.particle.mass.unit
        sph_particles.velocity = [0,0,0] | units.m/units.s
        (specific_internal_energy,), composition = self.interpolate_hydro_quantities(sph_particles.position.lengths())
        sph_particles.u = specific_internal_energy
        sph_particles.add_vector_attribute("composition", self.species_names)
        sph_particles.composition = composition
        return sph_particles
    
    def relax(self, particles):
        num_iterations = 40
        max_delta = 0.005  # maximum change to particle positions relative to its smoothing length
        
        result = []
        previous_acc = 0 | units.m / units.s**2
        unit_converter = self.compatible_converter(self.particle.radius, self.particle.mass, 1.0e-3 | units.s)
        hydro_legacy_code = self.sph_legacy_code(unit_converter)
        hydro_legacy_code.gas_particles.add_particles(particles)
        
        for i in range(1, num_iterations+1):
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
        return sph_particles
    

"""
Requests the internal structure of the star from a Stellar Evolution 
legacy code and converts it into an SPH model consisting of the 
specified number of particles. Useful for merging stars.

:argument particle: Star particle to be converted to an SPH model
:argument number_of_sph_particles: Number of gas particles in the resulting model
"""
def convert_stellar_model_to_SPH(particle, number_of_sph_particles, **keyword_arguments):
    converter = StellarModel2SPH(particle, number_of_sph_particles, **keyword_arguments)
    return converter.result
