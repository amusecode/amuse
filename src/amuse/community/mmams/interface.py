import numpy
from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.options import option
from amuse.units import units
from amuse.datamodel import Particle

import os.path

class MMAMSInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn,
        CodeWithDataDirectories):
    """
    MakeMeAMassiveStar is a computationally inexpensive method in which the 
    merger process is approximated, including shock heating, hydrodynamic 
    mixing and mass loss, with a simple algorithm based on conservation laws and 
    a basic qualitative understanding of the hydrodynamics of stellar mergers. 
    The algorithm relies on Archimedes' principle to dictate the distribution of 
    the fluid in the stable equilibrium situation. Without the effects of 
    microscopic mixing, the temperature and chemical composition profiles in a 
    collision product can become double-valued functions of enclosed mass. Such 
    an unphysical situation is mended by simulating microscopic mixing as a 
    post-collision effect.
    
    Relevant references:
        .. [#] Gaburov E., Lombardi J. C. & Portegies Zwart S., 2008, MNRAS, 383, L5
    """
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = "mmams_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar collision code. The particle is 
        initialized as an empty model (with zero shells). The model has to be 
        constructed shell by shell using `add_shell'. An index is returned that 
        can be used to refer to this particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            new empty particle was created
        -1 - ERROR
            particle could not be created"""
        return function
    
    @legacy_function
    def delete_particle():
        """
        Remove a particle from the stellar collision code.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def add_shell():
        """
        Add a new shell to an existing particle in the stellar collision code.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('d_mass', dtype='float64', direction=function.IN, description = "The mass of the current shell of this particle")
        function.addParameter('cumul_mass', dtype='float64', direction=function.IN, description = "The cumulative mass from the center to the current shell of this particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of this shell")
        function.addParameter('density', dtype='float64', direction=function.IN, description = "The density of this shell")
        function.addParameter('pressure', dtype='float64', direction=function.IN, description = "The pressure of this shell")
        function.addParameter('temperature', dtype='float64', direction=function.IN, description = "The temperature of this shell")
        function.addParameter('luminosity', dtype='float64', direction=function.IN, description = "The luminosity of this shell")
        function.addParameter('molecular_weight', dtype='float64', direction=function.IN, description = "The molecular_weight of this shell")
        function.addParameter('H1', dtype='float64', direction=function.IN, description = "The H1 fraction of this shell")
        function.addParameter('He4', dtype='float64', direction=function.IN, description = "The He4 fraction of this shell")
        function.addParameter('C12', dtype='float64', direction=function.IN, description = "The C12 fraction of this shell")
        function.addParameter('N14', dtype='float64', direction=function.IN, description = "The N14 fraction of this shell")
        function.addParameter('O16', dtype='float64', direction=function.IN, description = "The O16 fraction of this shell")
        function.addParameter('Ne20', dtype='float64', direction=function.IN, description = "The Ne20 fraction of this shell")
        function.addParameter('Mg24', dtype='float64', direction=function.IN, description = "The Mg24 fraction of this shell")
        function.addParameter('Si28', dtype='float64', direction=function.IN, description = "The Si28 fraction of this shell")
        function.addParameter('Fe56', dtype='float64', direction=function.IN, description = "The Fe56 fraction of this shell")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            new shell was added to the specified particle
        -1 - ERROR
            specified particle was not found"""
        return function
    
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT, description = "The number of currently defined particles")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_number_of_zones():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('number_of_zones', dtype='int32', direction=function.OUT, description = "The number of currently defined shells of this particle")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', unit=units.MSun, direction=function.OUT, description = "The total mass of this particle")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_stellar_model_element():
        """
        Return properties of the stellar model at a specific zone.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_zone', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('d_mass', dtype='float64', direction=function.OUT, description = "The mass of the current shell of this particle")
        function.addParameter('cumul_mass', dtype='float64', direction=function.OUT, description = "The cumulative mass from the center to the current shell of this particle")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The radius of this shell")
        function.addParameter('density', dtype='float64', direction=function.OUT, description = "The density of this shell")
        function.addParameter('pressure', dtype='float64', direction=function.OUT, description = "The pressure of this shell")
        function.addParameter('entropy', dtype='float64', direction=function.OUT, description = "The entropy of this shell")
        function.addParameter('temperature', dtype='float64', direction=function.OUT, description = "The temperature of this shell")
        function.addParameter('luminosity', dtype='float64', direction=function.OUT, description = "The luminosity of this shell")
        function.addParameter('molecular_weight', dtype='float64', direction=function.OUT, description = "The molecular_weight of this shell")
        function.addParameter('H1', dtype='float64', direction=function.OUT, description = "The H1 fraction of this shell")
        function.addParameter('He4', dtype='float64', direction=function.OUT, description = "The He4 fraction of this shell")
        function.addParameter('C12', dtype='float64', direction=function.OUT, description = "The C12 fraction of this shell")
        function.addParameter('N14', dtype='float64', direction=function.OUT, description = "The N14 fraction of this shell")
        function.addParameter('O16', dtype='float64', direction=function.OUT, description = "The O16 fraction of this shell")
        function.addParameter('Ne20', dtype='float64', direction=function.OUT, description = "The Ne20 fraction of this shell")
        function.addParameter('Mg24', dtype='float64', direction=function.OUT, description = "The Mg24 fraction of this shell")
        function.addParameter('Si28', dtype='float64', direction=function.OUT, description = "The Si28 fraction of this shell")
        function.addParameter('Fe56', dtype='float64', direction=function.OUT, description = "The Fe56 fraction of this shell")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def read_usm():
        """
        Define a new particle in the stellar collision code. Read the stellar model from a usm file.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('usm_file', dtype='string', direction=function.IN,
            description = "The path to the usm file.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def merge_two_stars():
        """
        Merge two previously defined particles. Returns an index to the new stellar model.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_product', dtype='int32', direction=function.OUT)
        function.addParameter('index_of_the_primary', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_secondary', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_dump_mixed_flag():
        """Set the dump_mixed flag: specifies whether the returned products must be mixed first."""
        function = LegacyFunctionSpecification()
        function.addParameter('dump_mixed_flag', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_dump_mixed_flag():
        """Get the dump_mixed flag: specifies whether the returned products must be mixed first."""
        function = LegacyFunctionSpecification()
        function.addParameter('dump_mixed_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_target_n_shells_mixing():
        """Set the target number of shells for mixed models."""
        function = LegacyFunctionSpecification()
        function.addParameter('target_n_shells_mixing', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_target_n_shells_mixing():
        """Get the target number of shells for mixed models."""
        function = LegacyFunctionSpecification()
        function.addParameter('target_n_shells_mixing', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_target_n_shells():
        """Set the target number of shells for unmixed models."""
        function = LegacyFunctionSpecification()
        function.addParameter('target_n_shells', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_target_n_shells():
        """Get the target number of shells for unmixed models."""
        function = LegacyFunctionSpecification()
        function.addParameter('target_n_shells', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_do_shock_heating_flag():
        """Set the flag that specifies whether to solve for shock-heating."""
        function = LegacyFunctionSpecification()
        function.addParameter('do_shock_heating_flag', dtype='int32', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_do_shock_heating_flag():
        """Get the flag that specifies whether to solve for shock-heating."""
        function = LegacyFunctionSpecification()
        function.addParameter('do_shock_heating_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function


class MMAMS(CommonCode):
    
    stellar_evolution_code_required = True
    gravity_code_required = False
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MakeMeAMassiveStarInterface(**options), **options)
        self.create_new_key = True
    
    def get_create_new_key(self):
        return self.create_new_key
    
    def set_create_new_key(self, value):
        self.create_new_key = value
    
    def define_properties(self, handler):
        handler.add_property("get_number_of_particles")
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_target_n_shells_mixing", 
            "set_target_n_shells_mixing",
            "target_n_shells_mixing", 
            "target number of shells for mixed models", 
            default_value = 200
        )
        
        handler.add_method_parameter(
            "get_target_n_shells", 
            "set_target_n_shells",
            "target_n_shells", 
            "target number of shells for unmixed models", 
            default_value = 10000
        ) 
        
        handler.add_boolean_parameter(
            "get_dump_mixed_flag",
            "set_dump_mixed_flag",
            "dump_mixed_flag",
            "dump_mixed flag: specifies whether the returned products must be mixed first",
            True
        )
        
        handler.add_boolean_parameter(
            "get_do_shock_heating_flag",
            "set_do_shock_heating_flag",
            "do_shock_heating_flag",
            "do_shock_heating flag: specifies whether to solve for shock-heating",
            True
        )
        
        handler.add_method_parameter(
            "get_create_new_key",
            "set_create_new_key",
            "create_new_key",
            "Flag to specify whether the merge product gets a new key, or the key of the most massive collider.",
            True
        )
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (units.MSun,),
            (handler.INDEX, handler.ERROR_CODE,)
        )
        handler.add_method(
            "delete_particle",
            (handler.INDEX,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "read_usm",
            (handler.NO_UNIT,),
            (handler.INDEX, handler.ERROR_CODE,)
        )
        handler.add_method(
            "add_shell",
            (handler.INDEX, units.MSun, units.MSun, units.RSun, units.g / units.cm**3, units.barye, 
                units.K, units.LSun, units.amu, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_stellar_model_element",
            (handler.INDEX, handler.INDEX,),
            (units.MSun, units.MSun, units.RSun, units.g / units.cm**3, units.barye, 
                handler.NO_UNIT, units.K, units.LSun, units.amu,
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_number_of_zones",
            (handler.INDEX, ),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_number_of_particles",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "merge_two_stars",
            (handler.INDEX, handler.INDEX,),
            (handler.INDEX, handler.ERROR_CODE,)
        )
        
        handler.add_method("get_target_n_shells_mixing", (), (handler.NO_UNIT, handler.ERROR_CODE,)) 
        handler.add_method("set_target_n_shells_mixing", (handler.NO_UNIT, ), (handler.ERROR_CODE,)) 
        
        handler.add_method("get_target_n_shells", (), (handler.NO_UNIT, handler.ERROR_CODE,)) 
        handler.add_method("set_target_n_shells", (handler.NO_UNIT, ), (handler.ERROR_CODE,)) 
        
        handler.add_method("get_number_of_particles", (), (handler.NO_UNIT, handler.ERROR_CODE,)) 
    
    def define_errorcodes(self, handler):
        handler.add_errorcode(-1, 'Unspecified, other error.')
        handler.add_errorcode(-2, 'Specified zone is undefined for this particle.')
        handler.add_errorcode(-3, 'A particle with the given index was not found.')
    
    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['native_stars', 'imported_stars', 'merge_products'], 
            index_to_default_set = 0)
        
        handler.define_set('native_stars', 'index_of_the_particle')
        handler.set_new('native_stars', 'new_particle')
        
        handler.define_set('imported_stars', 'index_of_the_particle')
        handler.set_new('imported_stars', 'read_usm')
        
        handler.define_set('merge_products', 'index_of_the_particle')
        handler.set_new('merge_products', 'merge_stars')
        
        for particle_set_name in ['native_stars', 'imported_stars', 'merge_products']:
            handler.set_delete(particle_set_name, 'delete_particle')
            handler.add_getter(particle_set_name, 'get_number_of_zones')
            handler.add_getter(particle_set_name, 'get_mass')
            handler.add_method(particle_set_name, 'add_shell') 
            handler.add_method(particle_set_name, 'get_stellar_model', 'get_internal_structure') 
    
    
    def get_stellar_model(self, index_of_the_particle):
        if hasattr(index_of_the_particle, '__iter__'):
            return [self._create_new_grid(self._specify_stellar_model, index_of_the_particle = x) for x in index_of_the_particle]
        else:
            return self._create_new_grid(self._specify_stellar_model, index_of_the_particle = index_of_the_particle)
    
    def get_range_in_zones(self, index_of_the_particle):
        """
        Returns the inclusive range of defined zones/mesh-cells of the star.
        """
        return (0, self.get_number_of_zones(index_of_the_particle)-1)
    
    def _specify_stellar_model(self, definition, index_of_the_particle = 0):
        definition.set_grid_range('get_range_in_zones')
        definition.add_getter('get_stellar_model_element', names=('d_mass', 'mass', 'radius', 
            'rho', 'pressure', 'entropy', 'temperature', 'luminosity', 'molecular_weight', 
            'X_H', 'X_He', 'X_C', 'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe'))
        definition.define_extra_keywords({'index_of_the_particle':index_of_the_particle})
    
    def _key_to_index_in_code(self, key):
        if self.native_stars.has_key_in_store(key):
            return self.native_stars._subset([key]).index_in_code[0]
        elif self.imported_stars.has_key_in_store(key):
            return self.imported_stars._subset([key]).index_in_code[0]
        else:
            return self.merge_products._subset([key]).index_in_code[0]
    
    def merge_stars(self, primary, secondary):
        indices_of_primaries = [self._key_to_index_in_code(particle.key) for particle in primary]
        indices_of_secondaries = [self._key_to_index_in_code(particle.key) for particle in secondary]
        result = self.merge_two_stars(
            indices_of_primaries, 
            indices_of_secondaries
        )
        return result
    
    def match_composition_to_mass_profile(self, stellar_model, mass_profile):
        new_composition = numpy.array([[0.0]*len(mass_profile)]*8)
        current_index = 0
        previous_mass = 0.0 | units.MSun
        for i, mass_i in enumerate(mass_profile):
            previous_index = current_index
            mass_fractions = [] | units.MSun
            while stellar_model.mass[current_index] < mass_i:
                mass_fractions.append(stellar_model.mass[current_index]-previous_mass)
                previous_mass = stellar_model.mass[current_index]
                current_index += 1
            if stellar_model.mass[current_index] > mass_i:
                mass_fractions.append(mass_i-previous_mass)
                previous_mass = mass_i
            weights = mass_fractions / mass_fractions.sum()
            next_index = previous_index+len(weights)
            new_composition[0, i] = (stellar_model.X_H[previous_index:next_index]*weights).sum()
            new_composition[2, i] = (stellar_model.X_He[previous_index:next_index]*weights).sum()
            new_composition[3, i] = (stellar_model.X_C[previous_index:next_index]*weights).sum()
            new_composition[4, i] = (stellar_model.X_N[previous_index:next_index]*weights).sum()
            new_composition[5, i] = (stellar_model.X_O[previous_index:next_index]*weights).sum()
            new_composition[6, i] = (stellar_model.X_Ne[previous_index:next_index]*weights).sum()
            new_composition[7, i] = ((stellar_model.X_Mg[previous_index:next_index] +
                stellar_model.X_Si[previous_index:next_index] + 
                stellar_model.X_Fe[previous_index:next_index])*weights).sum()
        return new_composition
    
    def key_for_merge_product(self, primary, secondary):
        if self.create_new_key:
            return None
        else:
            return primary.key if primary.mass > secondary.mass else secondary.key
    
    def handle_collision(self, primary, secondary, stellar_evolution_code=None, gravity_code=None):
        primary = self.add_particle_with_internal_structure(primary, stellar_evolution_code=stellar_evolution_code)
        secondary = self.add_particle_with_internal_structure(secondary, stellar_evolution_code=stellar_evolution_code)
        merge_product = Particle(
            key = self.key_for_merge_product(primary, secondary),
            primary = primary,
            secondary = secondary
        )
        return self.merge_products.add_particles(merge_product.as_set())
    
    def add_particle_with_internal_structure(self, particle, stellar_evolution_code=None):
        new = self.native_stars.add_particle(particle)
        if not stellar_evolution_code is None:
            particle = particle.as_set().get_intersecting_subset_in(stellar_evolution_code.particles)[0]
        number_of_zones     = particle.get_number_of_zones()
        density_profile     = particle.get_density_profile(number_of_zones = number_of_zones)
        radius_profile      = particle.get_radius_profile(number_of_zones = number_of_zones)
        if hasattr(particle, "get_mass_profile"):
            mass_profile        = particle.get_mass_profile(number_of_zones = number_of_zones)* particle.mass
        else:
            print "mass profile from density and radius"
            radii_cubed = radius_profile**3
            radii_cubed.prepend(0|units.m**3)
            mass_profile = (4.0/3.0 * numpy.pi) * density_profile * (radii_cubed[1:] - radii_cubed[:-1])
        if hasattr(particle, "get_cumulative_mass_profile"):
            cumul_mass_profile  = particle.get_cumulative_mass_profile(number_of_zones = number_of_zones) * particle.mass
        else:
            print "cumulative mass profile from mass profile"
            cumul_mass_profile  = mass_profile.accumulate()
        
        temperature_profile = particle.get_temperature_profile(number_of_zones = number_of_zones)
        #lum                 = particle.get_luminosity_profile(number_of_zones = number_of_zones)
        lum                 = numpy.zeros_like(temperature_profile) | units.LSun
        pressure_profile    = particle.get_pressure_profile(number_of_zones = number_of_zones)
        mu_profile          = particle.get_mu_profile(number_of_zones = number_of_zones)
        composition_profile = particle.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
        new.add_shell(mass_profile, cumul_mass_profile, radius_profile, density_profile, 
            pressure_profile, temperature_profile, lum, mu_profile, composition_profile[0], 
            composition_profile[1]+composition_profile[2], composition_profile[3], 
            composition_profile[4], composition_profile[5], composition_profile[6], 
            composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
        return new

MakeMeAMassiveStarInterface = MMAMSInterface
MakeMeAMassiveStar = MMAMS
