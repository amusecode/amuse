from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.options import option
import os.path


class MakeMeAMassiveStarInterface(CodeInterface, CommonCodeInterface, LiteratureRefs):
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
    include_headers = ['worker_mmams.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = "worker_mmams", **options)
        LiteratureRefs.__init__(self)
    
    @option(type="string")
    def data_directory(self):
        """
        Returns the root name of the directory for the MakeMeAMassiveStar
        application data files.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mmams', 'input')
    
    @option(type="string")
    def output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mmams', 'output')
    
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
    def add_shell():
        """
        Add a new shell to an existing particle in the stellar collision code.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('cumul_mass', dtype='float64', direction=function.IN, description = "The cumulative mass from the center to the current shell of this particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of this shell")
        function.addParameter('density', dtype='float64', direction=function.IN, description = "The density of this shell")
        function.addParameter('pressure', dtype='float64', direction=function.IN, description = "The pressure of this shell")
        function.addParameter('e_thermal', dtype='float64', direction=function.IN, description = "The thermal energy of this shell")
        function.addParameter('entropy', dtype='float64', direction=function.IN, description = "The entropy of this shell")
        function.addParameter('temperature', dtype='float64', direction=function.IN, description = "The temperature of this shell")
        function.addParameter('molecular_weight', dtype='float64', direction=function.IN, description = "The molecular_weight of this shell")
        function.addParameter('H1', dtype='float64', direction=function.IN, description = "The H1 fraction of this shell")
        function.addParameter('He4', dtype='float64', direction=function.IN, description = "The He4 fraction of this shell")
        function.addParameter('O16', dtype='float64', direction=function.IN, description = "The O16 fraction of this shell")
        function.addParameter('N14', dtype='float64', direction=function.IN, description = "The N14 fraction of this shell")
        function.addParameter('C12', dtype='float64', direction=function.IN, description = "The C12 fraction of this shell")
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
    def get_number_of_shells():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('number_of_shells', dtype='int32', direction=function.OUT, description = "The number of currently defined shells of this particle")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_shell():
        """
        Return properties of the stellar model at a specific shell.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_shell', dtype='int32', direction=function.IN)
        function.addParameter('cumul_mass', dtype='float64', direction=function.OUT, description = "The cumulative mass from the center to the current shell of this particle")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The radius of this shell")
        function.addParameter('density', dtype='float64', direction=function.OUT, description = "The density of this shell")
        function.addParameter('pressure', dtype='float64', direction=function.OUT, description = "The pressure of this shell")
        function.addParameter('e_thermal', dtype='float64', direction=function.OUT, description = "The thermal energy of this shell")
        function.addParameter('entropy', dtype='float64', direction=function.OUT, description = "The entropy of this shell")
        function.addParameter('temperature', dtype='float64', direction=function.OUT, description = "The temperature of this shell")
        function.addParameter('molecular_weight', dtype='float64', direction=function.OUT, description = "The molecular_weight of this shell")
        function.addParameter('H1', dtype='float64', direction=function.OUT, description = "The H1 fraction of this shell")
        function.addParameter('He4', dtype='float64', direction=function.OUT, description = "The He4 fraction of this shell")
        function.addParameter('O16', dtype='float64', direction=function.OUT, description = "The O16 fraction of this shell")
        function.addParameter('N14', dtype='float64', direction=function.OUT, description = "The N14 fraction of this shell")
        function.addParameter('C12', dtype='float64', direction=function.OUT, description = "The C12 fraction of this shell")
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
    
    def retrieve_stellar_structure(self, star, structure = dict()):
        structure['number_of_zones']   = star.get_number_of_zones().number
        structure['number_of_species'] = star.get_number_of_species().number
        structure['species_names']     = star.get_names_of_species(number_of_species = structure['number_of_species'])
        structure['species_IDs']       = star.get_IDs_of_species(number_of_species = structure['number_of_species'])
        structure['frac_mass']   = star.get_mass_profile(number_of_zones = structure['number_of_zones'])
        structure['density']     = star.get_density_profile(number_of_zones = structure['number_of_zones'])
        structure['radius']      = star.get_radius_profile(number_of_zones = structure['number_of_zones'])
        structure['temperature'] = star.get_temperature_profile(number_of_zones = structure['number_of_zones'])
        structure['mu']          = star.get_mu_profile(number_of_zones = structure['number_of_zones'])
        structure['composition'] = star.get_chemical_abundance_profiles(
            number_of_zones = structure['number_of_zones'], number_of_species = structure['number_of_species'])
        structure['specific_internal_energy'] = (1.5 * constants.kB * structure['temperature'] / structure['mu']).as_quantity_in(units.m**2/units.s**2)
        return structure
    

class MakeMeAMassiveStar(InCodeComponentImplementation):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MakeMeAMassiveStarInterface(**options), **options)
    
    def merge_stars(self, star_1, star_2):
        structure_1 = self.retrieve_stellar_structure(star_1)
        id_1 = self.add_star_1(structure_1)
        structure_2 = self.retrieve_stellar_structure(star_2)
        id_2 = self.add_star_2(structure_2)
        id_product = self.merge_two_stars(id_1, id_2)
        product = self.convert_to_stellar_structure(id_product)
        return product
    
    def merge_stars(self, star_1, star_2):
        pass
