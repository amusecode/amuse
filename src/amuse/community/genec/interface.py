from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function
from amuse.community import LiteratureReferencesMixIn
from amuse.community import CodeWithDataDirectories
from amuse.community import InCodeComponentImplementation
from amuse.community.interface.se import StellarEvolution
from amuse.community.interface.se import StellarEvolutionInterface
from amuse.community.interface.se import InternalStellarStructure
from amuse.community.interface.se import InternalStellarStructureInterface
from amuse.units import units


class GenecInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    StellarEvolutionInterface,
    InternalStellarStructureInterface,
    CodeWithDataDirectories,
):
    """
    GENEC is the Geneva Stellar Evolution Code

    References:
        .. [#] The Geneva Stellar Evolution Group
    """

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="genec_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)

    @legacy_function
    def set_genec_path():
        """
        set path to Genec data directories
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'genec_path', dtype='string', direction=function.IN,
            description="Path to Genec",
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function

    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_starname():
        """
        Set the star name (identical to AMUSE particle key?)
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'index_of_the_star', dtype='int32',
            direction=function.IN,
            description="The star's key"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            Unable to set.
        -2 - ERROR
            Cannot set at this point, already running.
        """
        return function

    @legacy_function
    def new_particle():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.OUT,
            description=(
                "The new index for the star. This index can be used to refer "
                "to this star in other functions"
            )
        )
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The initial mass of the star")
        function.addParameter(
            'metallicity', dtype='float64', direction=function.IN,
            description="The initial metallicity of the star")
        function.addParameter(
            'key', dtype='float64', direction=function.IN,
            description="The identifying key of the star")
        # function.addParameter(
        #     'age_tag', dtype='float64', direction=function.IN,
        #     description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function

    @legacy_function
    def get_luminosity_at_zone():
        """
        Retrieve the luminosity at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of"
        )
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of"
        )
        function.addParameter(
            'lum_i', dtype='float64', direction=function.OUT,
            description=(
                "The luminosity at the specified zone/mesh-cell of the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function

    @legacy_function
    def get_mass_fraction_at_zone():
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_star', dtype='int32', direction=function.IN,
            description="The index of the star to get the value of")
        function.addParameter(
            'zone', dtype='int32', direction=function.IN,
            description="The zone/mesh-cell of the star to get the value of")
        function.addParameter(
            'dq_i', dtype='float64', direction=function.OUT,
            description=(
                "The mass fraction at the specified zone/mesh-cell of "
                "the star."
            )
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function


class Genec(StellarEvolution, InternalStellarStructure):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(
            self,  GenecInterface(**options), **options)
        self.model_time = 0.0 | units.yr

    def define_parameters(self, handler):
        handler.add_method_parameter(
            None,
            "set_genec_path",
            "path_to_data",
            "Path to the data directory",
            default_value=self.data_directory
        )

    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_star')
        handler.set_new('particles', 'new_particle')

        handler.add_getter('particles', 'get_radius')
        handler.add_getter('particles', 'get_mass')
        handler.add_getter('particles', 'get_age')
        handler.add_getter('particles', 'get_luminosity')
        handler.add_getter('particles', 'get_temperature')
        # handler.add_method('particles', 'get_number_of_zones')
        # handler.add_method('particles', 'get_radius_profile')
        # handler.add_method('particles', 'get_temperature_profile')
        handler.add_method('particles', 'get_luminosity_profile')
        handler.add_method('particles', 'get_cumulative_mass_profile')

        handler.add_method('particles', 'evolve_one_step')
        handler.add_method('particles', 'evolve_for')
        InternalStellarStructure.define_particle_sets(
            self, handler, set_name='particles'
        )

        handler.set_delete('particles', 'delete_star')

    def define_state(self, handler):
        StellarEvolution.define_state(self, handler)
        # InternalStellarStructure.define_state(self, handler)

        # Only allow setting of starname in EDIT or UPDATE states
        # I.e. must do initialize_code and commit_parameters FIRST!

        # Initialized (initialize_code)

        # -> Edit (commit_parameters)
        handler.add_method('EDIT', 'set_starname')
        # handler.add_method('EDIT', 'new_particle')

        # -> Run (commit_particles)

        # -> Update

        # -> Run (recommit_particles)

        handler.add_method('UPDATE', 'set_starname')
        # handler.add_method('UPDATE', 'new_particle')

    def define_methods(self, handler):
        InternalStellarStructure.define_methods(self, handler)
        StellarEvolution.define_methods(self, handler)
        handler.add_method(
            "new_particle",
            (units.MSun, handler.NO_UNIT, handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "get_radius",
            (handler.INDEX,),
            (units.cm, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_number_of_zones",
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_radius_at_zone",
            (handler.INDEX, handler.NO_UNIT,),
            (units.cm, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_temperature_at_zone",
            (handler.INDEX, handler.NO_UNIT,),
            (units.K, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_luminosity_at_zone",
            (handler.INDEX, handler.NO_UNIT,),
            (units.erg/units.s, handler.ERROR_CODE,)
        )
    # def define_parameters(self, handler):

    def get_luminosity_profile(
            self,
            indices_of_the_stars,
            number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying luminosity profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_luminosity_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_mass_profile(
        self,
        indices_of_the_stars,
        number_of_zones=None,
    ):
        indices_of_the_stars = self._check_number_of_indices(
            indices_of_the_stars,
            action_string="Querying mass profiles"
        )
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_mass_fraction_at_zone(
            [indices_of_the_stars]*number_of_zones,
            list(range(number_of_zones)) | units.none
        )

    def get_cumulative_mass_profile(
            self, indices_of_the_stars, number_of_zones=None):
        frac_profile = self.get_mass_profile(
            indices_of_the_stars, number_of_zones=number_of_zones)
        return frac_profile.cumsum()
