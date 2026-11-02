from amuse.community.interface.chem import (
    ChemNumberDensityTemperature,
    ChemicalEvolution,
    ChemicalEvolutionInterface,
    ChemNumberDensityTemperatureInterface
)
from amuse.rfi.core import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.units import units as u


# (Grevesse & Sauval, 1998, Space Sci. Rev. 85, 161)
solar_abundances= dict(
    H = 1.0, HE = 0.085,
    C = 3.31e-4, N = 8.3e-5, O = 6.76e-4,
    Ne = 1.2e-4, SI = 3.55e-5, Fe = 3.2e-5
)


class KromeInterface(
    CodeInterface,
    ChemicalEvolutionInterface,
    ChemNumberDensityTemperatureInterface,
    LiteratureReferencesMixIn
):
    """
    KROME: a package to embed chemistry in astrophysical simulations.

    .. [#] Grassi, T.; Bovino, S.; Schleicher, D. R. G.; Prieto, J.; Seifried, D.; Simoncini, E.; Gianturco, F. A., MNRAS, 439, 3, p.2386-2419 [2014MNRAS.439.2386G]
    """

    def __init__(self, **options):
        CodeInterface.__init__(
            self, name_of_the_worker = self.name_of_the_worker(), **options
        )
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'krome_worker'

    @legacy_function
    def new_particle():
        """
        Add a new particle to Krome.

        Parameters
        ----------
        number_density : float
            Number density of the particle in units of cm**-3.
        temperature : float
            Temperature of the particle in units of K.
        ionrate : float
            Ionization rate of the particle in units of s**-1.

        Returns
        -------
        index_of_the_particle : int
            Id of the newly created particle, returned by Krome.
        int :
            0 on success.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.OUT)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_state():
        """
        Retrieve the state of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.

        Returns
        -------
        number_density : float
            Number density of the particle in units of cm**-3.
        temperature : float
            Temperature of the particle in units of K.
        ionrate : float
            Ionization rate of the particle in units of s**-1.
        int :
            0 on success.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_state():
        """
        Set the state of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle as returned by `new_particle`.
        number_density : float
            Number density of the particle in units of cm**-3.
        temperature : float
            Temperature of the particle in units of K.
        ionrate : float
            Ionization rate of the particle in units of s**-1.

        Returns
        -------
        int :
            0 on success.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for x in ['number_density','temperature','ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function


class Krome(ChemicalEvolution, ChemNumberDensityTemperature):
    """Krome is a package to embed chemistry in astrophysical simulations."""
    def __init__(self, unit_converter=None, **options):

        if unit_converter is not None:
            raise Exception('Krome uses predefined units')

        ChemicalEvolution.__init__(self, KromeInterface(**options))

        first, last = self.get_firstlast_species_index()
        self.species = dict()
        for i in range(first, last+1):
          self.species[self.get_species_name(i)] = i

    def define_properties(self, handler):
        handler.add_property('get_time', public_name='model_time')

    def define_methods(self, handler):
        ChemicalEvolution.define_methods(self, handler)
        ChemNumberDensityTemperature.define_methods(self, handler)
        handler.add_method(
            'evolve_model', (u.s,), (handler.ERROR_CODE,)
        )

        handler.add_method(
            'new_particle',
            (u.cm**-3, u.K, u.s**-1,),
            (handler.INDEX, handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_state',
            (handler.INDEX,),
            (u.cm**-3, u.K, u.s**-1, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_state',
            (handler.INDEX, u.cm**-3, u.K, u.s**-1,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_time', (), (u.s, handler.ERROR_CODE,)
        )

    def define_particle_sets(self, handler):
        ChemicalEvolution.define_particle_sets(self, handler)
        ChemNumberDensityTemperature.define_particle_sets(self, handler)

    def define_state(self, handler):
        ChemicalEvolution.define_state(self, handler)
        handler.add_method('EDIT', 'set_abundance')
        handler.add_method('EDIT', 'set_abundances')
        handler.add_method('EDIT', 'get_abundance')
