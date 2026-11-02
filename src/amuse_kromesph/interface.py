from amuse.community.interface.chem import (
    ChemDensityInternalEnergy,
    ChemDensityInternalEnergyInterface,
    ChemicalEvolution,
    ChemicalEvolutionInterface,
)
from amuse.rfi.core import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
    remote_function
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.units import units as u


# (Grevesse & Sauval, 1998, Space Sci. Rev. 85, 161)
solar_abundances= dict(
    H = 1.0, HE = 0.085,
    C = 3.31e-4, N = 8.3e-5, O = 6.76e-4,
    Ne = 1.2e-4, SI = 3.55e-5, Fe = 3.2e-5
)


class KromeSphInterface(
    CodeInterface,
    ChemicalEvolutionInterface,
    ChemDensityInternalEnergyInterface,
    LiteratureReferencesMixIn
):
    """
    KROME - a package to embed chemistry in astrophysical simulations

    .. [#] Grassi, T.; Bovino, S.; Schleicher, D. R. G.; Prieto, J.; Seifried, D.; Simoncini, E.; Gianturco, F. A., MNRAS, 439, 3, p.2386-2419 [2014MNRAS.439.2386G]
    """

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'kromesph_worker'

    @legacy_function
    def new_particle():
        """
        Add a new particle to Krome.

        Parameters
        ----------
        rho : float
            Density of the particle in units of g*cm**-3.
        u : float
            Internal energy of the particle in units of cm**2*s**-2.
        gamma : float
            Adiabatic index of the particle. Dimensionless.
        mu : float
            Mean molecular weight of the particle in units of g.
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
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.OUT
        )
        for x in ['rho', 'u', 'gamma', 'mu', 'ionrate']:
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
        rho : float
            Density retrieved from the particle, in units of g*cm**-3.
        u : float
            Internal energy retrieved from the particle in units of cm**2*s**-2.
        gamma : float
            Adiabatic index retrieved from the particle. Dimensionless.
        mu : float
            Mean molecular weight retrieved from the particle, in units of g.
        ionrate : float
            Ionization rate retrieved from the particle, in units of s**-1.
        int :
            0 on success.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.IN
        )
        for x in ['rho', 'u', 'gamma', 'mu', 'ionrate']:
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
        rho : float
            Density to set for the particle, in units of g*cm**-3.
        u : float
            Internal energy to set for the particle, in units of cm**2*s**-2.
        gamma : float
            Adiabatic index to set for the particle. Dimensionless.
        mu : float
            Mean molecular weight to set for the particle, in units of g.
        ionrate : float
            Ionization rate to set for the particle, in units of s**-1.

        Returns
        -------
        int :
            0 on success.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.IN
        )
        for x in ['rho', 'u', 'gamma', 'mu', 'ionrate']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @remote_function(can_handle_array=True)
    def get_gamma(index_of_the_particle='i'):
        """
        Retrieve the adiabatic index of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle as returned by `new_particle`.

        Returns
        -------
        gamma : float
            Adiabatic index retrieved from the particle. Dimensionless.
        int :
            0 on success.
        """
        returns (gamma='d')

    @remote_function(can_handle_array=True)
    def set_gamma(index_of_the_particle='i', gamma='d'):
        """
        Set the adiabatic index of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle as returned by `new_particle`.
        gamma : float
            Adiabatic index to set for the particle. Dimensionless.

        Returns
        -------
        int :
            0 on success.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_mu(index_of_the_particle='i'):
        """
        Retrieve the mean molecular weigth of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle as returned by `new_particle`.

        Returns
        -------
        mu : float
            Mean molecular weight retrieved from the particle, in units of g.
        int :
            0 on success.
        """
        returns (mu='d')

    @remote_function(can_handle_array=True)
    def set_mu(index_of_the_particle='i', mu='d'):
        """
        Set the mean molecular weigth of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle as returned by `new_particle`.
        mu : float
            Mean molecular weight retrieved from the particle, in units of g.

        Returns
        -------
        int :
            0 on success.
        """
        returns ()

    @remote_function
    def set_amu_in_g(amu_in_g='d'):
        """
        Set the value of the atomic mass unit in grams.
        This function is used internally by AMUSE to
        send the value of amu in grams to Krome at runtime.

        Parameters
        ----------
        amu_in_g : float
            Value of atomic mass unit in grams.
        """
        returns ()


class KromeSph(ChemicalEvolution, ChemDensityInternalEnergy):
    """
    Krome is a package to embed chemistry in astrophysical simulations.

    KromeSph is an interface for an alternate compilation of Krome
    more suitable for coupled hydrodynamic simulations.
    """
    def __init__(self, unit_converter=None, **options):

        if unit_converter is not None:
            raise Exception('Krome uses predefined units')

        ChemicalEvolution.__init__(self, KromeSphInterface(**options))

        first, last = self.get_firstlast_species_index()
        self.species = dict()
        for i in range(first, last+1):
          self.species[self.get_species_name(i)] = i

        amu_in_g = (1 | u.amu).value_in(u.g)
        self.set_amu_in_g(amu_in_g)

    def define_properties(self, handler):
        handler.add_property('get_time', public_name='model_time')

    def define_methods(self, handler):
        ChemicalEvolution.define_methods(self, handler)
        ChemDensityInternalEnergy.define_methods(self, handler)
        handler.add_method(
            'evolve_model', (u.s,), (handler.ERROR_CODE,)
        )

        handler.add_method(
            'new_particle',
            (
                u.g * u.cm**-3,
                u.cm**2 * u.s**-2,
                handler.NO_UNIT,
                u.g,
                u.s**-1,
            ),
            (handler.INDEX, handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_state',
            (handler.INDEX,),
            (
                u.g * u.cm**-3,
                u.cm**2 * u.s**-2,
                handler.NO_UNIT,
                u.g,
                u.s**-1,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            'set_state',
            (
                handler.INDEX,
                u.g * u.cm**-3,
                u.cm**2 * u.s**-2,
                handler.NO_UNIT,
                u.g,
                u.s**-1,
            ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_gamma',
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_gamma',
            (handler.INDEX, handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_mu',
            (handler.INDEX,),
            (u.g, handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_mu',
            (handler.INDEX, u.g,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_time', (), (u.s, handler.ERROR_CODE,)
        )

    def define_particle_sets(self, handler):
        ChemicalEvolution.define_particle_sets(self, handler)
        ChemDensityInternalEnergy.define_particle_sets(self, handler)
        handler.add_getter('particles', 'get_gamma')
        handler.add_setter('particles', 'set_gamma')
        handler.add_getter('particles', 'get_mu')
        handler.add_setter('particles', 'set_mu')

    def define_state(self, handler):
        ChemicalEvolution.define_state(self, handler)
        handler.add_method('EDIT', 'set_abundance')
        handler.add_method('EDIT', 'set_abundances')
        handler.add_method('EDIT', 'get_abundance')
