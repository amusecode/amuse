"""
Chemical Evolution Interface Definition
"""

from collections.abc import Sequence

from amuse.community.interface import common
from amuse.rfi.core import (
    LegacyFunctionSpecification,
    legacy_function,
    remote_function
)
from amuse.units import units as u
import numpy as np
from numpy.typing import NDArray


class ChemicalEvolutionInterface(common.CommonCodeInterface):

    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.

        Returns
        -------
        int :
            Return status:
            - 0: Model is initialized and evolution can start.
            - -1 : An error occurred during initialization. The error
              needs to be further specified by each code implementation.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def recommit_particles():
        """
        Let the code perform initialization actions
        after the number of particles have been updated
        or particle attributes have been updated from
        the script.

        Returns
        -------
        int :
            Return status:
            - 0: Model is initialized and evolution can start.
            - -1: An error occurred during initialization. The error
              needs to be further specified by each code implementation.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def evolve_model():
        """
        Evolve the model until the given time, or until a stopping
        condition is set. The model will be evolved until this time
        is reached exactly or just after.

        Parameters
        ----------
        time : float
            Time to evolve the system to.

        Returns
        -------
        int :
            Return status:
            - 0: Evolved model successfully.
            - -1 : Requested time could not be evolved to.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def delete_particle():
        """
        Remove the definition of particle from the code. After calling this
        function the particle is no longer part of the model evolution. It is
        up to the code if the index will be reused.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle to delete.

        Returns
        -------
        int :
            Return status:
            - 0 : Particle was removed from the model.
            - -1 : Particle could not be removed.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.IN
        )
        function.result_type = 'i'
        return function

    @remote_function(can_handle_array=True)
    def get_ionrate(index_of_the_particle='i'):
        """
        Retrieve the ionization rate of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.

        Returns
        -------
        ionrate : float
            Ionrate of the particle in units of s**-1.
        """
        returns (ionrate='d')

    @remote_function(can_handle_array=True)
    def set_ionrate(index_of_the_particle='i', ionrate='d'):
        """
        Set the ionization rate of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        ionrate : float
            Ionrate of the particle in units of s**-1.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_abundance(index_of_the_particle='i', species_index='i'):
        """
        Retrieve the chemical abundance of a species by index for a given particle.

        The `species_index` can be queried for using the methods `get_species_index`
        and `get_species_name`.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        species_index : int
            Index of the species in the abundance array.

        Returns
        -------
        abundance : float
            Abundance of the requested species.

        Examples
        --------
        # get H2O abundance of particle 0
        >>> get_species_index('H2O')
        2
        >>> get_abundance(0, 2)
        0.005
        """
        returns (abundance='d')

    @remote_function(can_handle_array=True)
    def set_abundance(
        index_of_the_particle='i', species_index='i', abundance='d'
    ):
        """
        Set the chemical abundance of a species by index for a given particle.

        The `species_index` can be queried for using the methods `get_species_index`
        and `get_species_name`.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        species_index : int
            Index of the species in the abundance array.
        abundance : float
            Abundance of the species to set.

        Examples
        --------
        # set H2O abundance of particle 0
        >>> get_species_index('H2O')
        2
        >>> set_abundance(0, 2, 0.005)
        """
        returns ()

    @remote_function(must_handle_array=True)
    def set_abundances(index_of_the_particle='i', abundances='d'):
        """
        Set all the abundances of a particle from an array.
        The abundance array must match the shape of the abundance
        array in the chemical evolution code.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        abundances : np.ndarray
            Array of abundances. Must match the shape of
            the abundance array inside the chemistry code.

        Examples
        --------
        >>> abundances = np.random.rand(100)
        >>> chem.set_abundances(0, abundances)

        Notes
        -----
        # to get the abundance array of a particle, use:
        >>> chem.particles[0].abundances
        [0.0005, 0.0004, ...]
        """
        returns ()

    @remote_function()
    def get_firstlast_species_index():
        """
        Retrieve the index bounds of the chemical abundance array.

        Used internally by the ChemicalEvolutionInterface.

        Returns
        -------
        first : int
            Index of the first species in the abundance array.
        last : int
            Index of the last species in the abundance array.
        """
        returns (first='i', last='i')

    @remote_function(can_handle_array=True)
    def get_species_index(name='s'):
        """
        Given the name of a chemical species in the
        chemical abundance array, retrieve its index.

        Chemical abundances for each particle are stored
        as a 1D array, where each element corresponds to
        the abundance of a particular species.

        Parameters
        ----------
        name : str
            Name of species.

        Returns
        -------
        species_index : int
            Index of the species in the abundance array.

        Examples
        --------
        >>> chem.get_species_index('H')
        0
        """
        returns (species_index='i')

    @remote_function(can_handle_array=True)
    def get_species_name(species_index='i'):
        """
        Given the index of a chemical species in the
        chemical abundance array, retrieve its name.

        Chemical abundances for each particle are stored
        as a 1D array, where each element corresponds to
        the abundance of a particular species.

        Parameters
        ----------
        species_index : int
            Index of the species.

        Returns
        -------
        name : str
            Species name corresponding to the input index.

        Examples
        --------
        >>> chem.get_species_name(0)
        'H'
        """
        returns (name='s')

    @remote_function
    def get_time():
        """
        Retrieve the model time. This time should be close to the end time
        specified in the evolve code.

        Returns
        -------
        time : float
            Current model time.
        """
        returns (time='d')

    @remote_function
    def get_number_of_particles():
        """
        Retrieve the total number of particles defined in the code.

        Returns
        -------
        number_of_particles : int
            Current number of particles in the model.

        Notes
        -----
        `number_of_particles` may not be updated if particles are
        deleted. To ensure the correct number is returned, call
        `commit_particles` or `recommit_particles` before calling
        `get_number_of_particles`.
        """
        returns (number_of_particles='i')


class ChemNumberDensityTemperatureInterface:
    """
    Mixin for number density and temperature setters
    and getters in ChemicalEvolutionInterface.
    """
    @remote_function(can_handle_array=True)
    def get_number_density(index_of_the_particle='i'):
        """
        Retrieve the number density of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.

        Returns
        -------
        number_density : float
            Number density of the particle in units of cm**-3.
        """
        returns (number_density='d')

    @remote_function(can_handle_array=True)
    def set_number_density(index_of_the_particle='i', number_density='d'):
        """
        Set the number density of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        number_density : float
            Number density of the particle in units of cm**-3.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_temperature(index_of_the_particle='i'):
        """
        Retrieve the temperature of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.

        Returns
        -------
        temperature : float
            Temperature of the particle in units of Kelvin.
        """
        returns (temperature='d')

    @remote_function(can_handle_array=True)
    def set_temperature(index_of_the_particle='i', temperature='d'):
        """
        Set the temperature of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle.
        temperature : float
            Temperature of the particle in units of Kelvin.
        """
        returns ()


class ChemDensityInternalEnergyInterface:
    """
    Mixin for density and internal energy setters
    and getters in ChemicalEvolutionInterface.
    """
    @remote_function(can_handle_array=True)
    def get_density(index_of_the_particle='i'):
        """
        Get the density of a particle.

        Parameters
        ----------
        index_of_the_particle : int
             Index of the particle to retrieve the density of.

        Returns
        -------
        rho : float
            Density retrieved from the particle, in units of g*cm**-3.
        """
        returns (rho='d')

    @remote_function(can_handle_array=True)
    def set_density(index_of_the_particle='i', rho='d'):
        """
        Set the density of a particle.

        Parameters
        ----------
        index_of_the_particle : int
             Index of the particle to set the density of.
        rho : float
            Density to set for the particle, in units of g*cm**-3.
        """
        returns ()

    @remote_function(can_handle_array=True)
    def get_internal_energy(index_of_the_particle='i'):
        """
        Get the internal energy of a particle.

        Parameters
        ----------
        index_of_the_particle : int
             Index of the particle to retrieve the internal energy of.

        Returns
        -------
        u : float
            Internal energy retrieved from the particle, in units of cm**2*s**-2.
        """
        returns (u='d')

    @remote_function(can_handle_array=True)
    def set_internal_energy(index_of_the_particle='i', u='d'):
        """
        Set the internal energy of a particle.

        Parameters
        ----------
        index_of_the_particle : int
             Index of the particle to set the internal energy of.
        u : float
            Internal energy to set for the particle, in units of cm**2*s**-2.
        """
        returns ()


class ChemicalEvolution(common.CommonCode):

    def __init__(self, legacy_interface, unit_converter=None, **options):
        self.unit_converter = unit_converter

        common.CommonCode.__init__(self, legacy_interface, **options)

    def get_abundances_by_name(
        self,
        index_of_the_particle: int | Sequence[int],
        species_names: str | Sequence[str]
    ) -> NDArray[np.float64]:
        """
        Get the abundances of one or more particles at the current
        simulation time, by species name.

        Both singular and array inputs are permitted for both
        `index_of_the_particle` and `species_names`.

        Parameters
        ----------
        index_of_the_particle : int | Sequence[int]
            Indices of the particle(s) as returned by `new_particle`.
        species_names : str | Sequence[str]
            Species name(s) to query. Each name must be a species
            tracked by the chemistry code network.

        Returns
        -------
        abundances : NDArray[np.float64]
            Array of shape (N_particle, N_species) containing the
            abundance of each queried species for each queried
            particle. Always 2D, regardless of whether scalar or
            sequence inputs were passed.

        Raises
        ------
        ValueError
            If `index_of_the_particle` or `species_names` is empty.

        Examples
        --------
        >>> chem = ChemicalEvolution()
        >>> chem.particles.add_particles(particles)
        >>> chem.get_abundances_by_names([0, 1], ['H', 'H2'])
        array([[1.00000000e-40, 1.00000000e-40],
              [1.00000000e-40, 1.00000000e-40]])

        >>> chem.get_abundances_by_names(0, ['H', 'H2'])
        array([[1.00000000e-40, 1.00000000e-40]])

        Notes
        -----
        To obtain a dictionary of each (species: index) in a chemistry
        code:
        >>> chem = ChemicalEvolution()
        >>> chem.species
        {'E': 0, 'H-': 1, 'H': 2, 'HE': 3, 'H2': 4, ...}
        """
        particle_ids = np.asarray(index_of_the_particle).ravel()
        species_names_arr = np.asarray(species_names).ravel()

        if particle_ids.size == 0:
            raise ValueError('index_of_the_particle must not be empty.')
        if species_names_arr.size == 0:
            raise ValueError('species_names must not be empty.')

        species_indeces = self.get_species_index(species_names_arr)

        particle_idx_grid, species_idx_grid = np.meshgrid(
            particle_ids, species_indeces, indexing='ij'
        )

        return self.get_abundance(
            particle_idx_grid.ravel(), species_idx_grid.ravel()
        ).reshape(particle_idx_grid.shape)

    def define_properties(self, handler):
        handler.add_property('get_time', public_name='model_time')

    def define_methods(self, handler):
        common.CommonCode.define_methods(self, handler)

        handler.add_method('evolve_model', (u.yr,), (handler.ERROR_CODE,))

        handler.add_method(
            'delete_particle', (handler.INDEX,), (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_ionrate',
            (handler.INDEX,),
            (u.s**-1, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_ionrate',
            (handler.INDEX, u.s**-1,),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_abundance',
            (handler.INDEX, handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_abundance',
            (handler.INDEX, handler.INDEX, handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_abundances',
            (handler.INDEX, handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_firstlast_species_index',
            (),
            (handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_species_index',
            (handler.NO_UNIT,),
            (handler.INDEX, handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_species_name',
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_time', (), (u.yr, handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_number_of_particles',
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,),
        )

    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_getter('particles', 'get_state')
        handler.add_setter('particles', 'set_ionrate')
        handler.add_getter('particles', 'get_ionrate')
        handler.add_gridded_getter(
            'particles',
            'get_abundance',
            'get_firstlast_species_index',
            names=('abundances',),
        )
        handler.add_gridded_setter(
            'particles',
            'set_abundance',
            'get_firstlast_species_index',
            names=('abundances',),
        )

    def define_state(self, handler):
        common.CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED', 'EDIT', 'commit_parameters')
        handler.add_transition('RUN', 'PARAMETER_CHANGE_A', 'invoke_state_change2')
        handler.add_transition('EDIT', 'PARAMETER_CHANGE_B', 'invoke_state_change2')
        handler.add_transition('PARAMETER_CHANGE_A', 'RUN', 'recommit_parameters')
        handler.add_transition('PARAMETER_CHANGE_B', 'EDIT', 'recommit_parameters')
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_method('UPDATE', 'set_abundance')
        handler.add_method('UPDATE', 'set_abundances')
        handler.add_method('UPDATE', 'get_abundance')
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_state')
        handler.add_method('RUN', 'get_abundance')


class ChemNumberDensityTemperature:
    """
    Mixin for number density and temperature setters
    and getters in ChemicalEvolution.
    """
    def define_methods(self, handler):
        handler.add_method(
            'get_number_density',
            (handler.INDEX,),
            (u.cm**-3, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_number_density',
            (handler.INDEX, u.cm**-3,),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_temperature',
            (handler.INDEX,),
            (u.K, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_temperature',
            (handler.INDEX, u.K,),
            (handler.ERROR_CODE,),
        )

    def define_particle_sets(self, handler):
        handler.add_setter('particles', 'set_number_density')
        handler.add_getter('particles', 'get_number_density')
        handler.add_setter('particles', 'set_temperature')
        handler.add_getter('particles', 'get_temperature')


class ChemDensityInternalEnergy:
    """
    Mixin for density and internal energy setters
    and getters in ChemicalEvolution.
    """

    def define_methods(self, handler):
        handler.add_method(
            'get_density',
            (handler.INDEX,),
            (u.g * u.cm**-3, handler.ERROR_CODE),
        )

        handler.add_method(
            'set_density',
            (handler.INDEX, u.g * u.cm**-3),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_internal_energy',
            (handler.INDEX,),
            (u.cm**2 * u.s**-2, handler.ERROR_CODE),
        )

        handler.add_method(
            'set_internal_energy',
            (handler.INDEX, u.cm**2 * u.s**-2),
            (handler.ERROR_CODE,),
        )

    def define_particle_sets(self, handler):
        handler.add_setter('particles', 'set_density')
        handler.add_getter('particles', 'get_density')
        handler.add_setter('particles', 'set_internal_energy')
        handler.add_getter('particles', 'get_internal_energy')
