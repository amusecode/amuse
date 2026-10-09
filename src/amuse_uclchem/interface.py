import numpy as np
from numpy.typing import NDArray
import uclchem
from uclchem.model import (
    get_species_names as _get_species_names, AbstractModel
)

from amuse.community.interface.chem import (
    ChemNumberDensityTemperature,
    ChemicalEvolution,
    ChemicalEvolutionInterface,
    ChemNumberDensityTemperatureInterface,
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.datamodel import Particle, Particles
from amuse.rfi.core import (
    LegacyFunctionSpecification,
    PythonCodeInterface,
    legacy_function,
    remote_function
)
from amuse.units import units as u


class UclchemImplementation(object):
    def __init__(self):
        """
        Implementation of Uclchem legacy interface functions.

        Parameters
        ----------
        current_time : float
            Current simulation time.
        chem_model : {'cloud', 'collapse', 'cshock', 'jshock', 'prestellarcore'}
            Uclchem model for the chemistry evolution.
        MODEL_MAP : dict[str, type[AbstractModel]]
            Dictionary to map `chem_model` to a Uclchem AbstractModel class.
            These are the actual models which compute the chemistry.
        model_class : type[AbstractModel]
            Current Uclchem AbstractModel class used when calling `evolve_model`.
        self.param_dict : dict
            Dictionary to hold any additional parameters sent to Uclchem before
            evolving. Parameter getters and setters should modify this dictionary.
        uclchem_particles : amuse.datamodel.Particles
            Particles datamodel for storing Uclchem particles.
        _ids : np.ndarray[int]
            Array of unique ids for each particle in `uclchem_particles`.
        _next_particle_id : int
            Next particle index. New ids are assigned by `_get_new_id`.
        """
        self.current_time: float = 0
        self.chem_model: str = 'cloud'
        self.MODEL_MAP = {
            'cloud': uclchem.model.Cloud,
            'collapse': uclchem.model.Collapse,
            'cshock': uclchem.model.CShock,
            'prestellarcore': uclchem.model.PrestellarCore,
            'jshock': uclchem.model.JShock,
        }
        self.model_class = self._validate_chemical_model(
            self.MODEL_MAP.get(self.chem_model, None)
        )
        self.param_dict: dict = {}
        self.uclchem_particles: Particles = Particles()
        self._ids: NDArray = np.empty(0, dtype=np.int64)
        self._next_particle_id: int = 0

    def initialize_code(self) -> int:
        return 0

    def cleanup_code(self) -> int:
        """Remove all the particles stored in Uclchem."""
        self.uclchem_particles.remove_particles(self.uclchem_particles)
        return 0

    def commit_parameters(self) -> int:
        """
        Convert the chemical model string name to its corresponding AbstractModel class.

        Validates that the chemical model specified by the user maps to a valid
        Uclchem model class.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `chem_model` is not a valid model name.
        """
        model = self.MODEL_MAP.get(self.chem_model, None)
        self.model_class = self._validate_chemical_model(model)
        return 0

    def commit_particles(self) -> int:
        """
        Initialize the abundance vector attribute of the particle set
        with an array of zeros.

        Returns
        -------
        int :
            0 on success.
        """
        species = tuple(_get_species_names())
        self.uclchem_particles.add_vector_attribute('abundances', species)
        self.uclchem_particles.abundances = np.zeros(len(species)) #+ 1-30
        return 0

    def recommit_parameters(self) -> int:
        """
        Convert the chemical model name to its corresponding AbstractModel class.

        Validates that the chemical model specified by the user maps to a valid
        Uclchem model class.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `chem_model` is not a valid model name.
        """
        self.commit_parameters()
        return 0

    def recommit_particles(self) -> int:
        """
        Recommit particles.

        Returns
        -------
        int :
            0 on success.
        """
        return 0

    def evolve_model(self, time) -> int:
        """
        Evolve the current model to a specified end time in years.
        The chemistry physics can be changed by calling
        `set_chemical_model`.

        Parameters
        ----------
        time : float
            Time to evolve the chemistry to in years.

        Returns
        -------
        int :
            Return status:
            - 0: Evolved model successfully.
            - -1 : Requested time is less than or equal to
              the current simulation time.
        """
        dt = float(time - self.current_time)
        if dt <= 0:
            return -1

        for particle in self.uclchem_particles:
            params = self._particle_to_dict(particle)
            params['finalTime'] = dt
            params.update(self.param_dict)

            starting_chem = (
                particle.abundances[np.newaxis, :] if
                np.any(particle.abundances) else None
            )
            new_model = self.model_class(
                param_dict=params,
                starting_chemistry=starting_chem
            )
            physics = new_model.physics_array
            particle.number_density = physics[-1, 0, 1]
            particle.temperature = physics[-1, 0, 2]
            particle.abundances = new_model.next_starting_chemistry_array[0, :]

        self.current_time = float(time)
        return 0

    def new_particle(
        self,
        index_of_the_particle,
        number_density,
        temperature,
        ionrate,
        radfield
    ) -> int:
        """
        Add a new particle to Uclchem.

        Parameters
        ----------
        index_of_the_particle : amuse.rfi.python_code.ValueHolder[int]
            Mutable container used to return the index of the new particle.
        number_density : float
            Number density of the particle in units of cm**-3.
        temperature : float
            Temperature of the particle in units of K.
        ionrate : float
            Ionization rate of the particle in units of s**-1.
        radfield : float
            Radiation field of the particle in units of habing.

        Returns
        -------
        int :
            0 on success.
        """
        p = Particle()
        p.number_density = number_density
        p.temperature = temperature
        p.ionrate = ionrate
        p.radfield = radfield

        id = self._get_new_id()
        index_of_the_particle.value = id
        self.uclchem_particles.add_particle(p)
        self._ids = np.append(self._ids, id)

        return 0

    def delete_particle(self, index_of_the_particle) -> int:
        """
        Delete a particle in Uclchem.

        Parameters
        ----------
        index_of_the_particle : int
            index of the particle to delete,
            as returned by 'new_particle`.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)
        self.uclchem_particles.remove_particles(self.uclchem_particles[i].as_set())
        self._ids = np.delete(self._ids, i)
        return 0

    def get_state(
        self,
        index_of_the_particle,
        number_density,
        temperature,
        ionrate,
        radfield
    ) -> int:
        """
        Retrieve the state of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.
        number_density : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the number density
            of the particle in units of cm**-3.
        temperature : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the temperature
            of the particle in units of K.
        ionrate : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the ionization
            rate of the particle in units of s**-1.
        radfield : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the radiation
            field of the particle in units of habing.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        number_density.value = p.number_density
        temperature.value = p.temperature
        ionrate.value = p.ionrate
        radfield.value = p.radfield
        return 0

    def set_state(
        self,
        index_of_the_particle,
        number_density,
        temperature,
        ionrate,
        radfield
    ) -> int:
        """
        Set the state of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.
        number_density : float
            Number density of the particle in units of cm**-3.
        temperature : float
            Temperature of the particle in units of K.
        ionrate : float
            Ionization rate of the particle in units of s**-1.
        radfield : float
            Radiation field of the particle in units of habing.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        p.number_density = number_density
        p.temperature = temperature
        p.ionrate = ionrate
        p.radfield = radfield
        return 0

    def get_number_density(self, index_of_the_particle, number_density) -> int:
        """
        Retrieve the number density of a particle by index.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.
        number_density : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the number density
            of the particle in units of cm**-3.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        number_density.value = p.number_density
        return 0

    def set_number_density(self, index_of_the_particle, number_density) -> int:
        """
        Set the number density of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        number_density: float
            Number density of the particle in units of cm**-3.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        p.number_density = number_density
        return 0

    def get_temperature(self, index_of_the_particle, temperature) -> int:
        """
        Retrieve the temperature of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        temperature : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the temperature
            of the particle in units of K.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        temperature.value = p.temperature
        return 0

    def set_temperature(self, index_of_the_particle, temperature) -> int:
        """
        Set the temperature of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        temperature : float
            Temperature of the particle in units of K.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        p.temperature = temperature
        return 0

    def get_ionrate(self, index_of_the_particle, ionrate) -> int:
        """
        Retrieve the ionization rate of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        ionrate : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the ionization
            rate of the particle in units of s**-1.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        ionrate.value = p.ionrate
        return 0

    def set_ionrate(self, index_of_the_particle, ionrate) -> int:
        """
        Retrieve the ionization rate of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        ionrate : float
            Ionization rate of the particle in units of s**-1.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        p.ionrate = ionrate
        return 0

    def get_radfield(self, index_of_the_particle, radfield) -> int:
        """
        Retrieve the radiation field of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        radfield : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the radiation
            field of the particle in units of habing.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        radfield.value = p.radfield
        return 0

    def set_radfield(self, index_of_the_particle, radfield) -> int:
        """
        Set the radiation field of a particle by index.

        Parameters
        ----------
        index_of_the_particle: int
            Index of the particle as returned by `new_particle`.
        radfield : float
            Radiation field of the particle in units of habing.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        p = self.uclchem_particles[i]
        p.radfield = radfield
        return 0

    def get_abundance(self, index_of_the_particle, species_index, abundance) -> int:
        """
        Retrieve the chemical abundance of a species for a given particle by index.

        The `species_index` can be queried for using the methods `get_species_index`
        and `get_species_name`.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.
        species_index : int
            Index of the species in the abundance array of the particle.
            The `species_index` can be queried for using `get_species_index`.
        abundance : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the abundance of a particle.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)
        abundance.value = self.uclchem_particles[i].abundances[species_index]
        return 0

    def set_abundance(self, index_of_the_particle, species_index, abundance) -> int:
        """
        Set the chemical abundance of a species for a given particle by index.

        The `species_index` can be queried for using the methods `get_species_index`
        and `get_species_name`.

        Parameters
        ----------
        index_of_the_particle : int
            Index of the particle as returned by `new_particle`.
        species_index : int
            Index of the species in the abundance array of the particle.
            The `species_index` can be queried for using `get_species_index`.
        abundance : np.ndarray[float]
            Abundance of the chemical species of the particle.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        """
        i = self._id_to_storage_index(index_of_the_particle)

        abundances = self.uclchem_particles[i].abundances
        abundances[species_index] = abundance
        self.uclchem_particles[i].abundances = abundances
        return 0

    def set_abundances(self, index_of_the_particle, abundances, N) -> int:
        """
        Set the full chemical abundance array for a given particle.

        Parameters
        ----------
        index_of_the_particle : np.ndarray[int]
            Index of the particle as returned by `new_particle`. If an array,
            only the first element is used.
        abundances : np.ndarray[float]
            Abundance values to assign to the particle, of length `N`.
        N : int
            Number of abundance values in `abundances`. Must match the
            particle's existing abundance array length.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not a valid particle id.
        ValueError :
            If `N` does not match the particle's existing number of abundances.
        """
        if not isinstance(index_of_the_particle, int):
            index_of_the_particle = index_of_the_particle[0]
        i = self._id_to_storage_index(index_of_the_particle)
        N_abundances = self.uclchem_particles[i].abundances.shape[0]
        if N != N_abundances:
            raise ValueError(
                f'abundances must have shape {N_abundances}, got {N}!'
            )
        self.uclchem_particles[i].abundances = abundances
        return 0

    def get_firstlast_species_index(self, first, last) -> int:
        """
        Get the index of the first and last species inside Uclchem.

        This is a helper method for accessing the abundance array as
        `instance.particles.abundances`.

        Parameters
        ----------
        first : amuse.rfi.python_code.ValueHolder[int]
            Index of the first species in the abundance array.
        last : amuse.rfi.python_code.ValueHolder[int]
            Index of the last species in the abundance array.

        Returns
        -------
        int :
            0 on success.
        """
        first.value = 0
        last.value = len(_get_species_names()) - 1
        return 0

    def get_species_index(self, name, species_index) -> int:
        """
        Given the name of a chemical species in the
        chemical abundance array, retrieve its index.

        Chemical abundances for each particle are stored
        as a 1D array, where each element corresponds to
        the abundance of a particular species.

        Parameters
        ----------
        name : str
            Name of chemical species. Must be one of the
            species tracked by Uclchem.
        species_index : amuse.rfi.python_code.ValueHolder[int]
            Mutable container used to return the index
            of the species.

        Returns
        -------
        int :
            Return status:
            - 0: Species index retrieved successfully.
            - -1: Species does not exist.

        Examples
        --------
        >>> chem = Uclchem()
        >>> chem.get_species_index('H2O')
        31
        """
        species_names = _get_species_names()

        try:
            idx = species_names.index(name)
        except ValueError:
            species_index.value = -1
            return -1

        species_index.value = idx
        return 0

    def get_species_name(self, species_index, name) -> int:
        """
        Given the index of a chemical species in the
        chemical abundance array, retrieve its name.

        Chemical abundances for each particle are stored
        as a 1D array, where each element corresponds to
        the abundance of a particular species.

        Parameters
        ----------
        species_index : int
            Index of the chemical species in the abundance array.
        name : amuse.rfi.python_code.ValueHolder[str]
            Mutable container used to return the name of
            the chemical species.

        Returns
        -------
        int :
            0 on success.

        Examples
        --------
        >>> chem = Uclchem()
        >>> chem.get_species_name(31)
        'H2O'
        """
        species_names = _get_species_names()
        if not 0 <= species_index < len(species_names):
            name.value = -1
            return -1

        name.value = species_names[species_index]
        return 0

    def get_time(self, time) -> int:
        """
        Retrieve current model time in years.

        Parameters
        ----------
        time : amuse.rfi.python_code.ValueHolder[float]
            Mutable container used to return the
            current time in units of years.

        Returns
        -------
        int :
            0 on success.
        """
        time.value = self.current_time
        return 0

    def get_number_of_particles(self, number_of_particles) -> int:
        """
        Retrieve the current number of particles in the code.

        Parameters
        ----------
        number_of_particles : amuse.rfi.python_code.ValueHolder[int]
            Mutable container used to return the
            current number of particles.

        Returns
        -------
        int :
            0 on success.
        """
        number_of_particles.value = len(self.uclchem_particles)
        return 0

    def get_chemical_model(self, chem_model) -> int:
        """
        Retrieve the chemical model type used by Uclchem.
        This is the physics model used internally by Uclchem
        to evolve the chemistry.

        Possible values are `'cloud'`, `'collapse'`, `'cshock'`,
        `'jshock'`, and `'prestellarcore'`.

        Parameters
        ----------
        chem_model : amuse.rfi.python_code.ValueHolder[str]
            Mutable container used to return the current
            model type.

        Returns
        -------
        int :
            0 on success.
        """
        chem_model.value = self.chem_model
        return 0

    def set_chemical_model(self, chem_model) -> int:
        """
        Set the chemical model used by Uclchem. This is
        the physics model used internally by Uclchem to
        evolve the chemistry.

        Possible values are `'cloud'`, `'collapse'`, `'cshock'`,
        `'jshock'`, and `'prestellarcore'`.

        Parameters
        ----------
        chem_model : {'cloud', 'collapse', 'cshock', 'jshock', 'prestellarcore'}
            Chemical model specifying the chemistry physics when evolving the particles.

        Returns
        -------
        int :
            0 on success.

        Raises
        ------
        ValueError :
            If `chem_model` is not one of the allowed models.
        """
        self.chem_model = chem_model
        self.model_class = self._validate_chemical_model(
            self.MODEL_MAP.get(self.chem_model, None)
        )
        return 0

    def _validate_chemical_model(self, model: type[AbstractModel] | None) -> type[AbstractModel]:
        """
        Validate chemical model class.

        Parameters
        ----------
        model : type[AbstractModel] | None
            Model class to validate.

        Returns
        -------
        type[AbstractModel] :
            Validated model class.

        Raises
        ------
        ValueError :
            If `model` is `None`.
        """
        if model is None:
            options = ", ".join(f"'{k}'" for k in self.MODEL_MAP.keys())
            raise ValueError(
                f"\nchem_model must be one of the following options: {options}! "
                f"Got {self.chem_model!r}."
            )
        return model

    def _id_to_storage_index(self, index_of_the_particle: int) -> int:
        """
        Resolve a particle id to its current index in internal storage.

        Ids are assigned once per particle by `_get_new_id` and never reused.
        `self._ids[k]` holds the id of the particle stored at position `k` in
        `self.uclchem_particles`; this mapping is kept in sync by `new_particle`
        and `delete_particle`. Because storage positions shift when a particle
        is deleted while ids remain fixed, this lookup must be used instead of
        assuming id and storage index are equal.

        Parameters
        ----------
        index_of_the_particle : int
            Id of the particle, as assigned in `new_particle` and returned via
            its `index_of_the_particle` output parameter.

        Returns
        -------
        int :
            Current index of the particle in internal storage
            (`self.uclchem_particles`).

        Raises
        ------
        ValueError :
            If no particle with the given id exists in `self._ids`.
        """
        idx = np.where(self._ids == index_of_the_particle)[0]
        if idx.size == 0:
            raise ValueError(f'Particle id: {index_of_the_particle} not found!')

        return int(idx[0])

    def _get_new_id(self) -> int:
        """Create a unique, monotonically increasing id for `new_particle`."""
        new_id = self._next_particle_id
        self._next_particle_id += 1

        return int(new_id)

    def _particle_to_dict(self, particle: Particle) -> dict:
        """
        Format an AMUSE Particle into a parameter dictionary
        readable by Uclchem.

        Unlike AMUSE, Uclchem does not work with particles but
        rather a dictionary of parameters. This is a helper function
        to be called before evolving a particle with Uclchem.

        Parameters
        ----------
        particle : amuse.datamodel.Particle
            Particle to evolve by Uclchem.

        Returns
        -------
        dict :
            Particle returned as a dictionary with key, value
            pairs taken from the particle attributes.
        """
        param_dict = {
            'initialDens': float(particle.number_density),
            'initialTemp': float(particle.temperature),
            'zeta': float(particle.ionrate) / 1.3e-17,
            'radfield': float(particle.radfield),
        }
        return param_dict


class UclchemInterface(
    ChemicalEvolutionInterface,
    ChemNumberDensityTemperatureInterface,
    PythonCodeInterface,
    LiteratureReferencesMixIn,
):
    """
    UCLCHEM: a gas-grain chemical code for astrochemical modelling.

    .. [#] ADS:2017AJ....154...38H (Holdship, J. ; Viti, S, ; Jiménez-Serra, I.; Makrymallis, A. ; Priestley, F. , 2017, AJ)
    """
    def __init__(self, **options):
        PythonCodeInterface.__init__(
            self,
            UclchemImplementation,
            'uclchem_worker',
            **options
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.OUT)
        for x in ['number_density', 'temperature', 'ionrate', 'radfield']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for x in ['number_density', 'temperature', 'ionrate', 'radfield']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='i', direction=function.IN)
        for x in ['number_density', 'temperature', 'ionrate', 'radfield']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @remote_function(can_handle_array=True)
    def get_radfield(index_of_the_particle='i'):
        returns (radfield='d')

    @remote_function(can_handle_array=True)
    def set_radfield(index_of_the_particle='i', radfield='d'):
        returns ()

    @remote_function
    def get_chemical_model():
        returns (chem_model='s')

    @remote_function
    def set_chemical_model(chem_model='s'):
        returns ()


class Uclchem(ChemicalEvolution, ChemNumberDensityTemperature):
    """Uclchem is a gas-grain chemical code for astrochemical modelling."""
    def __init__(self, unit_converter=None, **options):

        if unit_converter is not None:
            raise ValueError('Uclchem uses predefined units')

        chem_interface = UclchemInterface(**options)

        ChemicalEvolution.__init__(
            self,
            chem_interface
        )

        first, last = self.get_firstlast_species_index()
        self.species = dict()
        for i in range(first, last+1):
          self.species[self.get_species_name(i)] = i

    def define_methods(self, handler):
        ChemicalEvolution.define_methods(self, handler)
        ChemNumberDensityTemperature.define_methods(self, handler)
        handler.add_method(
            'new_particle',
            (u.cm**-3, u.K, u.s**-1, u.habing),
            (handler.INDEX, handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_state',
            (handler.INDEX,),
            (u.cm**-3, u.K, u.s**-1, u.habing, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_state',
            (handler.INDEX, u.cm**-3, u.K, u.s**-1, u.habing,),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            'get_radfield',
            (handler.INDEX,),
            (u.habing, handler.ERROR_CODE,),
        )

        handler.add_method(
            'set_radfield',
            (handler.INDEX, u.habing,),
            (handler.ERROR_CODE,),
        )

    def define_parameters(self, handler):
        handler.add_method_parameter(
            'get_chemical_model',
            'set_chemical_model',
            'chem_model',
            "'cloud', 'collapse', 'cshock', 'jshock', 'prestellarcore'",
            default_value='cloud'
        )
        handler.add_interface_parameter(
            'out_species', 'Array of molecules to use', default_value=['H', 'H2']
        )

    def define_particle_sets(self, handler):
        ChemicalEvolution.define_particle_sets(self, handler)
        ChemNumberDensityTemperature.define_particle_sets(self, handler)

    def define_state(self, handler):
        ChemicalEvolution.define_state(self, handler)
        handler.add_method('RUN', 'set_abundance')
        handler.add_method('RUN', 'set_abundances')
        handler.add_method('RUN', 'get_abundance')
