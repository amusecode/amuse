import numpy as np
from numpy.typing import NDArray
import tsunami

from amuse.community.interface.common import CommonCode
from amuse.community.interface.gd import (
    GravitationalDynamics,
    GravitationalDynamicsInterface,
    GravityFieldCode,
    GravityFieldInterface
)
from amuse.datamodel.particles import Particles
from amuse.rfi.core import (
    LegacyFunctionSpecification,
    PythonCodeInterface,
    legacy_function,
    remote_function
)
from amuse.rfi.python_code import ValueHolder
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.units import (
    generic_unit_system, nbody_system as ns, units as u
)


class TsunamiImplementation(object):
    """Implementation class for Tsunami Interface."""
    def __init__(self):
        # hard code N-Body units to: Lscale = 1 [AU], Mscale = 1 [MSun]
        self.tsunami = tsunami.Tsunami(1, 1)

        # temporary buffers for staging particles
        self._id_list: list[int] = []
        self._mass_list: list[float] = []
        self._radius_list: list[float] = []
        self._pos_list: list[list[float]] = []
        self._vel_list: list[list[float]] = []
        self._spin_list: list[list[float]] = []

        # commited particle arrays
        self._ids = np.empty(0, dtype=np.int32)
        self._mass = np.empty(0, dtype=np.float64)
        self._radius = np.empty(0, dtype=np.float64)
        self._pos = np.empty((0, 3), dtype=np.float64)
        self._vel = np.empty((0, 3), dtype=np.float64)
        self._spin = np.empty((0, 3), dtype=np.float64)
        self._stype = np.empty(0, dtype=np.int64)

        self._next_particle_id: int = 0

    def initialize_code(self) -> int:
        """
        Disable any TSUNAMI features that are not needed in AMUSE.

        Tsunami by default will recenter the system to its center of mass
        frame by subtracting it. External forces must be set to `True` so
        that Tsunami adds the center of mass back into the system.
        """
        self.tsunami.Conf.wExt = True
        return 0

    def cleanup_code(self) -> int:
        """Deallocate Tsunami particles."""
        self._clear_temporary_particle_buffers()
        return 0

    def commit_parameters(self) -> int:
        """
        Commit Tsunami parameters.

        Notes
        -----
        * `tsunami.commit_particles` should not be called
        here as it introduced subtle bugs in setting
        parameters. That method should be reserved for
        `recommit_parameters`.
        """
        return 0

    def commit_particles(self) -> int:
        """
        Add all particles stored in the particle buffers
        into the TSUNAMI code. TSUNAMI does not support
        adding individual particles and instead requires
        preallocated np.ndarrays of particles.

        This method grabs all pre-existing particles in TSUNAMI
        as well as any new particles inside the temporary buffers
        and formats them as numpy arrays to be read by TSUNAMI.

        This method should only be called if at least 2 particles exist
        between committed and buffered particles. Calling this method
        with no new particles in the buffers does nothing.

        Returns
        -------
         0 :
            If particles were added to TSUNAMI successfully or if no new
            particles are available to be added.
        -1 :
            If the temporary particle buffers have inconsistent lengths.

        Raises
        ------
        ValueError :
            If the total particle count is less than 2 particles between
            the pre-existing particles and the particles in the buffer.
        """
        buffer_lengths = {
            len(self._id_list),
            len(self._mass_list),
            len(self._radius_list),
            len(self._pos_list),
            len(self._vel_list),
            len(self._spin_list),
        }
        # if no new particles to add
        if len(buffer_lengths) == 0:
            return 0
        # if shape mismatch between particle buffers
        if len(buffer_lengths) != 1:
            return -1

        N_existing: int = self._pos.shape[0]
        N_new: int = buffer_lengths.pop()
        N_total: int = N_existing + N_new

        if N_new == 0:
            return 0

        if N_total < 2:
            raise ValueError(
                'Tsunami needs at least 2 particles when committing!'
            )

        if not (
            len(self._id_list) == N_new and
            len(self._mass_list) == N_new and
            len(self._radius_list) == N_new and
            len(self._vel_list) == N_new and
            len(self._spin_list) == N_new
        ):
            return -1

        # preallocate particleset for Tsunami
        ids = np.empty(N_total, dtype=np.int32)
        mass = np.empty(N_total, dtype=np.float64)
        radius = np.empty(N_total, dtype=np.float64)
        pos = np.empty((N_total, 3), dtype=np.float64)
        vel = np.empty((N_total, 3), dtype=np.float64)
        spin = np.empty((N_total, 3), dtype=np.float64)
        stype = np.ones(N_total, dtype=np.int64) * -1

        # add any existing Tsunami particles
        if N_existing != 0:
            ids[:N_existing] = self._ids
            mass[:N_existing] = self._mass
            radius[:N_existing] = self._radius
            pos[:N_existing, :] = self._pos
            vel[:N_existing, :] = self._vel
            spin[:N_existing, :] = self._spin

        # add new particles
        ids[N_existing:] = np.asarray(self._id_list, dtype=np.int32)
        mass[N_existing:] = np.asarray(self._mass_list, dtype=np.float64)
        radius[N_existing:] = np.asarray(self._radius_list, dtype=np.float64)
        pos[N_existing:, :] = np.asarray(self._pos_list, dtype=np.float64).reshape(-1, 3)
        vel[N_existing:, :] = np.asarray(self._vel_list, dtype=np.float64).reshape(-1, 3)
        spin[N_existing:, :] = np.asarray(self._spin_list, dtype=np.float64).reshape(-1, 3)

        self.tsunami.add_particle_set(pos, vel, mass, radius, stype, spin)

        self._ids = ids
        self._mass = mass
        self._radius = radius
        self._pos = pos
        self._vel = vel
        self._spin = spin
        self._stype = stype
        self._clear_temporary_particle_buffers()
        self.synchronize_model()

        return 0

    def recommit_parameters(self) -> int:
        """Recommit parameters after committing them."""
        self.tsunami.commit_parameters()
        return 0

    def recommit_particles(self) -> int:
        """Recommit particles after committing them."""
        self.commit_particles()
        return 0

    def evolve_model(self, time: float) -> int:
        """
        Evolve the system to a specified final time.

        Parameters
        ----------
        time : float
            Time to evolve to.

        Returns
        -------
         0 :
            If the evolution is successful.
        -1 :
            If `time` <= current model time.
        """
        if time <= self.tsunami.time:
            return -1

        self.tsunami.evolve_system(time)
        self.synchronize_model()

        return 0

    def evolve_model_dtmax(self, time: float) -> int:
        """
        Evolves the system for a single integration step,
        ensuring that the system time does not exceed the
        specified target `time`. This method is especially
        useful when precise control over the maximum timestep
        is needed, such as in visualizations or applications
        requiring fixed time intervals.

        Unlike `tsunami.evolve_system`, this method performs only
        one integration step and adapts the timestep to ensure it
        does not exceed `time`. If the integration requires a
        shorter timestep, it may stop well before `time`.

        Parameters
        ----------
        time : float
            Time to evolve to.

        Returns
        -------
         0 :
            If the evolution is successful.
        -1 :
            If `time` <= current model time.
        """
        if time <= self.tsunami.time:
            return -1

        self.tsunami.evolve_system_dtmax(time)
        self.synchronize_model()

        return 0

    def synchronize_model(self) -> int:
        """
        Synchronize AMUSE with TSUNAMI. This ensures that
        the current particle state in AMUSE matches the
        current state of TSUNAMI.

        Synchronizes the masses, radii, positions, velocities,
        and spins of the system.

        Call this function before accessing values to ensure
        that the returned values are up to date.
        """
        self.tsunami.sync_masses(self._mass)
        self.tsunami.sync_radii(self._radius)
        self.tsunami.sync_internal_state(self._pos, self._vel, self._spin)
        return 0

    def new_particle(
        self,
        index_of_the_particle: ValueHolder,
        mass: float,
        radius: float,
        x: float, y: float, z: float,
        vx: float, vy: float, vz: float,
        wx: float, wy: float, wz: float
    ) -> int:
        """
        Add a new particle to TSUNAMI.

        The particles are added to a temporary buffer,
        and are only added to TSUNAMI upon calling
        `commit_particles` or `recommit_particles`.

        Parameters
        ----------
        index_of_the_particle : ValueHolder[int]
             ValueHolder instance to return the index
             of the new particle.
        mass : float
            Particle mass.
        radius : float
            Particle radius.
        x, y, z : float
            Particle position.
        vx, vy, vz : float
            Particle velocity.
        wx, wy, wz : float
            Particle spin.

        Returns
        -------
        0 :
            Particle was created successfully.
        """
        self._mass_list.append(mass)
        self._radius_list.append(radius)
        self._pos_list.append([x, y, z])
        self._vel_list.append([vx, vy, vz])
        self._spin_list.append([wx, wy, wz])

        id = self._get_new_id()
        index_of_the_particle.value = id
        self._id_list.append(id)

        return 0

    def delete_particle(self, index_of_the_particle: int) -> int:
        """
        Delete a particle in TSUNAMI.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.

        Returns
        -------
        0 :
            Particle was deleted successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        ValueError :
            If deleting a particle would cause TSUNAMI
            to have less than 2 particles.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)

        if len(self._pos)-1 < 2:
            raise ValueError(
                'TSUNAMI does not support less than 2 particles!'
            )

        self._ids = np.delete(self._ids, i)
        self._mass = np.delete(self._mass, i)
        self._radius = np.delete(self._radius, i)
        self._pos = np.delete(self._pos, i, axis=0)
        self._vel = np.delete(self._vel, i, axis=0)
        self._spin = np.delete(self._spin, i, axis=0)
        self._stype = np.delete(self._stype, i)

        self.tsunami.remove_particle(i)
        self.synchronize_model()

        return 0

    def get_state(
        self,
        index_of_the_particle: int,
        mass: ValueHolder,
        radius: ValueHolder,
        x: ValueHolder, y: ValueHolder, z: ValueHolder,
        vx: ValueHolder, vy: ValueHolder, vz: ValueHolder,
        wx: ValueHolder, wy: ValueHolder, wz: ValueHolder
    ) -> int:
        """
        Retrieve the state of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        mass : ValueHolder[float]
            ValueHolder instance to return the mass value.
        radius : ValueHolder[float]
            ValueHolder instance to return the radius value.
        x, y, z : ValueHolder[float]
            ValueHolder instance to return the particle position.
        vx, vy, vz : ValueHolder[float]
            ValueHolder instance to return the particle velocity.
        wx, wy, wz : ValueHolder[float]
            ValueHolder instance to return the particle spin.

        Returns
        -------
        0 :
            State was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        mass.value = self._mass[i]
        radius.value = self._radius[i]
        x.value = self._pos[i, 0]
        y.value = self._pos[i, 1]
        z.value = self._pos[i, 2]
        vx.value = self._vel[i, 0]
        vy.value = self._vel[i, 1]
        vz.value = self._vel[i, 2]
        wx.value = self._spin[i, 0]
        wy.value = self._spin[i, 1]
        wz.value = self._spin[i, 2]

        return 0

    def set_state(
        self,
        index_of_the_particle: int,
        mass: float,
        radius: float,
        x: float, y: float, z: float,
        vx: float, vy: float, vz: float,
        wx: float, wy: float, wz: float
    ) -> int:
        """
        Set the state of a particle.

        TSUNAMI will rescale the positions and velocities
        to the center of mass frame.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        mass : float
            Particle mass.
        radius : float
            Particle radius.
        x, y, z : float
            Particle position.
        vx, vy, vz : float
            Particle velocity.
        wx, wy, wz : float
            Particle spin.

        Returns
        -------
        0 :
            State was set successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)

        self._mass[i] = mass
        self._radius[i] = radius
        self._pos[i, 0] = x
        self._pos[i, 1] = y
        self._pos[i, 2] = z
        self._vel[i, 0] = vx
        self._vel[i, 1] = vy
        self._vel[i, 2] = vz
        self._spin[i, 0] = wx
        self._spin[i, 1] = wy
        self._spin[i, 2] = wz

        self.tsunami.override_masses(self._mass)
        self.tsunami.override_position_and_velocities(self._pos, self._vel)

        return 0

    def get_mass(
        self,
        index_of_the_particle: int,
        mass: ValueHolder
    ) -> int:
        """
        Retrieve the mass of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        mass : ValueHolder[float]
            ValueHolder instance to return the mass value.

        Returns
        -------
        0 :
            Mass was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        mass.value = self._mass[i]
        return 0

    def set_mass(
        self,
        index_of_the_particle: int,
        mass: float
    ) -> int:
        """
        Set the mass of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        mass : float
            Particle mass.

        Returns
        -------
        0 :
            Mass was set successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)

        self._mass[i] = mass
        self.tsunami.override_masses(self._mass)

        return 0

    def get_radius(
        self,
        index_of_the_particle: int,
        radius: ValueHolder
    ) -> int:
        """
        Retrieve the radius of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        radius : ValueHolder[float]
            ValueHolder instance to return the radius value.

        Returns
        -------
        0 :
            Radius was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        radius.value = self._radius[i]
        return 0

    def get_position(
        self,
        index_of_the_particle: int,
        x: ValueHolder,
        y: ValueHolder,
        z: ValueHolder
    ) -> int:
        """
        Retrieve the position of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        x, y, z : ValueHolder[float]
            ValueHolder instance to return the position
            of the particle.

        Returns
        -------
        0 :
            Position was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        x.value = self._pos[i, 0]
        y.value = self._pos[i, 1]
        z.value = self._pos[i, 2]

        return 0

    def set_position(
        self, index_of_the_particle: int, x: float, y: float, z: float
    ) -> int:
        """
        Set the position of a particle. Ensure to update its mass first!

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        x, y, z : float
            Particle position.

        Returns
        -------
        0 :
            Position was set successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)

        self._pos[i, 0] = x
        self._pos[i, 1] = y
        self._pos[i, 2] = z

        self.tsunami.override_position_and_velocities(self._pos, self._vel)

        return 0

    def get_velocity(
        self,
        index_of_the_particle: int,
        vx: ValueHolder,
        vy: ValueHolder,
        vz: ValueHolder
    ) -> int:
        """
        Retrieve the velocity of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        vx, vy, vz : ValueHolder[float]
            ValueHolder instance to return the velocity
            of the particle.

        Returns
        -------
        0 :
            Velocity was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        vx.value = self._vel[i, 0]
        vy.value = self._vel[i, 1]
        vz.value = self._vel[i, 2]

        return 0

    def set_velocity(
        self, index_of_the_particle: int, vx: float, vy: float, vz: float
    ) -> int:
        """
        Set the velocity of a particle. Ensure to update its mass first!

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        vx, vy, vz : float
            Particle velocity.

        Returns
        -------
        0 :
            Velocity was set successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)

        self._vel[i, 0] = vx
        self._vel[i, 1] = vy
        self._vel[i, 2] = vz

        self.tsunami.override_position_and_velocities(self._pos, self._vel)

        return 0

    def get_spin(
        self,
        index_of_the_particle: int,
        wx: ValueHolder,
        wy: ValueHolder,
        wz: ValueHolder
    ) -> int:
        """
        Retrieve the spin of a particle.

        Parameters
        ----------
        index_of_the_particle : int
            Particle index as returned by `new_particle`.
        wx, wy, wz : ValueHolder[float]
            ValueHolder instance to return the spin
            of the particle.

        Returns
        -------
        0 :
            Spin was retrieved successfully.

        Raises
        ------
        ValueError :
            If `index_of_the_particle` is not valid.
        """
        i = self._get_particle_index_by_id(index_of_the_particle)
        self.synchronize_model()

        wx.value = self._spin[i, 0]
        wy.value = self._spin[i, 1]
        wz.value = self._spin[i, 2]

        return 0

    def set_units(self, Mscale: float, Lscale: float) -> int:
        """
        Set the mass and length scale of Tsunami.

        Parameters
        ----------
        Mscale : float
            Mass scale of Tsunami in MSun.
        Lscale : float
            Length scale of Tsunami in AU.

        Returns
        -------
        0 :
            Units were set successfully.
        """
        self.tsunami.set_units(float(Mscale), float(Lscale))
        return 0

    def get_Mscale(self, Mscale: ValueHolder) -> int:
        """
        Get mass unit of Tsunami in MSun.

        To set `Mscale`, use `set_units`.
        """
        Mscale.value = self.tsunami.Mscale
        return 0

    def get_Lscale(self, Lscale: ValueHolder) -> int:
        """
        Get length unit of Tsunami in AU.

        To set `Lscale`, use `set_units`.
        """
        Lscale.value = self.tsunami.Lscale
        return 0

    def get_Tscale(self, Tscale: ValueHolder) -> int:
        """
        Get time unit of Tsunami in years.
        Derived from `Mscale`, `Lscale`, and G=1.

        This value is read only.
        """
        Tscale.value = self.tsunami.Tscale
        return 0

    def get_Vscale(self, Vscale: ValueHolder) -> int:
        """
        Get velocity unit of Tsunami in km/s.
        Derived from `Mscale`, `Lscale`, and G=1.

        This value is read only.
        """
        Vscale.value = self.tsunami.Vscale
        return 0

    def get_speed_of_light(self, speed_of_light: ValueHolder) -> int:
        """
        Get Tsunami speed of light in N-body units.

        Parameters
        ----------
        speed_of_light : ValueHolder[float]
            ValueHolder instance to return the speed of light value.

        Returns
        -------
        0 :
            Speed of light was retrieved successfully.
        """
        speed_of_light.value = self.tsunami.speed_of_light
        return 0

    def get_equilibrium_tides(self, equilibrium_tides: ValueHolder) -> int:
        """
        Get equilibrium tides parameter. Enables or disables equilibrium tides.

        Parameters
        ----------
        equilibrium_tides : ValueHolder[bool]
            ValueHolder instance to return equilibrium tides value.
            If `True`, enables equilibrium tides.

        Returns
        -------
        0 :
            Equilibrium tides parameter was retrieved successfully.
        """
        equilibrium_tides.value = self.tsunami.Conf.wEqTides
        return 0

    def set_equilibrium_tides(self, equilibrium_tides: bool) -> int:
        """
        Set equilibrium tides parameter. Enables or disables equilibrium tides.

        Parameters
        ----------
        equilibrium_tides : ValueHolder[bool]
            ValueHolder instance to return equilibrium tides value.
            If `True`, enables equilibrium tides.

        Returns
        -------
        0 :
            Equilibrium tides parameter was set successfully.
        """
        self.tsunami.Conf.wEqTides = bool(equilibrium_tides)
        return 0

    def get_dynamical_tides(self, dynamical_tides: ValueHolder) -> int:
        """
        Get dynamical tides parameter. Enables or disables dynamical tides.

        Parameters
        ----------
        dynamical_tides : ValueHolder[bool]
            ValueHolder instance to return dynamical tides value.
            If `True`, enables dynamical tides.

        Returns
        -------
        0 :
            Dynamical tides parameter was retrieved successfully.
        """
        dynamical_tides.value = self.tsunami.Conf.wDynTides
        return 0

    def set_dynamical_tides(self, dynamical_tides: bool) -> int:
        """
        Set dynamical tides parameter. Enables or disables dynamical tides.

        Parameters
        ----------
        dynamical_tides : ValueHolder[bool]
            ValueHolder instance to return dynamical tides value.
            If `True`, enables dynamical tides.

        Returns
        -------
        0 :
            Dynamical tides parameter was set successfully.
        """
        self.tsunami.Conf.wDynTides = bool(dynamical_tides)
        return 0

    def initialize_tidal_parameters(
        self,
        kaps: NDArray[np.float64],
        taulag: NDArray[np.float64],
        polyt: NDArray[np.float64],
        gyrad: NDArray[np.float64],
        N: int
    ) -> int:
        """
        Initialize tidal parameters for all particles in TSUNAMI.

        Parameters
        ----------
        kaps : NDArray[np.float64]
            Apsidal motion constant (love number) for each particle
            in dimensionless units.
        taulag : NDArray[np.float64]
            Tidal lag time for each particle in units of time.
        polyt : NDArray[np.float64]
            Polytropic index for each particle in dimensionless units.
        gyrad : NDArray[np.float64]
            Gyration radius for each particle in dimensionless units.
        N : int
            Number of particles. Must match the length of `kaps`,
            `taulag`, `polyt`, and `gyrad`. Is used internally by AMUSE
            and is not set by the caller.

        Returns
        -------
        0 :
            Tidal parameters were initialized successfully.

        Raises
        ------
        ValueError :
            If `kaps`, `taulag`, `polyt`, or `gyrad` do not all have
            length `N`.
        """
        if any(len(arr) != N for arr in [kaps, taulag, polyt, gyrad]):
            raise ValueError(
                f'All tidal parameter arrays must have the same length ({N}) '
                f'with one value for each particle.'
            )
        self.tsunami.initialize_tidal_parameters(
            kaps, taulag, polyt, gyrad
        )
        self.synchronize_model()
        return 0

    def get_alpha(self, alpha: ValueHolder) -> int:
        """
        Get alpha regularization parameter.

        Parameters
        ----------
        alpha : ValueHolder[float]
            ValueHolder instance to return alpha value.

        Returns
        -------
        0 :
            Alpha parameter was retrieved successfully.
        """
        alpha.value = self.tsunami.Conf.alpha
        return 0

    def set_alpha(self, alpha: float) -> int:
        """
        Set alpha regularization parameter.

        Parameters
        ----------
        alpha : float
            alpha regularization value.

        Returns
        -------
        0 :
            Alpha parameter was set successfully.
        """
        self.tsunami.Conf.alpha = alpha
        return 0

    def get_beta(self, beta: ValueHolder) -> int:
        """
        Get beta regularization parameter.

        Parameters
        ----------
        beta : ValueHolder[float]
            ValueHolder instance to return beta value.

        Returns
        -------
        0 :
            Beta parameter was retrieved successfully.
        """
        beta.value = self.tsunami.Conf.beta
        return 0

    def set_beta(self, beta: float) -> int:
        """
        Set beta regularization parameter.

        Parameters
        ----------
        beta : float
            beta regularization value.

        Returns
        -------
        0 :
            Beta parameter was set successfully.
        """
        self.tsunami.Conf.beta = beta
        return 0

    def get_gamma(self, gamma: ValueHolder) -> int:
        """
        Get gamma regularization parameter.

        Parameters
        ----------
        gamma : ValueHolder[float]
            ValueHolder instance to return gamma value.

        Returns
        -------
        0 :
            Gamma parameter was retrieved successfully.
        """
        gamma.value = self.tsunami.Conf.gamma
        return 0

    def set_gamma(self, gamma: float) -> int:
        """
        Set gamma regularization parameter.

        Parameters
        ----------
        gamma : float
            gamma regularization value.

        Returns
        -------
        0 :
            Gamma parameter was set successfully.
        """
        self.tsunami.Conf.gamma = gamma
        return 0

    def get_pn(self, pn: ValueHolder) -> int:
        """
        Get pn parameter. Enables or disables
        post-Newtonian corrections.

        Parameters
        ----------
        pn : ValueHolder[bool]
            ValueHolder instance to return pn value. If `True`,
            enables post-Newtonian corrections. Ensure that `pn1`,
            `pn2`, `pn25`, `pn3` or `pn35` is also set.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn.value = self.tsunami.Conf.wPNs
        return 0

    def set_pn(self, pn: bool) -> int:
        """
        Set pn parameter. Enables or disables
        post-Newtonian corrections.

        Parameters
        ----------
        pn : bool
            If `True`, enables post-Newtonian corrections. Ensure
            that `pn1`, `pn2`, `pn25`, `pn3` or `pn35` is also set.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.wPNs = bool(pn)
        return 0

    def get_pn1(self, pn1: ValueHolder) -> int:
        """
        Get pn1 parameter. Enables or disables
        post-Newtonian corrections of order 1.

        Parameters
        ----------
        pn1 : ValueHolder[bool]
            ValueHolder instance to return pn1 value. If `True`,
            enables post-Newtonian corrections of order 1.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn1.value = self.tsunami.Conf.pn1
        return 0

    def set_pn1(self, pn1: bool) -> int:
        """
        Set pn1 parameter. Enables or disables
        post-Newtonian corrections of order 1.

        Parameters
        ----------
        pn1 : bool
            If True, enables post-Newtonian corrections of order 1.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.pn1 = bool(pn1)
        return 0

    def get_pn2(self, pn2: ValueHolder) -> int:
        """
        Get pn2 parameter. Enables or disables
        post-Newtonian corrections of order 2.

        Parameters
        ----------
        pn2 : ValueHolder[bool]
            ValueHolder instance to return pn2 value. If `True`,
            enables post-Newtonian corrections of order 2.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn2.value = self.tsunami.Conf.pn2
        return 0

    def set_pn2(self, pn2: bool) -> int:
        """
        Set pn2 parameter. Enables or disables
        post-Newtonian corrections of order 2.

        Parameters
        ----------
        pn2 : bool
            If True, enables post-Newtonian corrections of order 2.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.pn2 = bool(pn2)
        return 0

    def get_pn25(self, pn25: ValueHolder) -> int:
        """
        Get pn25 parameter. Enables or disables
        post-Newtonian corrections of order 2.5.

        Parameters
        ----------
        pn25 : ValueHolder[bool]
            ValueHolder instance to return pn25 value. If `True`,
            enables post-Newtonian corrections of order 2.5.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn25.value = self.tsunami.Conf.pn25
        return 0

    def set_pn25(self, pn25: bool) -> int:
        """
        Set pn25 parameter. Enables or disables
        post-Newtonian corrections of order 2.5.

        Parameters
        ----------
        pn25 : bool
            If True, enables post-Newtonian corrections of order 2.5.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.pn25 = bool(pn25)
        return 0

    def get_pn3(self, pn3: ValueHolder) -> int:
        """
        Get pn3 parameter. Enables or disables
        post-Newtonian corrections of order 3.

        Parameters
        ----------
        pn3 : ValueHolder[bool]
            ValueHolder instance to return pn3 value. If `True`,
            enables post-Newtonian corrections of order 3.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn3.value = self.tsunami.Conf.pn3
        return 0

    def set_pn3(self, pn3: bool) -> int:
        """
        Set pn3 parameter. Enables or disables
        post-Newtonian corrections of order 3.

        Parameters
        ----------
        pn3 : bool
            If True, enables post-Newtonian corrections of order 3.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.pn3 = bool(pn3)
        return 0

    def get_pn35(self, pn35: ValueHolder) -> int:
        """
        Get pn35 parameter. Enables or disables
        post-Newtonian corrections of order 3.5.

        Parameters
        ----------
        pn35 : ValueHolder[bool]
            ValueHolder instance to return pn35 value. If `True`,
            enables post-Newtonian corrections of order 3.5.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was retrieved successfully.
        """
        pn35.value = self.tsunami.Conf.pn35
        return 0

    def set_pn35(self, pn35: bool) -> int:
        """
        Set pn35 parameter. Enables or disables
        post-Newtonian corrections of order 3.5.

        Parameters
        ----------
        pn35 : bool
            If True, enables post-Newtonian corrections of order 3.5.

        Returns
        -------
        0 :
            Post-Newtonian correction parameter was set successfully.
        """
        self.tsunami.Conf.pn35 = bool(pn35)
        return 0

    def get_potential_energy(self, potential_energy: ValueHolder) -> int:
        """
        Get total potential energy of the system. This value is read only.

        Parameters
        ----------
        potential_energy : ValueHolder[float]
            ValueHolder instance to return the potential energy of the system.

        Returns
        -------
        0 :
            Potential energy was retrieved successfully.
        """
        potential_energy.value = self.tsunami.pot
        return 0

    def get_kinetic_energy(self, kinetic_energy: ValueHolder) -> int:
        """
        Get total kinetic energy of the system. This value is read only.

        Parameters
        ----------
        kinetic_energy : ValueHolder[float]
            ValueHolder instance to return the kinetic energy of the system.

        Returns
        -------
        0 :
            Kinetic energy was retrieved successfully.
        """
        kinetic_energy.value = self.tsunami.kin
        return 0

    def get_deltaE(self, deltaE: ValueHolder) -> int:
        """
        Get fractional energy error `deltaE = (E - E0) / E0` excluding
        energy losses due to dissipative terms (tides, post-Newtonians,
        external potentials)

        Parameters
        ----------
        deltaE : ValueHolder[float]
            ValueHolder instance to return the fraction energy error value.

        Returns
        -------
        0 :
            Energy error was retrieved successfully.
        """
        deltaE.value = self.tsunami.deltaE
        return 0

    def get_time(self, time: ValueHolder) -> int:
        """
        Get current model time in N-body units.

        Returns
        -------
        0 :
            Current model time was retrieved successfully.
        """
        time.value = self.tsunami.time
        return 0

    def get_number_of_particles(self, number_of_particles: ValueHolder) -> int:
        """
        Get total number of particles in TSUNAMI.

        Parameters
        ----------
        number_of_particles : ValueHolder[int]
            ValueHolder instance to return the number of particles in TSUNAMI.

        Returns
        -------
        0 :
            Number of particles was retrieved successfully.
        """
        number_of_particles.value = len(self._pos)
        return 0

    def _get_particle_index_by_id(self, index_of_the_particle: int) -> int:
        """
        Given an AMUSE particle id as returned by `new_particle`, find the
        index of that particle in TSUNAMI.

        This ensures that `index_of_the_particle` always points to the
        correct particle, even after particle deletion.

        Parameters
        ----------
        index_of_the_particle : int
            Particle id as returned by `new_particle`.
        """
        idx = np.where(self._ids == index_of_the_particle)[0]
        if idx.size == 0:
            raise ValueError(f'Particle id: {index_of_the_particle} not found!')

        return int(idx[0])

    def _get_new_id(self) -> int:
        """
        Create a unique, monotonically increasing id for `new_particle`.

        Returns
        -------
        int :
            New particle id to be assigned.
        """
        new_id = self._next_particle_id
        self._next_particle_id += 1

        return int(new_id)

    def _clear_temporary_particle_buffers(self) -> None:
        """Clear temporary particle buffers after committing particles."""
        self._id_list.clear()
        self._mass_list.clear()
        self._radius_list.clear()
        self._pos_list.clear()
        self._vel_list.clear()
        self._spin_list.clear()

class TsunamiInterface(
    PythonCodeInterface,
    GravitationalDynamicsInterface,
    GravityFieldInterface,
    LiteratureReferencesMixIn,
):
    """
    Tsunami is a direct few-body code with algorithmic regularization,
    tidal forces and post-Newtonian corrections.

    .. [#] ADS:2023IAUS..362..404T (Trani, A. A. & Spera, M., 2023, IAUS, 362, 404)
    """

    def __init__(self, **kwargs):
        PythonCodeInterface.__init__(
            self,
            TsunamiImplementation,
            'tsunami_worker',
            **kwargs
        )
        LiteratureReferencesMixIn.__init__(self)

    @remote_function
    def evolve_model_dtmax(time='d'):
        returns ()

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.OUT
        )
        for name in ['mass', 'radius', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'wx', 'wy', 'wz']:
            function.addParameter(name, dtype='d', direction=function.IN)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            particle was added to Tsunami
        """
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.IN
        )
        for name in ['mass', 'radius', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'wx', 'wy', 'wz']:
            function.addParameter(name, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            particle state was retrieved
        """
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='i', direction=function.IN
        )
        for name in ['mass', 'radius', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'wx', 'wy', 'wz']:
            function.addParameter(name, dtype='d', direction=function.IN)
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        """
        return function

    @remote_function(can_handle_array=True)
    def get_spin(index_of_the_particle='i'):
        returns (wx='d', wy='d', wz='d')

    @remote_function
    def get_speed_of_light():
        returns (speed_of_light='d')

    @remote_function
    def get_deltaE():
        returns (deltaE='d')

    @remote_function
    def get_equilibrium_tides():
        returns (equilibrium_tides='b')

    @remote_function
    def set_equilibrium_tides(equilibrium_tides='b'):
        returns ()

    @remote_function
    def get_dynamical_tides():
        returns (dynamical_tides='b')

    @remote_function
    def set_dynamical_tides(dynamical_tides='b'):
        returns ()

    @remote_function(must_handle_array=True)
    def initialize_tidal_parameters(
        kaps='d', taulag='d', polyt='d', gyrad='d'
    ):
        returns ()

    @remote_function
    def get_alpha():
        returns (alpha='d')

    @remote_function
    def set_alpha(alpha='d'):
        returns ()

    @remote_function
    def get_beta():
        returns (beta='d')

    @remote_function
    def set_beta(beta='d'):
        returns ()

    @remote_function
    def get_gamma():
        returns (gamma='d')

    @remote_function
    def set_gamma(gamma='d'):
        returns ()

    @remote_function
    def get_pn():
        returns (pn='b')

    @remote_function
    def set_pn(pn='b'):
        returns ()

    @remote_function
    def get_pn1():
        returns (pn1='b')

    @remote_function
    def set_pn1(pn1='b'):
        returns ()

    @remote_function
    def get_pn2():
        returns (pn2='b')

    @remote_function
    def set_pn2(pn2='b'):
        returns ()

    @remote_function
    def get_pn25():
        returns (pn25='b')

    @remote_function
    def set_pn25(pn25='b'):
        returns ()

    @remote_function
    def get_pn3():
        returns (pn3='b')

    @remote_function
    def set_pn3(pn3='b'):
        returns ()

    @remote_function
    def get_pn35():
        returns (pn35='b')

    @remote_function
    def set_pn35(pn35='b'):
        returns ()

    @remote_function
    def set_units(Mscale='d', Lscale='d'):
        returns ()

    @remote_function
    def get_Mscale():
        returns (Mscale='d')

    @remote_function
    def get_Lscale():
        returns (Lscale='d')

    @remote_function
    def get_Tscale():
        returns (Tscale='d')

    @remote_function
    def get_Vscale():
        returns (Vscale='d')


class Tsunami(GravitationalDynamics, GravityFieldCode, CommonCode):
    """
    Tsunami is a direct few-body code with algorithmic regularization,
    tidal forces and post-Newtonian corrections.
    """
    def __init__(self, convert_nbody=None, **options):
        """Create a Tsunami instance to run simulations with."""

        legacy_interface = TsunamiInterface(**options)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

        if convert_nbody is not None:
            M_s = convert_nbody.to_si(1 | generic_unit_system.mass).value_in(u.MSun)
            L_s = convert_nbody.to_si(1 | generic_unit_system.length).value_in(u.AU)
            self.set_units(M_s, L_s)

        Particles.add_global_vector_attribute('spin', ['wx', 'wy', 'wz'])


    def define_state(self, handler):
        """Define the state model of the code."""
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)


    def define_properties(self, handler):
        """Define properties of the code. These are read-only!"""
        GravitationalDynamics.define_properties(self, handler)
        handler.add_property('get_kinetic_energy', public_name='kinetic_energy')
        handler.add_property('get_potential_energy', public_name='potential_energy')
        handler.add_property('get_total_energy', public_name='total_energy')
        handler.add_property('get_deltaE', public_name='deltaE')
        handler.add_property('get_time', public_name='model_time')
        handler.add_property('get_Mscale', public_name='Mscale')
        handler.add_property('get_Lscale', public_name='Lscale')
        handler.add_property('get_Tscale', public_name='Tscale')
        handler.add_property('get_Vscale', public_name='Vscale')
        handler.add_property('get_speed_of_light', public_name='speed_of_light')


    def define_methods(self, handler):
        """Map legacy functions in TsunamiInterface into Tsunami user methods."""
        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            'evolve_model_dtmax',
            (ns.time,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'new_particle',
            (
                ns.mass,
                ns.length,
                ns.length, ns.length, ns.length,
                ns.speed, ns.speed, ns.speed,
                1 / ns.time, 1 / ns.time, 1 / ns.time,
            ),
            (handler.INDEX, handler.ERROR_CODE)
        )

        handler.add_method(
            'get_state',
            (handler.INDEX,),
            (
                ns.mass,
                ns.length,
                ns.length, ns.length, ns.length,
                ns.speed, ns.speed, ns.speed,
                1 / ns.time, 1 / ns.time, 1 / ns.time,
                handler.ERROR_CODE
            ),
        )

        handler.add_method(
            'set_state',
            (
                handler.INDEX,
                ns.mass,
                ns.length,
                ns.length, ns.length, ns.length,
                ns.speed, ns.speed, ns.speed,
                1 / ns.time, 1 / ns.time, 1 / ns.time,
            ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_velocity',
            (handler.INDEX,),
            (
                ns.speed,
                ns.speed,
                ns.speed,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            'set_velocity',
            (
                handler.INDEX,
                ns.speed,
                ns.speed,
                ns.speed,
            ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_spin',
            (handler.INDEX,),
            (
                1 / ns.time,
                1 / ns.time,
                1 / ns.time,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            'get_speed_of_light',
            (),
            (ns.speed, handler.ERROR_CODE)
        )

        handler.add_method(
            'initialize_tidal_parameters',
            (handler.NO_UNIT, ns.time, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )


    def define_parameters(self, handler):
        """
        Define model parameters.

        These have a native function for getting their value, another one for setting,
        and a name, description and default value. Functions with the appropriate names
        must be defined in the native wrapper code.
        """
        handler.add_method_parameter(
            'get_equilibrium_tides',
            'set_equilibrium_tides',
            'wEqTides',
            'enable equilibrium tides.',
            default_value=False
        )

        handler.add_method_parameter(
            'get_dynamical_tides',
            'set_dynamical_tides',
            'wDynTides',
            'enable dynamical tides.',
            default_value=False
        )

        handler.add_method_parameter(
            'get_alpha',
            'set_alpha',
            'alpha',
            'alpha regularization parameter.',
            default_value=1.0
        )

        handler.add_method_parameter(
            'get_beta',
            'set_beta',
            'beta',
            'beta regularization parameter.',
            default_value=0.0
        )

        handler.add_method_parameter(
            'get_gamma',
            'set_gamma',
            'gamma',
            'gamma regularization parameter.',
            default_value=0.0
        )

        handler.add_method_parameter(
            'get_pn',
            'set_pn',
            'wPNs',
            'enable post newtonian corrections',
            default_value=False
        )

        handler.add_method_parameter(
            'get_pn1',
            'set_pn1',
            'pn1',
            'enable post newtonian corrections of order 1',
            default_value=True
        )

        handler.add_method_parameter(
            'get_pn2',
            'set_pn2',
            'pn2',
            'enable post newtonian corrections of order 2',
            default_value=True
        )

        handler.add_method_parameter(
            'get_pn25',
            'set_pn25',
            'pn25',
            'enable post newtonian corrections of order 2.5',
            default_value=True
        )

        handler.add_method_parameter(
            'get_pn3',
            'set_pn3',
            'pn3',
            'enable post newtonian corrections of order 3',
            default_value=True
        )

        handler.add_method_parameter(
            'get_pn35',
            'set_pn35',
            'pn35',
            'enable post newtonian corrections of order 3.5',
            default_value=True
        )


    def define_particle_sets(self, handler):
        """Define any particle sets inside the model."""
        GravitationalDynamics.define_particle_sets(self, handler)
        handler.add_getter('particles', 'get_spin', names=('wx', 'wy', 'wz'))


    def define_converter(self, handler):
        """Handle unit conversion if an (optional) unit converter is specified."""
        if self.unit_converter is not None:
            handler.set_converter(
                self.unit_converter.as_converter_from_si_to_generic()
            )
