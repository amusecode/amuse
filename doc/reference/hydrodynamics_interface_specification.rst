==================================
Hydrodynamics Interface Definition
==================================


Introduction
~~~~~~~~~~~~

In this chapter we describe the common interface for all 
hydrodynamics, grid based codes. For particle based, SPH codes see
the gravitational dynamics specifications.

Parameters
~~~~~~~~~~
All parameters have to be accessed with functions following
the template of the ``get_timestep`` and ``set_timestep`` functions.
A parameter access function may only retrieve or 
update the value of a single parameter. After all parameters have been set, the 
``commit_parameters`` function should be called,
this gives the code the opportunity prepare the model.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_timestep, get_timestep, commit_parameters
   

Grid Management
~~~~~~~~~~~~~~~~~
Most hydrodynamical codes work on grids or a hierarchy of grids. The following 
methods define the functionality to setup and query the grids. 

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_index_range_inclusive, set_boundary, setup_mesh, get_position_of_index,  get_index_of_position
    
Grid state
~~~~~~~~~~~~
Grid points in a hydrodynamics code have a well known, *minimal* state. This state is is defined
by a density, momentum density and energy density. The state can 
be retrieved and updated with the following functions.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_grid_state, get_grid_state
   

Grid State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a grid point can be transferred with the fill_grid_state and get_grid_state functions. To
support other properties (like pressure or MHD properties), the code can define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)


.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_density, set_grid_density, get_momentum_density, set_grid_momentum_density, get_energy_density, set_grid_energy_density


Model evolution
~~~~~~~~~~~~~~~
The hydrodynamics codes evolve the properties all grid cells in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: initialize_grid, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. The following 
function can be used to query the exact model time.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_time

External fields
~~~~~~~~~~~~~~~
Some Hydrodynamics codes support external acceleration or potential fields. The following functions
can be used to enter a gravitational potential field.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_has_external_gravitational_potential, get_has_external_gravitational_potential, get_index_range_for_potential, set_potential, get_potential

