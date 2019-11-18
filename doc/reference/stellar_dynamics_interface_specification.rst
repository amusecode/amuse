=====================================
Stellar Dynamics Interface Definition
=====================================

Introduction
~~~~~~~~~~~~

In this chapter we describe the common interface for all stellar dynamics codes.

Parameters
~~~~~~~~~~
Gravity dynamics codes have at least one specified parameter. Other parameters need
to be specified on a per code basis. All parameters have to be accessed with functions following
the template of the ``get_eps`` and ``set_eps`` functions. A parameter access function may only
retrieve or update the value of a single parameter. After all parameters have been set, the 
``initialize_code`` function should be called, this gives the code the opportunity prepare the
model.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_eps2, set_eps2, initialize_code
   

Object Management
~~~~~~~~~~~~~~~~~
Most gravitational dynamics codes work on particles (stars, black holes or gas). The following 
methods define the functionality to create, remove and query the particles in the code. *Currently 
the interface does not handle different types of particles*

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: new_particle, delete_particle, get_number_of_particles, get_index_of_first_particle, get_index_of_next_particle

    
Object state
~~~~~~~~~~~~
Particles in gravitational dynamics have a well known, *minimal* state. This state is is defined
by a location, velocity and mass and radius. The state can be retrieved and updated with 
the following functions.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_state, set_state
   

Object State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a particle can be transferred with the get_state and set_state functions. To
support other properties (like acceleration), the code can define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)


.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_mass, set_mass, get_position, set_position, set_acceleration, get_acceleration, get_potential


Model evolution
~~~~~~~~~~~~~~~
The gravitational dynamics codes evolve the properties of the particles in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: commit_particles, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. The following 
function can be used to query the exact model time, the total energies and colliding particles.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_time, get_kinetic_energy, get_potential_energy, get_indices_of_colliding_particles, get_center_of_mass_velocity, get_center_of_mass_position, get_total_mass, get_total_radius

Services
~~~~~~~~
Some Gravitational Dynamics codes can provide services for other codes. Currently calculating
the gravity at a given point is the only specified function.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_gravity_at_point
