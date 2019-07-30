=======================================
Stellar Evolution Interface Definition
=======================================

Introduction
~~~~~~~~~~~~
In this chapter we describe the common interface for stellar evolution codes.
Currently the interface for stellar evolutions codes that store state of
a star is specified. The Stellar Evolution codes that get all the needed
state from the function interface (are stateless), are not yet described.     

Parameters
~~~~~~~~~~
Stellar Evolution codes have at least one specified parameter. Other parameters need
to be specified on a per code basis. All parameters have to be accessed with functions following
the template of the ``get_metallicity`` and ``set_metallicity`` functions. A parameter access function may only
retrieve or update the value of a single parameter. After all parameters have been set, the 
``initialize_code`` function should be called, this gives the code the opportunity prepare the
model.

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: get_metallicity, set_metallicity, initialize_code
   

Object Management
~~~~~~~~~~~~~~~~~
A number of stellar evolution codes work on star objects. The following 
methods define the functionality to create, remove and query the particles in the code. 
*Currently the interface does not specify query function for stellar evolution, see stellar evolution for possible direction*

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: new_particle, delete_star

Object State
~~~~~~~~~~~~~~
To support properties (like acceleration), the code must define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)
*Currently only get functions are specified*

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: get_mass, get_radius, get_luminosity, get_temperature, get_age, get_stellar_type, get_stellar_type


Model evolution
~~~~~~~~~~~~~~~
The stellar evolution codes evolve the properties of the star in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: commit_particles, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. 
*Currently no specific stellar evolution diagnostics functions have been defined*


Services
~~~~~~~~
Some stellar evolution codes can provide services for other codes. 
*Currently no specific stellar evolution service functions have been defined*
