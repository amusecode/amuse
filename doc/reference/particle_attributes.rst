===================
Particle Attributes
===================

Particle attributes are defined in the 
:py:mod:`amuse.datamodel.particle_attributes`. These attributes
can be accessed in two ways, on the particle(s) or as a function
imported from the module. When accessed on the particles, the
first parameter (usually called ```particles```) must not
be given.

To access these functions directly do:

.. code-block:: python

    from amuse.datamodel.particle_attributes import *
    from amuse.lab import new_plummer_sphere
    
    particles = new_plummer_sphere(100)
    print kinetic_enery(particles)

To access these functions as an attribute do:

.. code-block:: python

    from amuse.lab import new_plummer_sphere
    
    particles = new_plummer_sphere(100)
    
    print particles.kinetic_energy()
    

.. automodule:: amuse.datamodel.particle_attributes
    :members:
