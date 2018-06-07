==========================
Datamodel
==========================


Introduction
------------

Containers
~~~~~~~~~~
The AMUSE datamodel is based on objects stored in containers. The 
containers store all the relevant data of the objects. One can work 
with the data in a container as a whole or with an indiviual object 
in a container. When working with an individual object all data will 
be retrieved from and stored in the container of that object. *This 
is different from normal python objects and lists where the lists 
store references to the objects and the data is stored on the 
objects*.

Amuse containers::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(2)
    >>> stars.mass = [2, 3] | units.MSun
    >>> sun = stars[0]
    >>> sun.mass = 1 |  units.MSun
    >>> print stars.mass
    [1.0, 3.0] MSun
    >>> stars.mass = [0.5, 5] | units.MSun
    >>> print sun.mass
    0.5 MSun
    
Python lists::

    >>> from amuse.datamodel import Particle
    >>> sun = Particle(mass = 1 | units.MSun)
    >>> star = Particle(mass = 3 | units.MSun)
    >>> stars = [sun, star]
    >>> print stars.mass # cannot access attributes on the list
    Traceback (most recent call last):
    AttributeError: 'list' object has no attribute 'mass'
    

Set or grid
~~~~~~~~~~~
AMUSE supports two kinds of containers. One container kind 
implements a one dimensional set of objects. You can add and remove 
objects from the set and combine sets in different ways. This 
container kind is called **'Particles'** and the individual objects 
stored are called **'Particle'**. The other container kind 
implements a multidimensional grid of objects. You cannot change the 
shape of a grid and objects cannot be added or removed from a 
grid. This container kind is called **'Grid'** and individual 
objects stored in the container are called **'GridPoint'**.  

With particles::
    
    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(2)
    >>> stars.mass = [2, 3] | units.MSun
    >>> print stars.mass
    [2.0, 3.0] MSun
    
With grids::

    >>> from amuse.datamodel import Grid
    >>> from amuse.units import units
    >>> field = Grid(2,2)
    >>> field.density = [[1, 2], [3, 4]] | units.g / units.cm ** 3
    >>> point = field[1,0]
    >>> point.density = 5  | units.g / units.cm ** 3
    >>> print field.density
    [[ 1.  2.], [ 5.  4.]] g / (cm**3)
    

Memory, code or file
~~~~~~~~~~~~~~~~~~~~
The containers in AMUSE can save the data for individual objects in 
different kinds of **'Stores'**. AMUSE currently implements 3 kinds 
of stores. One kind of store saves all data in memory, this kind is 
used by default and provides the fastest way to work with the 
containers. Another kind of store saves and retrieves all data from 
a community code. It uses MPI messages to communicate this data. 
This kind is used by default for containers provided by a code and 
is the primary way by which you can interact with the data inside a 
running code. The last kind of container stores and retrieves data from
an HDF5 file. This kind is only used when saving or loading a container.

Varying attributes
~~~~~~~~~~~~~~~~~~
In memory containers can support a user defined set of attributes 
for the contained objects. You can define a new attribute for the 
container by assigning a value to the attribute. In code containers 
support a pre-defined, per code set of attributes for the contained 
objects. You cannot define a new attribute on these containers. Also 
as a code may not allow some attributes to be set individually, the 
container also cannot set some attributes individually. For example 
you often cannot set the X position of a particle, you must set the 
X, Y, and Z position in one go.

Adding attributes to a set or object::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(keys = [1,2])
    >>> # you can add an attribute by assigning to a name on a set
    >>> stars.mass = [1, 3] | units.MSun 
    >>> sun = stars[0]
    >>> print sun
    Particle(1, mass=1.0 MSun)
    >>> sun.radius = 1 | units.RSun
    >>> print sun
    Particle(1, mass=1.0 MSun, radius=1.0 RSun)

Objects not classes
~~~~~~~~~~~~~~~~~~~
Different kinds of particles or gridpoints are not modelled by 
implementing subclasses. Instead, different kinds of particles are 
defined ad-hoc, by variable naming (for example planets, stars or 
cloud_particles) and varying attributes (stars have position and 
mass, cloud_particles have position and density). This allows you to 
model your problem with the names and attributes that best fit your 
model. Unfortunately, this way of working does remove a level of 
self description in the system. To mitigate this problem the 
containers defined in community codes and in example scripts all 
follow the same conventions. We describe these conventions in a 
later section in this chapter.

Diffent attributes and names for different kinds::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> # stars and planets share some attributes (radius)
    >>> # but also have attributes that make sense only for
    >>> # the specific kind (population for planets, 
    >>> # luminosity for stars)
    >>> stars = Particles(keys = [1,2])
    >>> stars.luminosity = [1, 3] | units.LSun 
    >>> stars.radius = 1 | units.RSun
    >>> planets = Particles(keys = [3,4])
    >>> planets.radius = 6371  | units.km
    >>> planets.population = [7000000000, 0]
    
Identity
~~~~~~~~
Each object in a container has a unique identity, no two objects in 
a container can have the same identity. In a **'Particles'** 
container this identity is provided by a unique 64-bit key. In a 
**'Grid'** container this identity is provided by the n-dimensional 
index in the grid. Objects with the same identity can exists in 
multiple containers. These objects are considered as the same 
*conceptual* object in AMUSE. Different containers will provide 
different information about that object. For example the same *star* 
could live in a container in memory, in a container of a stellar 
evolution code and in a container of a stellar dynamics code. AMUSE 
provides several functions to link these objects and to transfer 
data between them.


Different sets store information on the same object::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(keys = [1,2])
    >>> stars.luminosity = [1, 3] | units.LSun 
    >>> bodies = Particles(keys = [1,2])
    >>> bodies.mass = [1, 3] | units.MSun 
    >>> print bodies[0] == stars[0] # the same 'conceptual' object in different sets
    True
    >>> print bodies[0] == stars[1] # not the same object
    False
    >>> print bodies[0]
    Particle(1, mass=1.0 MSun)
    >>> # finding the coresponding particle in the stars
    >>> print bodies[0].as_particle_in_set(stars)
    Particle(1, luminosity=1.0 LSun)
    
    
Ownership
~~~~~~~~~
Objects in a container are owned by that container, the container 
controls the data and the life-cycle of each object. 




Particle keys
-------------
All particles have a unique 64-bit key. This key is created using 
a random number generator. The chances of duplicate keys using 
64-bit integers are finite but very low. The chance of a duplicate 
key can be determined by a generalization of the birthday problem.

Duplicate keys::

    >>> # given n random integers drawn from a discrete uniform distribution 
    >>> # with range [1,highest_integer], what is the probability
    >>> # p(n;highest_integer) that at least two numbers are the same?
    >>> import math
    >>> number_of_bits = 64
    >>> highest_integer = 2**number_of_bits
    >>> number_of_particles = 1000000.0 # one million
    >>> probability = 1.0 - math.exp( (-number_of_particles * (number_of_particles - 1.0))/ (2.0* highest_integer) )
    >>> print probability
    2.71050268896e-08
    >>> # can also set the probablity and determine the set size
    >>> probability = 0.00001 # 0.001 percent
    >>> number_of_particles = math.sqrt(2 * highest_integer * math.log(1 / (1.0 - probability)))
    >>> print number_of_particles
    19207725.6894
    
If you use large sets or want to load a lot of simulations with 
different particles into a script the probablity of encountering a 
duplicate may be too high. You can check for duplicates in a set of 
particles by calling ``has_duplicates`` on a set. You can also change
the key generator to better match your problem.
    

Sets of particles
------------------
The AMUSE datamodel assumes all particles come in sets. The
data of a particle is stored in the set.

.. automodule:: amuse.datamodel

    .. autoclass:: AbstractParticleSet
        :members:
        
        .. automethod:: __add__
        .. automethod:: __sub__
        
    .. autoclass:: Particles
        :members:
        
    .. autoclass:: ParticlesSubset
        :members:
        
    .. autoclass:: ParticlesWithUnitsConverted
        :members:

    object
    ---------------

    .. autoclass:: Particle
        :members:
        
        .. automethod:: __add__
        .. automethod:: __sub__
    
    Methods to retreive physical properties of the particles set
    ------------------------------------------------------------
    
.. automodule:: amuse.datamodel.particle_attributes

    .. autofunction:: center_of_mass 
        :noindex:
    .. autofunction:: center_of_mass_velocity
        :noindex:
    .. autofunction:: kinetic_energy
        :noindex:
    .. autofunction:: potential_energy 
        :noindex:
    .. autofunction:: particle_specific_kinetic_energy
        :noindex:
    .. autofunction:: particle_potential
        :noindex:

    
Implementation
--------------


.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];
        
        "AttributeStorage" -> "AbstractSet";
        "AbstractSet" -> "AbstractParticlesSet";
        "AbstractSet" -> "AbstractGrid";
        
        "AbstractParticlesSet" -> "Particles";
        "AbstractParticlesSet" -> "ParticlesSubSet";
        "AbstractParticlesSet" -> "ParticlesSuperSet";
        "AbstractGrid" -> "Grid";
        "AbstractGrid" -> "SubGrid";
        
        "AttributeStorage" -> "InMemoryStorage" [label = "implements"];
        "AttributeStorage" -> "InCodeStorage";
    }



.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];
        edge [color="dodgerblue2", fontcolor="dodgerblue2"];

        "Particles" -> "AttributeStorage" [headlabel=" 1", taillabel="1",label = "store"];
        "Grid" -> "AttributeStorage" [headlabel=" 1", taillabel="1", label = "store"];
        "SubGrid" -> "Grid" [headlabel=" 1", taillabel="1", label="view on", color="dodgerblue2"];
        "ParticlesSubSet" -> "Particles" [headlabel=" 1", taillabel="1", label="view on", color="dodgerblue2"];
        "ParticlesSuperSet" -> "Particles" [headlabel=" *", taillabel="1", label="view on", color="dodgerblue2"];
    }


.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];

        "AbstractSet" -> "AttributeStorage" [headlabel=" 1", taillabel="1", label = "stored_attributes", color="dodgerblue2" , fontcolor="dodgerblue2"];
        "AbstractSet" -> "DerivedAttribute" [headlabel=" *", taillabel="1", label = "derived_attributes", color="dodgerblue2", fontcolor="dodgerblue2"]; 
        "DerivedAttribute" -> "CalculatedAttribue";
        "DerivedAttribute" -> "VectorAttribute";
        "DerivedAttribute" -> "FunctionAttribute";
        "DerivedAttribute" -> "DomainAttribute";
    }

.. image:: derived_attributes.png


.. autoclass:: amuse.datamodel.AttributeStorage
    :members:



