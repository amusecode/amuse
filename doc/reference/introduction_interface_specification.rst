============
Introduction
============

In this document we will specify the low-level interfaces of the comminity codes. 
These interfaces should follow a pattern. This pattern is described bellow. This
pattern can be used to keep the interfaces of the codes consistent
between codes of different physical domains and to help in learning about
the functionality of the interfaces. In later chapters we define the specific
interfaces for the modules of different physical domains. The actual interface of
a code should derive from these interfaces.

These interfaces only specify the low-level interfaces. The low-level interfaces
do not support unit handling, data conversion and attributes. This functionality
is handled by the next-level interfaces described in nextlevel_.

Data Types
----------
The exact size (in bytes) and type of the values sent between the python 
layer and the community codes is very important. In this document we will 
specify the type of each argument and return value. 
Currently AMUSE supports 4 types, these are described in the following table. 

========== ======== ================ ================
Type name  Size     Fortran          C               
..         bytes    type             type             
========== ======== ================ ================ 
int32      4        integer          long             
float64    8        double precision double
float32    4        real             float
string     n                         char *
========== ======== ================ ================ 

Function template
------------------
All functions in the interface should follow the same template.
Each function returns an error or status code. Results are returned through 
the function arguments. In C these arguments need to be pointers
to valid memory locations.

.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: example_function
        
The error codes all have the same general form. Zero stands for no error, a negative
value indicates some error happened, a positive value is returned when the function
ends in a special state.

0 - OK
    Function encountered no error or special state
<0 - ERROR
    Something went wrong in the execution of the function
>0 - STATE
    Function has encountered an expected special state. For example the code
    has detected a collision between two stars.
    

Function categories
-------------------

Parameters
~~~~~~~~~~
Codes can have a number of parameters. Some parameters are domain specific, these
parameters are found for all codes in a specific domain (for example the smoothing
length in gravitational dynamics). Other parameters are only defined and used by
a specific code. The domain specific parameters are defined on the domain specific
interfaces, when a code supports the parameter it should implement the
specified functions. Other parameters have to be accessed with functions following
the template of the :meth:`~ amuse.community.interface.example.ExampleInterface.get_example_parameter` 
and :meth:`~amuse.community.interface.example.ExampleInterface.set_example_parameter` functions. 

.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: get_example_parameter, set_example_parameter
   

A function used to access (set or get) a parameter may only retrieve
or update the value of a single parameter. Functions setting two or more
parameters in one go are not supported by the next-level interfaces. After all 
parameters have been set, the  :meth:`~amuse.community.interface.example.ExampleInterface.initialize_code`
function should be called, this gives the code the opportunity prepare the model.


.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: initialize_code

   
Object Management
~~~~~~~~~~~~~~~~~
Codes can work on particles or grids (stars, black holes or gas). The methods 
in the *Object Management* category define the functionality to create, remove 
and query the particles or gridpoints in the codes. 

When a code supports objects, the code is responsible for managing these objects. 
The code needs to assign a unique index to a particle so that the particle can be 
referred to in other function calls. This is a *major* difference with the 
``MUSE`` code. Where the user of the code was responsible for assigning unique ids
to the particles. This change makes the implementation of the code simpler and
allows the code to support creation of new objects during simulation. For example
a hydrocode can add or delete gridpoints during the evolution of the model.

    
Object state
~~~~~~~~~~~~
Particles in the same physical domain can have a well known, *minimal* state. 
For example, in the gravitational dynamics domain the state of a particle can be
defined by a location vector, velocity vector, a mass and a radius. The methods
in the *Object State* category provide a way to access this state in one function.
  

Object State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a particle can be transfered with the get_state
and set_state functions. Some codes may support other properties of a particle,
the code can define ``get_`` and ``set_`` functions for these properties. These 
functions follow the pattern defined in the *Parameters* category. The functions
must either get or set a scalar property (1 argument) or
a vector property (3 arguments).

Model evolution
~~~~~~~~~~~~~~~
The main function of a code is often evolving a model in time or solving a steady
state solution. The methods that control model evolution or start and stop the 
model calculations all belong to the *Model evolution* category. At this time, 
no pattern is defined for the functions in this category.

Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model
calculations. All the functions in this category follow 
the 'get_name' pattern. The state of code should not change during a function 
call to a function in this category. The functions must either get
a scalar property (1 argument) or a vector property (3 arguments).

Services
~~~~~~~~
Some codes can provide services for other codes in the same or other physical
domains. For example, gravitational dynamics code might provide a function to
calculate the gravity force at a point. The methods that provide these
services all belong to this category.
