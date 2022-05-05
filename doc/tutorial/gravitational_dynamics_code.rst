====================================
Adding a Gravitational Dynamics Code
====================================

In the previous tutorial we have seen how to create an AMUSE interface to
a C++ or a Fortran code from scratch. In this tutorial, we will expand on this
with a number of additional features. We will do this through the implementation
of a simple gravitational dynamics code.

As this tutorial deals only with features of the Python interface, we restrict our
legacy code example to C++. The concepts described can just as easily be applied to
an interface to a Fortran code. 

We will also assume that we are working with the same environment settings as in
the previous example. Thus, we can create the initial directory structure of our new code 'SimpleGrav' with the following terminal command:

.. code-block:: bash

    > $AMUSE_DIR/build.py --type=c --mode=dir SimpleGrav


The Legacy Code
---------------
For this example we will use a very simple forward Eulerian integration routine as our
legacy code. We simply have a function that takes the initial condition of a set of particles (in the form of seven dynamic arrays), an integer containing the number of particles, and double precision scalars containing the time step, smoothing length, and gravitational constant. It outputs the new particle positions and velocities by updating the input arrays.

.. literalinclude:: simplegrav/SimpleGrav.cc
    :language: C++


.. note::

    Just like in the previous tutorial, the example algorithm is by no means
    particularly good. In fact, forward Eulerian integration might be the worst 
    method of integrating gravitational systems due to the tendency of the energy
    to diverge. Again, we use it as a readable example.

Simply put this script in the **src** directory as **code.cc**, and again change
the Makefile line::

    CODEOBJS = test.o

to::

    CODEOBJS = test.o code.o


Interface Templates
-------------------
In order to promote uniformity among interfaces, and to make it easier to create
an interface, a number of types of codes have interface templates. These classes
have a number of functions common to that type of code already defined. They can be
found in the folder **src/amuse/community/interface** of the main AMUSE github
repository. 

Gravitational dynamics is one type of code having a template interface, defined
in **gd.py**. We can use this template by inheriting from it:

.. code-block:: python

    from amuse.community.interface.gd import GravitationalDynamics
    from amuse.community.interface.gd import GravitationalDynamicsInterface

    class SimpleGravInterface(CodeInterface, GravitationalDynamicsInterface):
        ...

    class SimpleGrav(GravitationalDynamics):
        ...

The legacy interface now contains all legacy functions defined in **gd.py**.
This includes adding and deleting particles, getting and setting particle properties
(mass, position, velocity, acceleration, and radius), getting and setting 
parameters (smoothing length and begin time), and a host of other functions. Our
definition of the ``SimpleGravInterface`` is then simply:

.. code-block:: python

    class SimpleGravInterface(CodeInterface, GravitationalDynamicsInterface):

        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="SimpleGrav_worker",
                **keyword_arguments)

Note that when compiling an interface inheriting from
``GravitationalDynamicsInterface``, all legacy functions defined there
(and in its parent class, ``CommonCodeInterface``) must be
defined in the C++ interface. For now, you can simply let them return an error 
value (-1 or -2).

The object oriented interface similarly contains definitions for methods, properties,
particle sets, parameters, and states (the latter we will discuss later). Its definition
is simply:

.. code-block:: python

    class SimpleGrav(GravitationalDynamics):

        def __init__(self, convert_nbody = None, **keyword_arguments):
            
            GravitationalDynamics.__init__(self,
                                           SimpleGravInterface(**keyword_arguments),
                                           convert_nbody,
                                           **keyword_arguments)


Note the ``convert_nbody`` keyword that appears in the initializer. This is used
to handle the conversion of units between the legacy code and the object-oriented
interface. The gravitational dynamics interface template assumes that the legacy
codes work in dimensionless N-body units, but by passing a ``convert_nbody``
instance the unit conversions can be handled automatically.

The only part left to write for a functional version is the C++ interface. As mentioned before,
this must contain the corresponding functions of all legacy functions. Also note
that ``GravitationalDynamicsInterface`` inherits from yet another interface,
``CommonCodeInterface``, which contains four more legacy functions to be defined.
These are very general functions, mostly concerned with code logistics. Their 
implementation, and whether they are even necessary, depends on the legacy code. If
they are not needed they can simply return 0 (in the case of our legacy code, only
``cleanup_code`` has actual functionality). They are often called automatically,
when needed, through the state model (discussed below). The four functions from
``CommonCode``, plus three similar functions from
``GravitationalDynamics``, are described in the table below:

=========================== ===================================
definition                  description
=========================== ===================================
initialize_code             any code that must be executed when the code is
                            initialized, but before parameters are set, particles
                            are added, a grid is initialized, etc. 
cleanup_code                any code that must be executed when the code is
                            stopped. This typically frees up memory allocated
                            by the code.
commit_parameters           any initialization that must be done that is dependent
                            on parameters, for example the allocation of memory
                            for a grid.
recommit_parameters         similar to commit_parameters, but after 
                            commit_parameters has been called once. This could
                            be a redefinition of a grid.
commit_particles            any initialization that must be done that is dependent
                            on particles, for example the construction of a tree
                            of some kind from a set of particles.
recommit_particles          similar to commit_particles, but after 
                            commit_particles has been called once. This could
                            the addition or removal of particles from a tree.
synchronize_model           evolve all particles to the same time, if they
                            are at different times.
=========================== ===================================

The following code contains definitions for all legacy functions, although some
non-essential functions only return error values:


.. literalinclude:: simplegrav/interface_1.cc
    :language: C++



Properties & Parameters
-----------------------
A property of an AMUSE code is a read-only quantity of the code that contains
information about the (often physical) state of the code. Examples include the
current model time and the potential and kinetic energies of the total system.

A parameter of an AMUSE code is a quantity that is used in the internal workings
of the code. Examples include the time step in most types of evolution codes and
the metallicity in stellar evolution codes. These are often, but not always,
writable. 

Both are defined in the object-oriented interface, in the ``define_properties``
and ``define_parameters`` functions. 

Properties are added as follows:

.. code-block:: python

    def define_properties(self, handler):
        handler.add_property("get_property", public_name="property_name")

Here, ``get_property`` is a legacy interface function that returns a single value.
``public_name`` is an optional argument that defines what the property is called;
if it is not given, the name of the parameter of the legacy function is used.

Properties are accessed directly as properties of the code object:

.. code-block:: python

    gravity = SimpleGrav()
    the_current_time = gravity.model_time

Parameters are added as follows:

.. code-block:: python

   def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_parameter",
            "set_parameter",
            "parameter",
            "what is this parameter?",
            default_value = 1. | some_unit,
            is_vector = BOOL,
            must_set_before_get = BOOL,
        )

Here, ``get_parameter`` and ``set_parameter`` are two legacy interface functions
that respectively have one output and one input. ``parameter`` is the name of the
parameter in the code. The fourth parameter is a documentation string describing
what the parameter signifies. Finally, there are three optional parameters. 
``default_value`` is the value the parameter is set to by default.
If ``is_vector == True``, the parameter is a vector of a fixed length. This
length is defined as the number of parameters of the getters and setters
(scalar parameters have only one parameter). 
If ``must_set_before_get == True``, the parameter must be set before it can be returned. 
In most cases, this can be False.


Note that there are methods to define other types of parameters:

=========================== ============================
function                    purpose
=========================== ============================
add_boolean_parameter       the parameter is a boolean,
                            instead of an integer/float.
add_caching_parameter       the parameter is only set
                            once ``commit_parameters`` is
                            called.
add_array_parameter         the parameter is an array;
                            in contrast with a normal
                            parameter with ``is_vector=True``,
                            the getters and setters have
                            only a single argument, which
                            must be able to handle arrays.
                            Between the setter and the name,
                            an additional function must be
                            passed that specifies the size
                            of the array.
=========================== ============================

``add_alias_parameter`` is also available to add an alias to an existing parameter.
It is of the following form: ``add_alias_parameter("alias_parameter", "original_parameter", 
"parameter documentation")``.


Parameters are accessed through the ``parameters`` property of the code object:

.. code-block:: python

    gravity = SimpleGrav()
    param = gravity.parameters.parameter




.. note::

    With properties, particles, and parameters, we have three ways of communicating
    data between the legacy code and AMUSE. With complex codes that are not
    necessarily similar to other codes in AMUSE, it might be nontrivial in which
    of these manners data is to be communicated. As a rule of thumb, any quantity 
    relating to individual particles should be communicated as particle properties, 
    even if it is conceptually closer to a parameter. To demonstrate the thin line
    between these, if a stellar evolution code has a single metallicity for all
    stars, it will be parameter, whereas if different stars can have different
    metallicities, it will be a particle property.

    Parameters and properties, on the other hand, both apply to the system as a whole.
    Of these, parameters should not be altered by the legacy code itself, whereas
    properties are read-only. 

In our SimpleGrav example we have introduced a number of parameters that the template
interface does not have. For one, we want the time step to be a settable parameter
(note that the definition in **gd.py** has ``None`` as a setter function). Additionally,
we want to be able to set a smoothing length, and for didactic reasons we include an
entirely new parameter, ``gravity_strength``, which functions as a scaling of the
magnitude of the gravitational force. 

We have to define new legacy functions to set the time step and get and set the gravity
strength (the others are defined in **gd.py**). The legacy interface then looks as follows:

.. code-block:: python

    class SimpleGravInterface(CodeInterface,
                       GravitationalDynamicsInterface):

        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="SimpleGrav_worker",
                **keyword_arguments)

        @legacy_function
        def get_grav_fac ():
            function = LegacyFunctionSpecification()
            function.addParameter('grav_fac', dtype='float64', direction=function.OUT)
            function.result_type = 'int32'
            return function

        @legacy_function
        def set_grav_fac ():
            function = LegacyFunctionSpecification()
            function.addParameter('grav_fac', dtype='float64', direction=function.IN)
            function.result_type = 'int32'
            return function

        @legacy_function
        def set_time_step ():
            function = LegacyFunctionSpecification()
            function.addParameter('time_step', dtype='float64', direction=function.IN)
            function.result_type = 'int32'
            return function

In the object-oriented interface, we have to overload the ``define_methods`` and the 
``define_parameters`` functions:

.. code-block:: python

    def define_methods (self, handler):

        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            "get_grav_fac",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE)
        )

        handler.add_method(
            "set_grav_fac",
            (handler.NO_UNIT),
            (handler.ERROR_CODE)
        )

        handler.add_method(
            "set_time_step",
            (nbody_system.time),
            (handler.ERROR_CODE)
        )

        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE)
        )

        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length),
            (handler.ERROR_CODE)
        )


    def define_parameters (self, handler):

        GravitationalDynamics.define_parameters(self, handler)

        handler.add_method_parameter(
            "get_grav_fac",
            "set_grav_fac",
            "gravity_strength",
            "constant factor by which G is multiplied",
            default_value = 1.
        )

        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "time_step",
            "constant integrator timestep",
            default_value = 0.01 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length


We also call the equivalent functions of **gd.py** so we still inherit
from the template interface. 




State Model
-----------
Many legacy codes have pieces of logistical code that must be executed between
physically relevant pieces of code. For example, in a gravitational tree code,
the tree must be constructed after all particles have been added, but before
the system is evolved, and not necessarily between evolve calls. In order to automate this, a state model can be defined
for the code. This defines what functions can be run in what state of the code,
but also how to transition between states, and what functions trigger that
transition. A particularly convenient function is that it allows a transition, 
and thus
the function associated with that, to be triggered automatically. This allows the
aforementioned tree code to automatically build its tree between the definition
of the particle set and the evolution of the system. 

The state of the code can be handled automatically by defining a state model that 
describes the states the code can be in as a graph with the nodes as the state and 
transitions mediated by interface calls. The state model of a code is defined in
the ``define_state`` interface function by calling the following methods on its 
handler argument: 

- ``set_initial_state(state)``  this defines the initial entry state, given as a string label (usually 'UNINITIALIZED').

- ``add_method(state, method_name)`` this defines that the given method is allowed in state ``state``. Again the state is a string, the string can be prepended with '!' to indicate a method is forbidden in the state. If state is not yet defined it will be added to the state model.

- ``add_transition(state1, state2, method_name, is_auto=True)`` this adds a transition between states state1 and state2 (and adds these states if not previously defined) which is triggered by a call to the interface function method_name. The is_auto argument determines whether this transition is allowed to  be triggered automatically. 

A method that is not mentioned in any add_method, is allowed in any state (and
doesn't trigger any transitions). If the state model detects that an interface call 
needs a state change it tries to hop to a state where the interface call is allowed
by calling transitions that are added with is_auto=True. (this is only possible if
they don't have non-default arguments)

The state model can be built up by the above methods. The state model of a code can be printed:

- ``interface.state_machine.to_table_string()`` or 

- ``interface.state_machine.to_plantuml_string()``. 

Note that function arguments in the above methods are strings! They are evaluated later to methods of the low level code interface class!
(so they must be defined there, they need not be remote function but can be ordinary class methods) 

Note that it can be convenient to use a number of sentinel methods which do not (by default) do anything, these are:

- before_get_parameter, before_set_parameter, before_set_interface_parameter, 

- before_new_set_instance, before_get_data_store_names.

(e.g. before_get_parameter is called before each parameter is retrieved.)

State models have been defined for our code's parent class, ``GravitationalDynamics``, and its grandparent class, ``CommonCode``,
and our code is simple enough that we do not need to expand this. Thus we do not need to define it.

To complete this example we will take a look at the state model of ``GravitationalDynamics``:

.. code-block:: python

    def define_state(self, handler): 
        common.CommonCode.define_state(self, handler)   
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        ...
        
        
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_method('UPDATE', 'new_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        handler.add_method('EVOLVED', 'evolve_model')
        handler.add_transition('EVOLVED','RUN', 'synchronize_model')
        ...

Note that the entry state is defined in ``CommonCode``, and we have shortened it
for readability.

Let's focus on what happens when we add a particle after evolving for some time.
We can see that ``evolve_model`` is only allowed in the ``EVOLVED`` state, which
means that that is the state after evolving. ``new_particle``, then, is allowed in
the ``EDIT`` and ``UPDATE`` states, but also when transitioning from ``RUN`` to 
``UPDATE``. From ``EVOLVED`` to ``RUN``, we can transition automatically with
``synchronize_model``. Thus, after evolving, if we add a particle, the function
``synchronize_model`` is automatically called, and we end up in the ``UPDATE`` 
state. If we call ``evolve_model`` again, we first go from ``UPDATE`` to ``RUN`` 
with an automatic call of ``recommit_particles``, and from there ``evolve_model``
leads us from ``RUN`` back to ``EVOLVED``.




Literature References
---------------------
When adding your code to AMUSE you of course want your work to be recognized.
AMUSE actively provides the references of every code used in an AMUSE script,
at the end of every run of the script. The references are defined in the Python
interface, as in the following code snippet:

.. code-block:: python

    class SimpleGravInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface):
        """
        SimpleGrav is a forward Euler code to dynamically evolve a Newtonian
        gravity N-body system. 

        .. [#] Alvarez, Etienne; Monthly Notices of the Astronomical and Astrophysical Review Letters Z, Vol. 42 (2020)
        .. [#] Adams, Douglas; Hitchhiker's Guide to the Galaxy (1979)
        """
        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="simplegrav_worker",
                **keyword_arguments)
            LiteratureReferencesMixIn.__init__(self)

Upon finishing a script using ``SimpleGrav`` we will get the following warning:

.. code-block::
    /home/.../amuse/src/amuse/support/literature.py:78: AmuseWarning: 

    You have used the following codes, which contain literature references:

	    "SimpleGravInterface"
		    Alvarez, Etienne; Monthly Notices of the Astronomical and Astrophysical Review Letters Z, Vol. 42 (2020)
		    Adams, Douglas; Hitchhiker's Guide to the Galaxy (1979)


	    "AMUSE"
		    Portegies Zwart, S. & McMillan, S.L.W., 2018, "Astrophysical Recipes: the art of AMUSE", 
                AAS IOP Astronomy publishing (411 pages) [2018araa.book.....P]
		    ** Portegies Zwart, S. et al., 2013, Multi-physics Simulations Using a Hierarchical Interchangeable 
                Software Interface, Computer Physics Communications 183, 456-468 [2013CoPhC.183..456P]
		    ** Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software Environment, 
                Astronomy and Astrophysics 557, 84 [2013A&A...557A..84P]
		    Portegies Zwart, S. et al., 2009, A multiphysics and multiscale software environment 
                for modeling astrophysical systems, *New Astronomy*, **Volume 14**, **Issue 4**, 369-378 [2009NewA...14..369P]

      warnings.warn(prefix + self.all_literature_references_string(), exceptions.AmuseWarning)

Note how the ``.. [#]`` denotes each entry for the literature list.
