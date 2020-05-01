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
For this example we will use a very simple Eulerian integration routine as our
legacy code. We simply have a function that takes the initial condition of a set of particles (in the form of seven dynamics arrays), an integer containing the number of particles, and double precision scalars containing the time step, smoothing length, and gravitational constant. It outputs the new particle positions and velocities by updating the input arrays.

.. literalinclude:: simplegrav/code.cc
    :language: c++


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

    from amuse.community import *
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
definition of the **SimpleGravInterface** is then simply:

.. code-block:: python

    class SimpleGravInterface(CodeInterface, GravitationalDynamicsInterface):

        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="SimpleGrav_worker",
                **keyword_arguments)

Note that when compiling an interface inheriting from
**GravitationalDynamicsInterface**, all legacy functions defined there must be
defined in the C++ interface. For now, you can let them return error values.

The object oriented interface contains definitions for methods, properties,
particle sets, parameters, and states (the latter we will discuss later).
We can define these by calling the functions that define them in the
gravitational dynamics interface:

.. code-block:: python

    class SimpleGrav(GravitationalDynamics):

        def __init__(self, convert_nbody = None, **keyword_arguments):
            
            GravitationalDynamics.__init__(self,
                                           SimpleGravInterface(**keyword_arguments),
                                           convert_nbody,
                                           **keyword_arguments)


        def define_state (self, handler):

            GravitationalDynamics.define_state(self, handler)


        def define_methods (self, handler):

            GravitationalDynamics.define_methods(self, handler)


        def define_parameters (self, handler):

            GravitationalDynamics.define_parameters(self, handler)


        def define_particle_sets (self, handler):

            GravitationalDynamics.define_particle_sets(self, handler)


        def define_properties(self, handler):

            GravitationalDynamics.define_properties(self, handler)

Note the **convert_nbody** keyword that appears in the initializer. This is used
to handle the conversion of units between the legacy code and the object-oriented
interface. The gravitational dynamics interface template assumes that the legacy
codes work in dimensionless N-body units, but by passing a **convert_nbody**
instance the unit conversions can be handled automatically.

The only part left to write now is the C++ interface. As mentioned before,
this must contain the corresponding functions of all legacy functions. Also note
that **GravitationalDynamicsInterface** inherits from yet another interface,
**CommonCodeInterface**, which contains four more legacy functions to be defined.
These are very general functions, mostly concerned with code logistics. Their 
implementation, and whether they are even necessary, depends on the legacy code. If
they are not needed they can simply return 0 (in the case of our legacy code, only
**cleanup_code** has actual functionality). They are often called automatically,
when needed, through the state model (discussed below). The four functions from
**CommonCode**, plus three similar functions from
**GravitationalDynamicsInterface**, are described in the table below:

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

.. literalinclude:: simplegrav/interface.cc
    :language: c++


State Model
-----------
Many legacy codes have pieces of logistical code that must be executed between
physically relevant pieces of code. For example, in a gravitational tree code,
the tree must be constructed after all particles have been added, but before
the system is evolved. In order to automate this, a state model can be defined
for the code. This defines what functions can be run in what state of the code,
but also how to transition between states, and what functions trigger that
transition. A particularly convenient function is to allow a transition, and thus
the function associated with that, to be triggered automatically. This allows the
aforementioned tree code to automatically build its tree between the definition
of the particle set and the evolution of the system. 

The state of the code can be handled automatically by defining a state model that 
describes the states the code can be in as a graph with the nodes as the state and 
transitions mediated by interface calls. The state model of a code is defined in
the **define_state** interface function by calling the following methods on its 
handler argument: 

- **set_initial_state(state)**  this defines the initial entry state, given as a
string label (usually 'UNINITIALIZED').
- **add_method(state, method_name)** this defines that the given method is allowed
in state **state**. Again the state is a string, the string can be prepended with
'!' to indicate a method is forbidden in the state. If state is not yet defined it
will be added to the state model.
- **add_transition(state1, state2, method_name, is_auto=True)** this adds a
transition between states state1 and state2 (and adds these states if not
previously defined) which is triggered by a call to the interface function
method_name. The is_auto argument determines whether this transition is allowed to 
be triggered automatically. 

A method that is not mentioned in any add_method, is allowed in any state (and
doesn't trigger any transitions. If the state model detects that an interface call 
needs a state change it tries to hop to a state where the interface call is allowed
by calling transitions that are added with is_auto=True. (this is only possible if
they don't have non-default arguments)

The state model can be built up by the above methods. The state model of a code can be printed:
- **interface.state_machine.to_table_string()** or 
- **interface.state_machine.to_plantuml_string()**. 

Note that function arguments in the above methods are strings! They are evaluated later to methods of the low level code interface class!
(so they must be defined there, they need not be remote function but can be ordinary class methods) 

Note that it can be convenient to use a number of sentinel methods which do not (by default) do anything, these are:
- before_get_parameter, before_set_parameter, before_set_interface_parameter, 
- before_new_set_instance, before_get_data_store_names.
(e.g. before_get_parameter is called before each parameter is retrieved.)

State models have been defined for our code's parent class, **GravitationalDynamics**, and its grandparent class, **CommonCode**,
and our code is simple enough that we do not need to expand this. Thus we can again simply define our states to be those of the parent:

.. code-block:: python

    class SimpleGrav(GravitationalDynamics):

        def __init__(self, **options):
            GravitationalDynamics.__init__(self, SimpleGravInterface(**options), convert_nbody, **options)

        def define_methods(self, handler):
            ...

        def define_particle_sets(self, handler):
            ...

        def define_state(self, handler):
            GravitationalDynamics.define_state(self, handler)


To complete this example we will take a look at the state model of **GravitationalDynamics**:

.. code-block:: python

    def define_state(self, handler): 
        common.CommonCode.define_state(self, handler)   
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        handler.add_method('EVOLVED','before_get_parameter')
        
        
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
        handler.add_method('RUN', 'synchronize_model')
        handler.add_method('RUN', 'get_state')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_velocity')
        handler.add_method('RUN', 'get_potential')
        handler.add_method('RUN', 'get_potential_energy')
        handler.add_method('RUN', 'get_kinetic_energy')
        handler.add_transition('RUN', 'UPDATE', 'set_mass', False)
        handler.add_transition('RUN', 'UPDATE', 'set_position', False)
        handler.add_transition('RUN', 'UPDATE', 'set_velocity', False)
        handler.add_transition('RUN', 'UPDATE', 'set_radius', False)

Note that the entry state is defined in **CommonCode**.

Let's focus on what happens when we add a particle after evolving for some time. We can see that 
**evolve_model** changes the state to **EVOLVED**. The function **new_particle**, then, can only
be run in the **EDIT** or **UPDATE** states. Luckily, we can see that there is a path from
**EVOLVED** to **RUN**, and then from **RUN** to **UPDATE**, by calling **synchronize_model** and
**new_particle**, respectively. **synchronize_model** is just a function that makes sure all particles are
evolved up to the same time, which our example code does by construction. We have defined it, 
but the only thing it does is to **return 0**. If we then want to evolve
the system further, we need to go back to the **RUN** state, which we can access from **UPDATE** through
the **recommit_particles** function.


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

Upon finishing a script using **SimpleGrav** we will get the following warning:

.. code-block::
    /home/martijn/venvs/amuse_dev/amuse/src/amuse/support/literature.py:78: AmuseWarning: 

    You have used the following codes, which contain literature references:

	    "SimpleGravInterface"
		    Alvarez, Etienne; Monthly Notices of the Astronomical and Astrophysical Review Letters Z, Vol. 42 (2020)
		    Adams, Douglas; Hitchhiker's Guide to the Galaxy (1979)


	    "AMUSE"
		    Portegies Zwart, S. & McMillan, S.L.W., 2018, "Astrophysical Recipes: the art of AMUSE", AAS IOP Astronomy publishing (411 pages) [2018araa.book.....P]
		    ** Portegies Zwart, S. et al., 2013, Multi-physics Simulations Using a Hierarchical Interchangeable Software Interface, Computer Physics Communications 183, 456-468 [2013CoPhC.183..456P]
		    ** Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software Environment, Astronomy and Astrophysics 557, 84 [2013A&A...557A..84P]
		    Portegies Zwart, S. et al., 2009, A multiphysics and multiscale software environment for modeling astrophysical systems, *New Astronomy*, **Volume 14**, **Issue 4**, 369-378 [2009NewA...14..369P]

      warnings.warn(prefix + self.all_literature_references_string(), exceptions.AmuseWarning)

Note how the ``.. [#]`` denotes each entry for the literature list.
