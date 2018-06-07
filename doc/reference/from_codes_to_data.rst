==================
From Codes to Data
==================

Introduction
------------

The framework containts two distinct levels on wich interaction with the codes takes place.
The first level interacts directly with the code and is implemented as a set of functions
on a class. The second level is built on top of the first level and abstracts the function calls
to handling objects and datasets. Often multiple function-calls on the first level can be 
abstracted to a single statement (assignment or query) on the second level.

The first level
---------------
The first level is a direct interface to the code, in this chapter the kind of functions
supported by the first level interface code will be breifly described. All functions that can
be defined on the first level fall in two categories, those that handle scalars (a single
value for each parameter) or those that handle 1-D vectors (a list of values for each
parameters). For functions that handle 1-D vectors each vector must be of the same length. 

.. note:: 

    Not every function will fit in the two categories, but it is usually possible
    to rewrite a function or create a interfacing function that do fit into
    one of the two categories. Supporting only these two categories keeps the 
    communication layer simpler and allows for some optimizations in the 
    communication between python and C/Fortran codes.

An exampe of using the first level with scalars and vectors::

    from amuse.community.codes.athena.interface import AthenaInterface
    
    # create an instance of the code (will start an application 
    # in the background to handle all requests)
    hydro = AthenaInterface()
    
    # set parameters needed by the code
    # these are functions handling one scalar input 
    # parameter
    hydro.set_gamma(1.6666666666666667)
    hydro.set_courant_friedrichs_lewy_number(0.8)
    
    # define a grid having 5 cells in all directions and 
    # with the total length of the grid in each direction
    # is 1.0
    # this is a function handling multiple
    # scalar input parameters
    hydro.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
    
    # setup boundary conditions 
    # (can be periodic, reflective, outflow)
    hydro.set_boundary(
        "periodic","periodic",
        "periodic","periodic",
        "periodic","periodic"
    )
    
    # let the code do some work
    # (athena will allocate the grid)
    # this is a function handling no 
    # scalar input parameters and having 
    # 1 scalar output parameter
    hydro.commit_parameters()
    
    # lets print the center position
    # of one grid point
    print hydro.get_position_of_index(1,2,3)
    
    # all calls so far have been to functions handling scalar values
    # the next calls will be to functions handling vectors of values
    
    # lets print the center positions
    # of all grid points on one line
    print hydro.get_position_of_index(range(0,5), [0] * 5, [0] * 5)

In the previous example we used functions with scalar parameters and
vector parameters. The functions handling vectors often can also
handle scalars, the framework will take care of the necessary 
conversions.

All first level functions are not actual python functions, these
functions are instances of a special Python class that
implements function call handling. To continue our example::

    # let's take a look at the kind of functions
    # on the first level
    print hydro.get_position_of_index
    
    # you can ask the specificition of a
    # first level function
    print hydro.get_position_of_index.specification


Adding units
~~~~~~~~~~~~
The first step after defining a first level function is to
specify the units of the in- and output-parameters of the first
level function. 

    
