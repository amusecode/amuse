====================
Integrate a C++ code
====================


In this tutorial we will create an AMUSE interface to a C++
code. We will first define the legacy interface, then implement the 
code and finally build an object oriented interface on top of the 
legacy interface.

The legacy code will be a very simple and naive implementation to find
3 nearest neighbors of a particle.

The legacy code interface supports methods that transfer values to 
and from the code. The values do not have any units and no error 
handling is provided by the interface. We can add error handling, unit 
handling and more functions to the legacy interface by defining a 
subclass of a InCodeComponentImplementation (this is the objected oriented interface).

.. graphviz::

   digraph layers4 {
      fontsize=10.0;
        rankdir="LR";
        node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow];
        
        "Legacy Code" -> "Legacy Interface" -> "Object Oriented Interface" -> "Script";
    }

The legacy code in this tutorial will be a very simple and naive
implementation to find 3 nearest neighbors of a particle.

Two paths
---------
When defining the interface will walk 2 paths:

1. Management of particles in AMUSE (python)
2. Management of particles in the code (C or Fortran)

The first path makes sense for legacy codes that perform a 
transformation on the particles, or analyse the particles state or 
do not store any internal state between function calls (all data is 
external). For every function of the code, data of every particle is 
send to the code. If we expect multiple calls, the code would incur a 
high communication overhead and we are better of choosing path 2.

The second path makes sense for codes that already have management 
of a particles (or grid) or were we want to call multiple functions 
of the code and need to send the complete model to code for every 
function call. The code is first given the data, then calls are made 
to the code to evolve it's model or perform reduction steps on the 
data, finally the updated data is retrieved from the code. 

Procedure
---------

The suggested procedure for creating a new interface is as follows:

0. **Legacy Interface.** Start with creating the legacy 
   interface. Define functions on the interface to input and
   output relevant data.
   The InCodeComponentImplementation code depends on the legacy interface code.   
1. **Make a Class.** Create a subclass of the InCodeComponentImplementation class
2. **Define methods.** In the legacy interface we have defined functions
   with parameters. In the code interface we need to define the
   units of the parameters and if a parameter or return value
   is used as an errorcode.
3. **Define properties.** Some functions in the legacy interface can
   be better described as a property of the code. These are read only 
   variables, like the current model time.
4. **Define parameters.** Some functions in the legacy interface provide
   access to parameters of the code. Units and default values
   need to be defined for the parameters in this step
5. **Define sets or grids.** A code usually handles objects or gridpoints with
   attributes. In this step a generic interface is defined for these
   objects so that the interoperability between codes increases.

Before we start
---------------

This tutorial assumes you have a working amuse environment. Please 
ensure that amuse is setup correctly by running 'nosetests' in the 
amuse directory.


Environment variables
~~~~~~~~~~~~~~~~~~~~~
To simplify the work in the coming sections, we first define the 
environment variable 'AMUSE_DIR'. This environment variable must 
point to the root directory of AMUSE (this is the directory 
containing the build.py script).

.. code-block:: bash

    > export AMUSE_DIR=<path to the amuse root directory>
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR <path to the amuse root directory>

After building the code, we want to run and test the code. Check if 
amuse is available in your python path by running the following code 
on the command line.

.. code-block:: bash

    > python -c "import amuse"
    Traceback (most recent call last):
    File "<string>", line 1, in <module>
    ImportError: No module named amuse
    
If this code ends in a "ImportError" as shown in the example, the 
PYTHONPATH environment variable must be extended with the src directory
in AMUSE_DIR. 
We can do so by using one of the following commands.

.. code-block:: bash

    > export PYTHONPATH=${PYTHONPATH}:${AMUSE_DIR}/src
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR ${PYTHONPATH}:${AMUSE_DIR}/src
    
    
The name of our project
~~~~~~~~~~~~~~~~~~~~~~~
We will be writing a code to find the nearest neighbors of a particle, 
so let's call our project 'NearestNeighbor'.

Creating the initial directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First we need to create a directory for our project and put some 
files in it to help build the code. The fastest method to setup the 
directory is by using the build.py script.

.. code-block:: bash

    > $AMUSE_DIR/build.py --type=c --mode=dir NearestNeighbor

The script will generate a directory with all the files needed to 
start our project. It also generates a very small example legacy code 
with only one function ```echo_int```. We can build and test our new 
module::

    > cd nearestneighbor/
    > make all
    > $AMUSE_DIR/amuse.sh -c 'from interface import NearestNeighbor; print NearestNeighbor().echo_int(10)' 
    OrderedDictionary({'int_out':10, '__result':0})
    > nosetests -v
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.556s

    OK
    
    
.. note::

    The build.py script can be used to generate a range of files. To
    see what this file can do you can run the script with a ```--help```
    parameter, like so::

        > $AMUSE_DIR/build.py --help

The Legacy Code
---------------
Normally the legacy code already exists and our task is limited to 
defining and implementing an interface so that AMUSE scripts can 
access the code. For this tutorial we will implement our legacy code.

When a legacy code is integrated all interface code is put in one 
directory and all the legacy code is put in a **src** directory 
placed under this directory. The build.py script created a **src** 
directory for us, and we will put the nearest neighbor algorithm in 
this directory.

Go to the **src** directory and create a **code.cc** file, open 
this file in your favorite editor and copy and paste this code into 
it:
    
.. literalinclude:: nearestneighbor/code.cc
    :language: c++
    

.. note::

    This algorithm is un-optimized and has N*N order. It is not meant
    as very efficient code but as a readable example.

Before we can continue we also need to alter the **Makefile** in the 
**src** directory, so that our **code.cc** file is included in the 
build. To do so, open an editor on the Makefile and change the line::

    CODEOBJS = test.o

to::

    CODEOBJS = test.o code.o
    

Test if the code builds. As we have not coupled our algorithm to the 
interface we (we have not even defined an interface) we do not have 
any new functionality. In the legacy interface directory (not the 
**src** directory) do:

.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works, if the test fails for any reason please check that the C++
code is correct and that ``worker_code`` exists in your directory.

Path 1
------

Defining the legacy interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will first define a legacy interface so that we can call the 
**find_nearest_neighbors** function from python. AMUSE can interact 
with 2 classes of functions:


1. A function with all scalar input and output variables. All variables
   are simple, non-composite variables (like int or double).
   For example:
   
   .. code-block:: c++
    
        int example1(double input, double *output)
        {
            *output = input;
            return 0;
        }
        
2. A function with all vector (or list) input and output variables and
   a length variable. The return value is a scalar value.
   For example:
   
   .. code-block:: c++
    
        int example2(double *input, double *output, int n)
        {
            for(int i = 0; i < n; i++)
            {
                output[i] = input[i];
            }
            return 0;
        }

If you have functions that don't follow this pattern you need to 
define a convert function in C++ that provides an interface 
following one of the two patterns supported by AMUSE.

In our case the **find_nearest_neighbors** complies to pattern 2 and
we do not have to write any code in C++ to convert the function
to a compliant interface. We only have to specify the function in 
python. We do so by adding a **find_nearest_neighbors** method to
the **NearestNeighborInterface** class in the **interface.py** file.
Open and editor on the interfaces.py file and add the following method
to the NearestNeighborInterface class:

.. code-block:: python

    class NearestNeighborInterface(InCodeComponentImplementation):    
        #...
        
        @legacy_function
        def find_nearest_neighbors():
            function = LegacyFunctionSpecification()  
            function.must_handle_array = True 
            function.addParameter(
                'npoints', 
                dtype='int32', 
                direction=function.LENGTH)
            function.addParameter(
                'x', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'y', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'z', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'n0', 
                dtype='int32', 
                direction=function.OUT)
            function.addParameter(\
                'n1', 
                dtype='int32', 
                direction=function.OUT)
            function.addParameter(
                'n2', 
                dtype='int32', 
                direction=function.OUT)
            function.result_type = 'int32'
            return function
            
In the *find_nearest_neighbors* method we specify every parameter of 
the C++ function and the result type. For each parameter we need 
to define a name, data type and whether we will input, output (or 
input and output) data using this parameter. AMUSE knows only a 
limited amount of data types for parameters: float64, float32, int32, \
string, bool and int64. We also have a special parameter, with LENGTH as 
direction. This parameter is needed for all functions that follow 
pattern 2, it will be filled with the length of the input arrays. We 
also must specify that the function follows pattern 2 by setting 
```function.must_handle_array = True```.

Save the file and recompile the code.

.. code-block:: bash

    > make clean
    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works! But, how do we know the find_nearest_neighbors method really works?
Let's write a test and find out. Open an editor on the *test_nearestneighbor.py*
file and add the following method to the **NearestNeighborInterfaceTests** class:

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        x = [0.0, 1.0, 4.0, 7.5]
        y = [0.0] * len(x)
        z = [0.0] * len(x)
        
        n0, n1, n2, error = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(error[0], 0)
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()


This test calls the *find_nearest_neighbors* method with 4 positions
and checks if the nearest neighbors are determined correctly.
Let's run the test, and see if everything is working:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 test in 0.491s

    OK

We now have a simple interface that works, but we have to do our own 
indexing after the call and we could send data of any unit to the
method, also we have to do our own error checking after the method. Let's
define a object oriented interface to solve these problems

Defining the Object Oriented Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The object oriented interface sits on top of the legacy interface.
It decorates this interface with sets, unit handling, state engine
and more. We start creating the object oriented interface by inheriting
from InCodeComponentImplementation and writing the __init__ function. The build script
has added this class to the interface.py file for us. Open an editor
on interface.py and make sure this code is in the file (at the end
of the file):

.. code-block:: python

    
    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)

Configuring the handlers
~~~~~~~~~~~~~~~~~~~~~~~~

We configure the object oriented interface by implementing several 
methods. The object oriented interface is implement by several 
"handlers". Each handler provides support for a specific aspect of 
the interface. AMUSE defines a handler for the unit conversion, a 
handler for the interfacing with sets of particles, a handler to 
ensure the methods are called in the right order, etc. Each handler 
is very generic and needs to be configured before use. The handler 
are configured using the "Visitor" pattern. The following 
pseudo-code shows how the handlers are configured

.. code-block:: python

    class InCodeComponentImplementation(object):
        #...
        
        def configure_handlers(self):
            #...
            for handler in self.get_all_handlers():
                handler.configure(self)
        
        def define_converter(self, handler):
            """ configure the units converter handler """
            
            handler.set_nbody_converter(...)
            
        def define_particle_sets(self, handler):
            """ configure sets of particles """
            
            handler.define_incode_particle_set(...)
            handler.set_getter(...)
            
    class HandleConvertUnits(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_converter(self)
            
    class HandleParticles(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_particle_sets(self)
            
            
Configuration of the handlers is optional, we only have to define
those handler that we need in our interface. In our example we need
to configure the "HandleMethodsWithUnits" handler (to define units and
error handling) and the "HandleParticles" to define a particle set.

Defining methods with units
+++++++++++++++++++++++++++

We first want to add units and error handling to the 
**find_nearest_neighbors**. We do this by creating a 
**define_methods** function on the **NearestNeighbor** class. Open 
an editor on *interface.py* and add this method to the class:

.. code-block:: python

    def define_methods(self, handler):
        
        handler.add_method(
            "find_nearest_neighbors",
            (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                handler.INDEX,
                handler.INDEX,
                handler.INDEX,
                handler.ERROR_CODE
            )
        )
        
The **add_method** call expects the name of the function in the
legacy interface as it's first, next it expects a list of the units
of the input parameters and a list of the units of the output parameters.
The return value of a function is always the last item in the list
of output parameters. We specify a **generic_unit_system.length** unit
for the x, y and z parameters. The output parameters are indices and 
an errorcode. The errorcode will be handled by the AMUSE interface (0
means success and < 0 means raise an exception).

Let's write a test to see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test3(self):
        instance = NearestNeighbor()
        x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        y = [0.0] * len(x) | generic_unit_system.length
        z = [0.0] * len(x) | generic_unit_system.length
        
        n0, n1, n2 = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()

.. note::
    This test looks a lot like test2, but we now have to
    define a unit and we do not need to handle the errorcode.
    
Now build and test the code:

.. code-block:: bash
    
    > make clean; make all
    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.650s

    OK
    
.. note::
    
    Although we only edited python code we still need to run make. The
    code will check if the "worker_code" executable is up to date on 
    every run. It cannot detect if the update broke the code but it
    will still demand that the code is rebuilt.
    

Defining the particle set
+++++++++++++++++++++++++

Particle sets in AMUSE can be handled by python (we call these 
"inmemory") and by the legacy code (we call these "incode"). In our 
case the code does not handle the particles and we need to 
configure the particles handler to manage an inmemory particle set. 
Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    def define_particle_sets(self, object):
        object.define_inmemory_set('particles')

That's all we now have a "particles" attribute on the class and 
we can add, remove, delete particles from this set. But we are
still missing a connection between the particles and the nearest 
neighbors. AMUSE provides no handler for this, instead, we 
will write a method to run the find_nearest_neighbors function and
set the indices on the particles set.

Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    
    def run(self):
        indices0, indices1, indices2 = self.find_nearest_neighbors(
            self.particles.x,
            self.particles.y,
            self.particles.z
        )
        self.particles.neighbor0 = list(self.particles[indices0])
        self.particles.neighbor1 = list(self.particles[indices1])
        self.particles.neighbor2 = list(self.particles[indices2])
        
This function gets the "x", "y" and "z" attributes from the particles
set and sends these to the "find_nearest_neighbors" method. This methods
returns 3 lists of indices and we need to find the particles with
these indices. 

.. note::
    
    Particle sets have no given sequence, deletion and addition of 
    particles will change the order of the particles in the set. It
    is therefor never a good idea to use the index of the particle in the
    set as a reference to that particle. However, in the "run" method
    we "own" the particle set, it cannot change between the find_nearest_neighbor
    call and the moment we find the particles in the set by index (using
    self.particles[indices0]), and in this case it is save to use
    index as a valid reference.

Let's write a test and see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test4(self):
        instance = NearestNeighbor()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()

Now, make and run the tests:

.. code-block:: bash

    > make clean; make all
    > nosetests
    ....
    ----------------------------------------------------------------------
    Ran 4 tests in 0.797s

    OK
 
We are done, we have defined an object oriented interface on the
legacy interface. Only, if we look at our tests, the code seems
to be more rather than less complex. But, remember we now
have units and we are compatible with other parts of amuse. And
we can make more complex scripts easier.

Let's make a plummer model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash

    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png



Path 2
------

Defining the legacy interface
-----------------------------
We define our code interface so that a user can add, update and 
delete particles, start the nearest neighbors finding
algorithm and retrieve the ids of the nearest neighbors.

To define the interface, open interface.py with your favorite
editor and replace the contents of this file with:

.. literalinclude:: nearestneighbor/nn1.py

We can generate a stub from the interface code with::

    > $AMUSE_DIR/build.py --type=c --mode=stub interface.py NearestNeighborInterface -o interface.cc

The generated **interface.cc** replaces the original file generated in
the previous section.     

The code builds, but does not have any functionality yet::

    > make clean
    > make all
    
.. note::
    
    Compiling the interface code will result in a lot of warnings about
    unused dummy arguments. These warnings can be safely ignored for now.
    
The tests are broken (the echo_int function has been removed):

.. code-block:: bash

    > nosetests
    E
    ======================================================================
    ERROR: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "../src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "nearestneighbor/test_nearestneighbor.py", line 11, in test1
        result,error = instance.echo_int(12)
    AttributeError: 'NearestNeighborInterface' object has no attribute 'echo_int'

    ----------------------------------------------------------------------
    Ran 1 test in 0.315s

    FAILED (errors=1)

Let's create a working test by calling the new_particle method, open an editor
on the test_nearestneighbor.py file and replace the test1 method with:

.. code-block:: python

    def test1(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()

As this is python code we do not need to rebuild the code, instead
we can run the tests right after saving the code. Unfortunately, when
we run the test, it still fails.

.. code-block:: bash

    > nosetests
    F
    ======================================================================
    FAIL: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "/src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "/nearestneighbor/test_nearestneighbor.py", line 13, in test1
        self.assertEquals(result, 1)
      File "/src/amuse/test/amusetest.py", line 62, in failUnlessEqual
        self._raise_exceptions_if_any(failures, first, second, '{0} != {1}', msg)
      File "/src/amuse/test/amusetest.py", line 49, in _raise_exceptions_if_any
        raise self.failureException(msg or err_fmt_string.format(first, second, *args))
    AssertionError: 0 != 1
    -------------------- >> begin captured logging << --------------------
    legacy: INFO: start call 'NearestNeighborInterface.new_particle'
    legacy: INFO: end call 'NearestNeighborInterface.new_particle'
    --------------------- >> end captured logging << ---------------------

    ----------------------------------------------------------------------
    Ran 1 test in 0.319s

When you look closely at the output of the test you see that the result from the
method is 0 and not the expected 1. We need to edit the c code to make this 
test work. Open an editor on interface.cc and go to the *new_particle* function.

.. code-block:: c++

    int new_particle(int * index_of_the_particle, double x, double y, 
        double z)
    {
        *index_of_the_particle = 1;
        return 0;
    }

.. note ::
    
    In AMUSE all interface functions return an errorcode. Any other return
    values must be passed through the arguments of the functions.

We need to rebuild the code, and after building run the tests.

.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

The are tests work again! Only, we do not have any real working legacy code.

Filling the stubs
-----------------
The implementation of the algorithm does not match the interface
we defined and created. We need to write some glue code to
connect the code with the interface. To do so we fill in the
stubs generated earlier.

Open the **interface.cc** file in your favorite editor and
change its contents to:

.. literalinclude:: nearestneighbor/interface1.cc
    :language: c++

Test if the code builds and try it out. In the legacy interface
directory do:
    
.. code-block:: bash

    > make clean
    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.311s

    OK
    
Let's check some more functionality by adding another test

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 2)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, -1)
        error = instance.delete_particle(1)
        self.assertEquals(error, 0)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()
        

The tests should succeed:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 tests in 0.448s

    OK


We now have done everything in Step 0 *Legacy Interface*. We have
a legacy code and can access it in our python script. But, our
interface is not very friendly to work with. We have to think about
errorcodes and we have not information about units. To make our
interface easier to works with we start defining methods, properties
and parameters.

Defining methods
----------------
The object oriented interface is also defined in the **interface.py**.
So, we continue by opening an editor on this file. We will be 
writing methods for the **NearestNeighbor** class, in your editor
seek this code (at the end of the file)::

    class NearestNeighbor(CodeInterface):

        def __init__(self, **options):
            CodeInterface.__init__(self,  NearestNeighborInterface(), **options)
            
We will start by defining methods, we will do this by implementing
the **define_methods** function, like so::

    class NearestNeighbor(CodeInterface):

        def __init__(self, **options):
            CodeInterface.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            
            builder.add_method(
                "new_particle", 
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.INDEX, builder.ERROR_CODE)
            )
            
            builder.add_method(
                "delete_particle", 
                (builder.INDEX,),
                (builder.ERROR_CODE)
            )
            
            builder.add_method(
                "get_state", 
                (builder.INDEX,),
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, builder.ERROR_CODE),
                public_name = "get_position"
            )
            
            builder.add_method(
                "set_state", 
                (builder.INDEX, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.ERROR_CODE),
                public_name = "set_position"
            )
            
            builder.add_method(
                "run", 
                (),
                (builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_close_neighbors", 
                (builder.INDEX,),
                (builder.LINK('particles'), builder.LINK('particles'), builder.LINK('particles'), builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_nearest_neighbor", 
                (builder.INDEX,),
                (builder.LINK('particles'), generic_unit_system.length, builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_number_of_particles", 
                (),
                (builder.NO_UNIT, builder.ERROR_CODE),
            )

        
With this code, we define the methods and specify how to interpret the
arguments and return values. We get a special object (the **builder**
object) that provides us with the **add_method** function to be able
to this. The definition of the **add_method** function is as follows::
    
    add_method(
        name of the original function in the legacy interface,
        list of arguments (unit or type),
        list of return values,
        public_name = name for the user of the class (optional)
    )
    
For every argument or return value we can specify if it has a unit or
if it is special. The special arguments are:

=========================== ===================================
definition                  description
=========================== ===================================
builder.ERROR_CODE          the value is interpreted as an errorcode,
                            zero means no error for all other values and
                            Exception will be raise, only valid for one return value.
builder.NO_UNIT             the value has not unit (for example for the
                            number of items in a list)
builder.INDEX               the value is interpreted as an index for
                            object identifiers
builder.LINK('particles')   the value is interpreted as a link to
                            objects in the set with the given name
=========================== ===================================

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test:

.. code-block:: python

    def test3(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(2)
        instance.commit_parameters()
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length, 
            3.0 | generic_unit_system.length
        )
        self.assertEquals(result, 1)
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length
        )
        self.assertEquals(result, 2)
        x,y,z = instance.get_position(1)
        self.assertEquals(1.0 | generic_unit_system.length, x)
        self.assertEquals(2.0 | generic_unit_system.length, y)
        self.assertEquals(3.0 | generic_unit_system.length, z)
        instance.stop()

And run the tests:

.. code-block:: bash

    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.664s

    OK

As you can see our script is now a little simpler and we support units.
We do not have to think about the errorcodes in this script, AMUSE will
interpret the errorcodes and raise the right exceptions if needed. The 
units are also automatically converted to the right units for the 
code. But the script is still not very easy and we have to manage
all the ids we get from the code. To make our code even easier to
handle we will continue by defining a **set**.

.. note::

    We skip defining parameters and properties, we will come back to 
    this later in this tutorial.

Defining a set
--------------
We have made our interface a little easier but we still have to
do a some management work in our script. We would like to work
with objects and adding or removing these objects from the code.
AMUSE supports this by defining **sets**. Each set is capable of
storing specific attributes of the objects in the set. Our code
is capable of storing the x, y and z position of an object. An object
in AMUSE is called a *Particle* and the sets that contain these
particles are called *ParticleSets* or shorter *Particles*.

We define our particle set by implementing a **define_particle_sets** 
function on our **NearestNeighbor** class like so::

    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            ...
            
        def define_particle_sets(self, builder):
            builder.define_set('particles', 'index_of_the_particle')
            builder.set_new('particles', 'new_particle')
            builder.set_delete('particles', 'delete_particle')
            builder.add_setter('particles', 'set_position')
            builder.add_getter('particles', 'get_position')
            builder.add_getter('particles', 'get_close_neighbors', names=('neighbor0', 'neighbor1', 'neighbor2') )


That's all, we now have defined a set called **particles**. Again, we 
get a builder object to use in defining our set. All methods have
the name of the set as their first argument, this name can be any
name you want, but in AMUSE most codes provide a set called 
**particles**. For the **add_setter**, **add_getter**, **set_new** 
and **set_delete** functions, the second argument is the name of
the method we defined in the previous step. Finally you can set
the name of the attribute in the particles set with the **names**
argument. This is optional for legacy functions, if not given the names 
of the attributes will be derived from the names of the arguments in
the original calls. For example, the **get_position** call we specified
earlier has parameter name **x**, **y** and **z**, these names are 
also used in the particles set.

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test::
    
    def test4(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(100)
        instance.commit_parameters()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()
    
    
The support for a particle set means we can now also interact
with other parts of AMUSE. Let's make a plummer
model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash

    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png

