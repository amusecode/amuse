=======================================
Create an low-level Interface to a Code
=======================================

In this tutorial we will create an interface to a code. 

.. warning::
    
    This tutorial does not use the build script provided with 
    AMUSE. All files and directories need to be created "by hand".
    Use this tutorial if you want to get a deeper understanding of
    how the build process works and which steps are involved in
    creating a low level interface. To learn how to create
    an interface to a code we recommend :doc:`c_code` and 
    :doc:`fortran_code`.


Work directory
~~~~~~~~~~~~~~
We start by making a directory for our code. This directory should
be a subdirectory of the "src/amuse/community" directory. It also will be
a python package and we need to create the file "__init__.py" in 
this directory. So, let's open a shell and go to the AMUSE 
root directory. To create the code directory we then do:

.. code-block:: bash
    
    >> cd src/amuse/community
    >> mkdir mycode
    >> cd mycode
    >> touch __init__.py
    >> pwd
    ../src/amuse/community/mycode
    >> ls
    __init__.py
    

The code
~~~~~~~~
To create an interface we first need the code. For
this example we will use a very simple code do some calculations 
on two numbers.

The contents of the code.c file:

.. code-block:: cpp

  #include "code.h"
   
  double sum(double x, double y) {
    return x + y;
  }
  
  int divide(double x, double y, double * result) {
    if(y == 0.0) {
        return -1;
    } else {
        *result = x / y;
        return 0;
    }
  }

We need to access these function from another C file, so we need to
define a header file. 

The contents of the code.h file:

.. code-block:: cpp

    double sum(double x, double y);
    
    int divide(double x, double y, double * result);


The interface code
~~~~~~~~~~~~~~~~~~
Now we can define the interface class for our code in python. An 
interface needs to inherit from the class "CodeInterface".

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
In this example we first import names from the ``amuse.community`` 
module on line 1. The ``amuse.community`` module defines the typical 
classes and function needed to write a legacy interface. On line 3 
we define our class and inherit from ``CodeInterface``. The class 
will be used to generate a C++ file. In this file we will need the 
definitions of our legacy functions. So, on line 4 we specify the 
necessary include files in a array of strings. Each string will be
converted to an include statement.

Building the code
~~~~~~~~~~~~~~~~~
Building the code takes a couple of steps, first generating the C file
and then compiling the code. We will construct a makefile to automate
the build process.

.. code-block:: make
    :linenos:
    
    ifndef AMUSE_DIR
        AMUSE_DIR=../../../..
    endif

    CODE_GENERATOR = $(AMUSE_DIR)/build.py

    CXXFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
    LDFLAGS = -lm $(MUSE_LD_FLAGS)

    OBJS = code.o

    all: worker_code

    cleanall: clean
        $(RM) worker_code *~
        
    clean:
        rm -f *.so *.o *.pyc worker_code.cc

    worker_code.cc: interface.py
        $(CODE_GENERATOR) --type=c interface.py MyCode -o $@

    worker_code: worker_code.cc $(OBJS)
        mpicxx $@.cc $(OBJS) -o $@

    .cc.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<
        
    .c.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<

.. compound:

    You need to convert the spaces into tabs, 
    if you copy the above text to a new file.
    

Let's start make and build the ``worker_code`` application

.. code-block:: bash
    
    >> make clean
    >> make
    ...
    mpicxx worker_code.cc code.o -o worker_code
    >> ls
    ... worker_code ...
    
Running the code
~~~~~~~~~~~~~~~~
Before we run the code we need to start the MPI daemon process ''mpd''. 
This daemon process manages the start of child processes.

.. code-block:: bash
    
    >> mpd &

We can use ``amuse.sh`` and try the interface.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance
    <amuse.community.mycode.interface.MyCode object at 0x7f57abfb2550>
    >>> del instance
    
We have not defined any methods and our interface class is not
very useful. We can only create an instance of the code. When we 
create this instance the "worker_code" application will start 
to handle all the function calls. We can see the application 
running when we do "ps x | grep worker_code"

Implementing a method
~~~~~~~~~~~~~~~~~~~~~~
Now we will define the ``sum`` method. We will add the definition to
the MyCode class.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function
            
The new code starts from line 9. On line 9 we specify a decorator. This
decorator will convert the following function into a specification that
can be used to call the function and generate the C++ code. On line 10
we give the function the same name as the function in our code. This
function may not have any arguments. On line 11 we create an instance of
the "LegacyFunctionSpecification" class, this class has methods to specify the intercase.
On line 12 and 13 we specify the parameters for out functions. Parameters have
a name, type and direction. The type is specified with a single character
*type code*. The following type codes are defined:
            
            
=========  ===========  ================
Type code  C type       Fortran type  
=========  ===========  ================
'i'        int          integer
'd'        double       double precision
'f'        float        single precision
=========  ===========  ================

The direction of the parameter can be ``IN``, ``OUT`` or ``INOUT``. On line 
14 we define the return type, this can be a *type code* or ``None``. The default
value is ``None``, specifying no return value (void function).

Let's rebuild the code.

.. code-block:: bash
    
    >> make clean
    >> make
    ...
    mpicxx worker_code.cc code.o -o worker_code

We can now start ```amuse.sh``` again and do a simple sum

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance.sum(40.5, 10.3)
    50.799999999999997
    >>> 40.5 + 10.3
    50.799999999999997
    >>> del instance

And we see that our interface correctly sums two numbers.

A method with an OUT parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We can complete out interface by defining the ``divide`` function.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function

        @legacy_function
        def divide():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.addParameter('result', 'd', function.OUT)
            function.result_type = 'i'
            return function
            
On line 22 we define the parameter "result" as an OUT parameter. In python
we do not have to provide this parameter as an argument to our function. It
After rebuilding we can try this new function.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> (result, error) =  instance.divide(10.2, 30.2)
    >>> result
    0.33774834437086093
    >>> error
    0
    >>> del instance

We see that the function returns two values, the OUT parameter and also
the return value of the function.

Working with arrays
~~~~~~~~~~~~~~~~~~~
Some functions will be called to perform on the elements of an array. 
For example:

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> x_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    >>> y_values = [10.3, 20.3, 30.4 , 40.4, 50.6]
    >>> results = []
    >>> for x , y in map(None, x_values, y_values):
    ...     results.append(instance.sum(x,y))
    ...
    >>> print results
    [11.300000000000001, 22.300000000000001, 33.399999999999999, 
    44.399999999999999, 55.600000000000001]
    
    
The MPI message passing overhead is incurred for every call on 
the method. We can change this by specifing the function to be able
to handle arrays.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            function.can_handle_array = True
            return function

On line 15 we specify that the function can be called with an array of
values. The function will be called for every element of the array. The 
array will be send in one MPI message, reducing the overhead.

Let's rebuild the code and run an example.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> x_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    >>> y_values = [10.3, 20.3, 30.4 , 40.4, 50.6]
    >>> results = instance.sum(x_values, y_values)
    >>> print results
    [ 11.3  22.3  33.4  44.4  55.6]
    

Other interfaces
~~~~~~~~~~~~~~~~
The community codes directory contains a number of codes. Please look at
these codes to see how the interfaces are defined.


             


