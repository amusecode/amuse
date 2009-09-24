Tutorial - Create an Interface to a Legacy Code
================================================

In this tutorial we will follow all steps to create an interface to 
a legacy code.

We start by making a directory for out legacy code. First, go 
to the AMUSE root directory and add a directory to the
legacy codes directory

.. code-block:: bash
    
    >> cd src/amuse/legacy
    >> mkdir mycode
    >> cd mycode
    >> touch __init__.py
    >> pwd
    ../src/amuse/legacy/mycode
    >> ls
    __init__.py

To create an interface to a legacy code we need a legacy code. For
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



Now we can define the interface class for our code in python.

.. code-block:: python

    from amuse.legacy import *
    
    class MyCode(LegacyInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             LegacyInterface.__init__(self)

We construct a makefile to build all the code. This makefile can run
as is in our directory.

.. code-block:: make

    CXXFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
    LDFLAGS = -lm $(MUSE_LD_FLAGS)

    OBJS = code.o

    all: muse_worker

    cleanall: clean
        $(RM) muse_worker *~
        
    clean:
        rm -f *.so *.o *.pyc muse_worker.cc

    muse_worker.cc: interface.py
        ../../../../bin/create_c_worker.py interface.py MyCode > $@

    muse_worker:	muse_worker.cc $(OBJS)
        mpicxx $@.cc $(OBJS) -o $@

    .cc.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<
        
    .c.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<
        
We can now start ```amuse.sh``` and try-out the interface

.. code-block:: pycon

    >>> from amuse.legacy.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance
    <amuse.legacy.mycode.interface.MyCode object at 0x7f57abfb2550>
    >>> del instance
    
We have not defined any methods and our interface class is not
very useful. We can only create an instance. When we create an instance
the "muse_worker" application will start to handle all the function
calls. We can see the application running when 
we do "ps -x | grep muse_worker"

Now we will define the sum method. We must create a method in the MyCode
class.

.. code-block:: python

    from amuse.legacy import *
    
    class MyCode(LegacyInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             LegacyInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = RemoteFunction()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function
            
Rebuild the code

.. code-block:: bash
    
    >> make clean
    >> make
    ...
    mpicxx muse_worker.cc code.o -o muse_worker


We can now start ```amuse.sh``` again and try-out our new interface

.. code-block:: pycon

    >>> from amuse.legacy.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance.sum(40.5, 10.3)
    50.8
    >>> del instance


.. code-block:: python

    from amuse.legacy import *
    
    class MyCode(LegacyInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             LegacyInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = RemoteFunction()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function

        @legacy_function
        def divide():
            function = RemoteFunction()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.addParameter('result', 'd', function.OUT)
            function.result_type = None
            return function
            
.. code-block:: pycon

    >>> from amuse.legacy.mycode import interface
    >>> instance = interface.MyCode()
    >>> (result, error) =  instance.divide(10.2, 30.2)
    >>> result
    0.33774834437086093
    >>> error
    0
    >>> del instance

             


             


