===========================================
Simplified interface function specification
===========================================

In addition to the specification of remote interface functions by using 
the ``@legacy_function`` decorator, a simplified method is available, 
which can handle most cases of interface functions. Here we will 
describe the simplified interface function specification.

Example
~~~~~~~
Let's start with a simple example specification:

.. code-block:: python
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.addParameter('sum', 'd', function.OUT)
            function.result_type = 'i'
            function.can_handle_array = True
            return function

This can be converted to:

.. code-block:: python
             
        @remote_function(can_handle_array=True)
        def sum(x='d',y='d'):
            returns (sum='d')

As can be seen the parameters are specified in keyword/value style. 
``can_handle_array`` and ``must_handle_array`` are supplied, if 
necessary, as keywords to the decorator. A default integer return value 
(for the error code) is implied (but can be overridden, see below). The 
following table lists the options for the parameter specifications:

=========  ================  ===============================
data type  no default value  default value (for input) 
=========  ================  ===============================
boolean    "b", "bool"       True, False
integer    "i", "int32"      <int>, numpy.int32(<value>)
long       "l", "int64"      numpy.int64(<value>)
float      "f", "float32"    numpy.float32(<value>)
double     "d", "float64"    <float>, numpy.float64(<value>)
string     "s", "string"     "any other string"
=========  ================  ===============================

A unit specification can be added. Remember that parameters without 
default cannot follow parameters with default. So the following will 
generate an error:

.. code-block:: python
             
        @remote_function(can_handle_array=True)
        def sum(x=0.,y='d'):
            returns (sum='d')

Below are some more examples of valid specifications:

.. code-block:: python
             
        @remote_function
        def initialize_code():
            pass

        @remote_function
        def get_time():
            returns (time=0. | units.s)

        @remote_function(must_handle_array=True)
        def inout_sum(x='d',y='d', sum=0.):
            returns (sum=0.)

        @remote_function
        def function_with_float_return():
            return (__result=0.)

One limitation of this type of specification is that they won't
work if generated dynamically (so don't try, for example, to generate
a bunch of functions based on a list of parameter names).
