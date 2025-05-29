====================
Directory Structure
====================

The amuse source-code is separated into 3 directories:

* ``src`` - source code, implementation of the environment.
* ``test`` - applications, examples and unittests.
* ``support`` - build system, test system.

Under the ``src`` directories all code needed to run AMUSE can be found.
One can view this code as an *library* that can be used to create
*applications* to do numerical astrophysical experiments. This code will
contain the building blocks needed to interface with codes,
import and export data, do unit conversions, and all other AMUSE
functionality.

Under the ``test`` directories all application and test code can be
found. This directory tree will contain scripts to do a complete
astrophysical experiment. Also all unit-tests can be found here. These
unit tests each cover only a small part (unit) of the functionality of
AMUSE. For example a  test to check the import of a file to AMUSE data
format.

Under the ``support`` directories all support code for the building
system can be found.

The ``src`` directories
~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``src`` directory are further split into:

* ``amuse`` - contains the AMUSE framework, the *glue* that makes it
  possible to access the community codes from a Python script and
  combine them together

* ``amuse_<code>`` - contains the community codes and the AMUSE wrappers
  that allow them to be used from AMUSE.


The ``test`` directories
~~~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``test`` directory are further split into:

* ``unit_tests`` - All unit testing code. These tests are coded using
  the standard unit testing framework that is included in the Python
  distribution (``unittest``). See python module documentation for further
  information: http://docs.python.org/library/unittest.html.
* ``application`` - contains the source code of published applications.
* ``examples`` - contains documented example codes.

The ``support`` directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``support`` directory are further split into:

* ``test`` - Scripts to support the testing of AMUSE code.
* ``build`` - Scripts used by the building system of AMUSE.



