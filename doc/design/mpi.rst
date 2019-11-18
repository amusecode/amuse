=============
MPI interface
=============

The interface between the AMUSE python core and the legacy-codes is 
based on the MPI framework. Choosing MPI and not SWIG (or any other direct
coupling method) has several advantages:

* MPI is a well-known framework in the astrophysics community. 
  Other coupling methods are less well known (like SWIG)
* Legacy code does not run in the python space (memory usage, names)
* Multiple instances of the same legacy code can easily be supported (not so
  in SWIG / f2py couplings)
* Multi-process support taken into account at the start of 
  the project.
* Coupling is much looser.

There are also be some disadvantages:

* Need to define a protocol over MPI
* More "hand-work" needed to couple code. Other frameworks, like SWIG and f2py,
  generate an interface based on the application code.
* More overhead for every call, slower calls

These disadvantages are mitigated by creating a library that handles
most of the coupling details. This library has a Python, C++ and
Fortran version. It implements the protocol and generates
hooks to connect with the legacy codes.

The overhead per call may be an important factor in the speed of the
framework. This will be tested during development of the first codes. It
should be possible to limit the overhead by sending a lot of data per call. For
example, setting the properties of a lot of stars in one call. Calling a lot of methods
with limited data will be compared to sending one method with a lot of data.


