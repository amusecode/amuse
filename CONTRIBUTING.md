Contributing to AMUSE
=====================

So you're interested in contributing code to AMUSE? Excellent! 

Reporting Issues
----------------

When opening an issue to report a problem, please try and provide a minimal
code example that reproduces the issue, and also include details of the
operating system, compiler, and the Python, Numpy, and AMUSE versions you are using.

Contributing code
-----------------

Most contributions to AMUSE are done via pull requests from GitHub users'
forks of the [amuse repository](https://github.com/amusecode/amuse).

Once you open a pull request (which should be opened against the ``master``
branch, not against any of the other branches), please make sure that you
include the following:

- **Code**: the code you are adding

- **Tests**: these are usually tests to ensure that code that previously
  failed now works (regression tests) or tests that cover as much as possible
  of the new functionality to make sure it doesn't break in future, and also
  returns consistent results on all platforms (since we run these tests on many
  platforms/configurations). 


Checklist for Contributed Code
------------------------------

A pull request for a new feature will be reviewed to see if it meets the
following requirements.  For any pull request, an AMUSE maintainer can
help to make sure that the pull request meets the requirements for inclusion
in the package.

**Scientific Quality**
(when applicable)
  * Is the submission relevant to AMUSE?
  * Are references included to the origin paper for the simulation code?
  * Does the code perform as expected?
  * Has the code been tested against previously existing codes in the same domain?

**Code Quality**
  * Is the code compatible with Python >=2.7?
  * Are there dependencies other than AMUSE, MPI, the Python Standard
    Library, and NumPy?
  * For compatibility reasons we prefer code that also works on older 
    versions of Numpy, matplotlib etc.
  * Are additional dependencies handled appropiately? If possible, factor out 
    additional dependencies or make them optional.
  * Does the code follow the AMUSE Style Guide (http://www.amusecode.org/doc/reference/style_guide.html)?

**Testing**
  * Are the inputs to the functions sufficiently tested?
  * Are there tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Does python setup.py test run without failures?

**Documentation**
  * Is there any information needed to be added to the docs to describe the code?

**License**
  * Does the code require a specific license other than the AMUSE license?
  * Are there any conflicts with this code and AMUSE?

**AMUSE requirements**
  * Can you checkout the pull request and repeat the examples and tests?
