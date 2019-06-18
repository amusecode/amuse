This directory contains the AMUSE software. With AMUSE you can write
scripts to simulate astrophysical problems in different domains.

The documentation and more info can be found at:

* http://www.amusecode.org

Getting Started
===============

To build amuse you need a working build environment with python and 
install some prerequisites. This document contains the quick install
instructions, if these fail please look for options at the detailed 
descriptions of the installation procedure in the documents in the 
'doc/install' directory.

Compilers
=========

To build AMUSE from source you need to have a working build
environment. The AMUSE build system needs a C++ and a fortan 90
compiler. Please check first if you have a working build environment
on your system.

In Ubuntu you can setup the environment with (as root):

```bash
apt-get install build-essential curl g++ gfortran gettext zlib1g-dev
```

In Fedora you can setup the environment with (as root)::

```bash
yum groupinstall "Development Tools" "Development Libraries"
```

In OS X you can install homebrew or macports package managers (both
need the Apple Developer Tools). If you do not want to use any of
these package managers you will need to install a fortran compiler
as the Apple Developer Tools do not include a fortran compiler, you
can find one at:

* http://hpc.sourceforge.net/

Finally, AMUSE needs Python 2 version >2.7 or Python3 version >=3.5 installed
preferably with pip and virtualenv.

Installing Prerequisites
========================

The following libraries need to be installed:

* HDF (version 1.6.5 - 1.8.x)
* MPI (OpenMPI or MPICH)

The following are needed for some codes:
* FFTW (version >= 3.0)
* GSL
* CMake (version >= 2.4)
* GMP (version >= 4.2.1)
* MPFR (version >= 2.3.1)

Installing+building AMUSE
=========================

AMUSE can be installed through pip:

```bash
pip install [--user] amuse
```

if necessary this will also install some required Python packages:

* Numpy (version >= 1.3.0)
* h5py (version >= 1.2.0)
* mpi4py (version >= 1.0)
* nose (version >= 0.11)
* docutils (version >= 0.6)

AMUSE Development 
=================

A install for AMUSE development can also be handled through pip, by executing
in the root of a clone of the repository

```bash
pip install -e .
```
after this the codes need to be build:

```bash
python seteup.py develop_build
```

Running the tests
=================
AMUSE comes with a large set of tests, most can be run automatically.
To run these tests start the nosetests command from the main
amuse directory (directory this README file lives in).

To run these tests do:

1. install the tests

```bash
pip install [--user] amuse-tests
```

2. Run the automatic tests

```bash
nosetests -v amuse.tests.suite
```
