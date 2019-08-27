This directory contains the AMUSE software. With AMUSE you can write
scripts to simulate astrophysical problems in different domains.

The documentation and more info can be found at:

* http://www.amusecode.org

Getting Started
===============

In short, most probably

```bash
pip install amuse
```
should get you going if you have a linux or Mac were you compile 
codes on (HDF5 and an MPI libraries must be installed). 

Below are some hints for a quick install, if these fail please 
look for options at the detailed descriptions of the installation 
procedure in the documents in the 'doc/install' directory.

Compilers
=========

To build AMUSE from source you need to have a working build
environment. The AMUSE build system needs C/C++ and fortan 90
compilers, we recommend a recent version of GCC. 

In Ubuntu you can setup the environment with (as root):

```bash
apt-get install build-essential curl g++ gfortran gettext zlib1g-dev
```

Other distributions have similar package or package groups available.

In OS X you can use the homebrew or macports package manager (both
require the Apple Developer Tools and Xcode to be installed).

Python
======

AMUSE needs Python 2, version >2.7, or Python3 version >=3.5 installed
preferably with pip and virtualenv. It may be necessary to update pip
to a recent version.

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

This will build and install AMUSE with an extensive set of codes.
If necessary this will also install some required Python packages:

* Numpy (version >= 1.3.0)
* h5py (version >= 1.2.0)
* mpi4py (version >= 1.0)
* nose (version >= 0.11)
* docutils (version >= 0.6)

If you are not using pip these must be installed by hand.

It is possible to install the minimal framework by:

```bash
pip install [--user] amuse-framework
```

This does not include any codes. These can be added
```bash
pip install [--user] amuse-<code name>
```

AMUSE Development 
=================

If you are using Python 2, an AMUSE development install can also 
be handled through pip by executing (in the root of a clone of the 
repository)

```bash
pip install -e .
```

after this the codes need to be build:

```bash
python setup.py develop_build
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
