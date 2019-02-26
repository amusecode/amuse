This package contains the AMUSE software. With AMUSE you can write scripts to
simulate astrophysical problems in different domains.

This release (12.0a4) does not contain the following codes, for space-saving reasons:
- MESA
- SimpleX
- EVTwin

The documentation and the software can be found at:

* https://github.com/amusecode/amuse
* http://www.amusecode.org (currently outdated)

Getting Started
===============

To install AMUSE you need to first install some prerequisites. An MPI
installation is needed, both OpenMPI and MPICH are known to work.
The installation method below assumes a Python2 environment, but Python3 has
been tested with some success.

## Installing prerequisites

### macOS:

Use a package manager (e.g. macports) to install the prerequisites:
```bash
sudo port install gcc7 openmpi-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
sudo port install python27 py27-virtualenv 
```
Next, set the just installed compilers to be the default.
```bash
sudo port select --set mpi openmpi-gcc7-fortran
sudo port select --set gcc mp-gcc7
sudo port select --set python2 python27
```

### Linux (Ubuntu):

For openmpi:

```bash
sudo apt-get install build-essential gfortran python-dev \
libopenmpi-dev openmpi-bin \
libgsl0-dev cmake libfftw3-3 libfftw3-dev \
libgmp3-dev libmpfr4 libmpfr-dev \
libhdf5-serial-dev hdf5-tools \
python-nose python-numpy python-setuptools python-docutils \
python-h5py python-setuptools git
```

For mpich:

```bash
sudo apt-get install build-essential gfortran python-dev \
mpich libmpich-dev \
libgsl0-dev cmake libfftw3-3 libfftw3-dev \
libgmp3-dev libmpfr4 libmpfr-dev \
libhdf5-serial-dev hdf5-tools \
python-nose python-numpy python-setuptools python-docutils \
python-h5py python-setuptools git
```

## Installing AMUSE
The preferred way of installing AMUSE is in a clean virtual environment.
In this environment, first install wheel and mpi4py:
```bash
pip install wheel mpi4py
```

Then, install AMUSE and its prerequisite packages (this will take a long time):
```bash
pip install amuse
```

Done!
