Installing on Ubuntu
====================

Installing on Ubuntu version > 10.10
------------------------------------

In this section we assume a default Ubuntu desktop installation.

All
---
The prerequisites can be installed with a couple of commands
on Ubuntu. The only choice to make is between openmpi and mpich. Most
of our testing is done with MPICH but openmpi should also work.

For openmpi do::

	> sudo apt-get install build-essential gfortran python-dev \
	  libopenmpi-dev openmpi-bin \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	> sudo easy_install mpi4py

For mpich do::
	
	> sudo apt-get install build-essential gfortran python-dev \
	  mpich libmpich-dev \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	> sudo easy_install mpi4py

.. note::
	
	Please make sure not to install mpich2 and openmpi together. 
	When both openmpi and mpich2 are installed strange errors
	will occur and AMUSE will not work. If you see both installed
	please remove both and install one. On older Ubuntu versions 
	the above package names for mpich should be mpich2 and libmpich2.

Installing on Ubuntu 9.04
~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a default Ubuntu desktop installation.

Python
------
Ubuntu comes with python2.6 pre-installed, you can check if
python is installed by doing:

.. code-block:: sh

	> python --version
	Python 2.6.2

If this fails with an error or a version before 2.6, please install 
python first(the package is called ``python2.6``). You also need 
the ``python2.6-dev`` development package.
To install it, do::

    > sudo apt-get install python2.6-dev
    

GCC
---
By default, Ubuntu does not install a fortran 90 or a C++ compiler. We
suggest using gfortran and g++. These compilers are installed with
the ``build-essential`` and the ``gfortran`` package. 
To install these, do::

    > sudo apt-get install build-essential gfortran

MPI2
----
Ubuntu does not provide installation packages for MPICH. You can 
build MPICH by hand (a good HOWTO can be found at 
https://wiki.ubuntu.com/MpichCluster). Or, you can download and install
pre-build packages from the MPICH site (http://www.mpich.org).

If you prefer OpenMpi over MPICH, you can install openmpi
from the Ubuntu packages. To install
the openmpi packages, do::

     > sudo apt-get install libopenmpi-dev openmpi-bin 

HDF5
----
Amuse can work with HDF5 versions 1.6.* and 1.8.3. Ubuntu 9.04 comes
with HDF5 version 1.6.6. To install it, do::

    > sudo apt-get install libhdf5-serial-dev hdf5-tools 

FFTW
----
On Ubuntu, FFTW can be installed with::

    > sudo apt-get install libfftw3 libfftw3-dev libfftw3-doc

GSL
-------
On Ubuntu, GSL can be installed with::

    > sudo apt-get install libgsl0 libgsl0-dev

CMake
-------
CMake is used to build EVTwin. On Ubuntu, CMake can be installed with::

    > sudo apt-get install cmake

GMP
-------
GMP is required for Adaptb. On Ubuntu, GMP can be installed with::

    > sudo apt-get install libgmp3 libgmp3-dev

MPFR
-------
MPFR is required for Adaptb. On Ubuntu, MPFR can be installed with::

    > sudo apt-get install libmpfr4 libmpfr-dev

Python packages in Ubuntu
-------------------------
Ubuntu comes with python packages for nose and numpy. You also need 
the setuptools package to be able to install the ``mpi4py`` and ``h5py`` 
software. To install these , do::

    > sudo apt-get install python-nose python-numpy python-setuptools python-docutils

Python packages with easy_install
---------------------------------
The ``mpi4py`` and ``h5py`` can be installed with the ``easy_install``
command::

    > sudo easy_install mpi4py
    > sudo easy_install h5py
    
Installing on Ubuntu 9.10
~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a default Ubuntu desktop installation. This 
installation is for the most part the same as for Ubuntu 9.04, see 
previous section. 

The development packages of python are needed, to install these do::

    > sudo apt-get install python-dev 

FFTW
-------
For 9.10 the FFTW package name is fftw3 and not libfftw3, FFTW can be installed with::

    > sudo apt-get install fftw3 fftw3-dev fftw3-doc

