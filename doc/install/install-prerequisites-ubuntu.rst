Installing on Ubuntu
====================

Installing on Ubuntu (up to date for version 18.04)
---------------------------------------------------

In this section we assume a default Ubuntu desktop installation.

All
---
The prerequisites can be installed with a couple of commands
on Ubuntu. The only choice to make is between openmpi and mpich2. 

For openmpi do::

	> sudo apt-get install build-essential gfortran python-dev \
	  libopenmpi-dev openmpi-bin \
	  libgsl-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr6 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	>  [sudo] pip install mpi4py
	or alternatively setuptools easy_install (deprecated):
	>  [sudo] easy_install mpi4py


For mpich do::
	
	> sudo apt-get install build-essential gfortran python-dev \
	  mpich libmpich-dev \
	  libgsl-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr6 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	>  [sudo] pip install mpi4py
	or alternatively setuptools easy_install (deprecated):
	>  [sudo] easy_install mpi4py

.. note::
	
	Please make sure not to install mpich2 and openmpi together. 
	When both openmpi and mpich2 are installed strange errors
	will occur and AMUSE will not work. If you see both installed
	please remove both and install one.
