Installing on RedHat (CentOS)
=============================

Installing on CentOS 6
~~~~~~~~~~~~~~~~~~~~~~~

In this section we asume a minimal CentOS 6 installation.

All
---
The prerequisites can be installed with a couple of commands
on CentOS 6. 

To install the prerequisites do (for base-devel select *all* members)::

    > sudo yum install make gcc gcc-c++ gcc-gfortran \
	cmake zlib-devel\
	openmpi openmpi-devel \
	fftw fftw-devel \
	gsl gsl-devel gmp
	
After installing openmpi, you need to activate it using the 'module'
command::
    
    > module load openmpi-$(uname -i)

.. note::

    We recommend to put the openmpi module activation script
    in your .bashrc or .cshrc file.

Install python and dependencies::

    > sudo yum install python-devel \
	docutils python-nose \
	numpy numpy-f2py\
	python-docutils

To install hdf5 and docutils first install an additional rpm forge.
For documentation see http://wiki.centos.org/AdditionalResources/Repositories/RPMForge

After installing an rpm forge do::

    > sudo yum install hdf5 hdf5-devel
    
To install h5py do::

    > sudo easy_install h5py
    
Last, you need to install mpi4py with::

    > su -
    > module load openmpi-$(uname -i)
    > easy_install mpi4py

.. note::
    
    The default CentOS sudo policy resets the environments variables and
    thereby removes the openmpi settings. So for the last step
    you cannot use ```sudo easy_install mpi4py``` but must install
    under root directly.

