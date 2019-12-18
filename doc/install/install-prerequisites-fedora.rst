

Installing on Fedora 18
~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a basic install of Fedora 18 installation.


All in One
----------
The prerequisites can be installed with a couple of commands
on Fedora.

For mpich2 do::
	
	> sudo yum install make gcc gcc-c++ gcc-gfortran\
		cmake zlib-devel\
		mpich2 mpich2-devel\
		hdf5 hdf5-devel\
		fftw fftw-devel\
		gsl gsl-devel\
		gmp gmp-devel\
		mpfr mpfr-devel\
		python-nose numpy numpy-f2py\
		h5py\
		python-setuptools python-setuptools-devel\
		mpi4py-mpich2\
		python-matplotlib
    


.. note::
    This line will also install `matplotlib`, this package is
    used for all plotting in AMUSE. If you do not need any plotting
    you can leave it out.

After installing mpich2, you need to activate it using the 'module'
command::

	> module load mpi/mpich2-$(uname -i)
	
.. note::

    We recommend to put the module activation script
    in your .bashrc or .cshrc file.

For openmpi do::

	> sudo yum install make gcc gcc-c++ gcc-gfortran\
		cmake zlib-devel\
		openmpi openmpi-devel\
		hdf5 hdf5-devel\
		fftw fftw-devel\
		gsl gsl-devel\
		gmp gmp-devel\
		mpfr mpfr-devel\
		python-nose numpy numpy-f2py\
		h5py\
		python-setuptools python-setuptools-devel\
		mpi4py-openmpi\
		python-matplotlib
	
After installing openmpi, you need to activate it using the 'module'
command::

	> module load mpi/openmpi-$(uname -i)


.. note::

    On Fedora you can install both mpich2 and openmpi, the module
    command will keep manage these separate installation, so
    no conflict will exists. If you change between implementation, you
    will need to recompile the amuse community codes with::
    
	> make clean; make
	


Installing on Fedora 11
~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a live-cd install of Fedora 11 installation.

Python
------
Fedora comes with python2.6 pre-installed, you can check if
python is installed by doing:

.. code-block:: sh

	> python --version
	Python 2.6.2

If this fails with an error or a version before 2.6, please install 
python first(the package is called ``python``). You also need 
the ``python-devel`` development package.
To install it, do::

    > sudo yum install python-devel
    

GCC
---
By default, Fedora does not install a fortran 90 or a C++ compiler. We
suggest using gfortran and g++. These compilers are installed with
the ``gcc``, ``gcc-c++`` and the ``gcc-gfortran`` packages. 
To install these, do::

    > sudo yum install gcc gcc-c++ gcc-gfortran

MPI2
----
Fedora comes with packages for MPICH2 and Openmpi.

To install MPICH2, do::
    
    > sudo yum install mpich2 mpich2-devel

If you prefer OpenMpi over MPICH2, you can install openmpi
from the Fedora yum database. 
To install the openmpi packages, do::

     > sudo yum install openmpi openmpi-devel

HDF5
----
Amuse can work with HDF5 versions 1.6.* and 1.8.3. Fedora 11 has a package
with HDF5 version 1.8.3. To install it, do::

    > sudo yum install hdf5 hdf5-devel

FFTW
-------
On Fedora, FFTW can be installed with::

    > sudo yum install fftw fftw-devel

GSL
-------
On Fedora, GSL can be installed with::

    > sudo yum install gsl gsl-devel

CMake
-------
CMake is used to build EVTwin. On Fedora, CMake can be installed with::

    > sudo yum install cmake

GMP
-------
GMP is required for Adaptb. On Fedora, GMP can be installed with::

    > sudo yum install gmp

MPFR
-------
MPFR is required for Adaptb. On Fedora, MPFR is currently included in the gmp 
package. So, if you have not already done so, MPFR can be installed with::

    > sudo yum install gmp

Python packages in Fedora
-------------------------
Fedora comes with python packages for nose and numpy. You also need 
the setuptools package to be able to install the ``mpi4py`` and ``h5py`` 
software. To install these , do::

    > sudo yum install python-nose numpy numpy-f2py \
        python-setuptools python-setuptools-devel

Python packages with easy_install
---------------------------------
The ``mpi4py``, ``h5py`` and ``docutils`` can be 
installed with the ``easy_install`` command::

    > sudo easy_install mpi4py
    > sudo easy_install h5py
    > sudo easy_install docutils 
