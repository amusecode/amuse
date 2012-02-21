Installing on Suse Linux
========================

Installing on OpenSuse 11
~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we asume a normal desktop install of OpenSuse 11. Not
all packages are available in the default OpenSuse package repository.
We recommend to add the **Packman Repository** to the list of 
configured sofware reporistories (To do so, open Yast and go to 
*Software Repositories*).

Python
------
OpenSuse comes with python2.6 pre-installed, you can check if
python is installed by doing:

.. code-block:: sh

	> python --version
	Python 2.6

If this failes with an error or a version before 2.6, please install 
python first(the package is called ``python``). You also need 
the ``python-devel`` development package.
To install it, do::

    > sudo zypper install python-devel
    

GCC
---
By default, OpenSuse does not install a fortran 90 or a C++ compiler. We
suggest using gfortran and g++. These compilers are installed with
the ``gcc``, ``gcc-c++`` and the ``gcc-fortran`` packages. 
To install these, do::

    > sudo zypper install gcc gcc-c++ gcc-fortran

MPI2
----
The Packman Repository provides an OpenMPI package.
To install the openmpi packages, do::

    > sudo zypper install openmpi openmpi-devel

Unfortunately the openmpi installation does not work out
of the box, you need to set the  **LD_LIBRARY_PATH** variable
and edit a configuration file first.

Setting the LD_LIBRARY_PATH
****************************

The LD_LIBRARY_PATH must be set so that mpi4py can find the
openmpi libraries. To set the variable we must first find out
where the openmpi libs can be found, to do so execute::

    > mpicxx -showme:link
    -pthread -L/usr/lib/mpi/gcc/openmpi/lib -lmpi_cxx -lmpi 
    -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
    

We need to set LD_LIBRARY_PATH variable to the path after the **-L**
in the output (so in this example case '/usr/lib/mpi/gcc/openmpi/lib',
this may be a different path if you system is 64-bits or if the
opensuse version is different).

In bash do::
    
    > export LD_LIBRARY_PATH=/usr/lib/mpi/gcc/openmpi/lib
    
We recommend you add this line to your '.bashrc' file so that
the variable is set correctly for all sessions. If you have a
C shell you need to do a *setenv* and edit the .cshrc file.

Editing the configuration file
*******************************

It seems that the default openmpi installation has some problems
with loading an LDAP library. To check if your installation has 
this problem do::

    > python -c "from mpi4py import MPI; print MPI.get_vendor()"
    ...
    WARNING: ....
    ...
    DAT: library load failure: libdaplscm.so.2: cannot open shared object file: No such file or directory
    ...

If you get a long list of warings about DAT providers not found, you
need to edit the configuration file and turn off ldap. To do so, 
open an editor (as root) on the file 
**/etc/openmpi-mca-params.conf** 
and add this line to the bottom of the file::

    btl = ^udapl

After saving the file, you can rerun the python statement::

    > python -c "from mpi4py import MPI; print MPI.get_vendor()"
    ('Open MPI', (1, 2, 8))
    
    
HDF5
----
Amuse can work with HDF5 versions 1.6.* and 1.8.*. The Packman Repository
has a package with HDF5 version 1.8.1. To install it, do::

    > sudo zypper install hdf5 hdf5-devel

FFTW
-------
Some codes in AMUSE need FFTW 3, FFTW can be installed with::

    > sudo zypper install fftw3 fftw3-devel

GSL
-------
On OpenSuse (10.2 and newer), GSL can be installed with::

    > sudo zypper install gsl gsl-devel

CMake
-------
CMake is used to build EVTwin. On OpenSuse, CMake can be installed with::

    > sudo zypper install cmake

GMP
-------
GMP is required for Adaptb. On OpenSuse, GMP can be installed with::

    > sudo zypper install gmp-devel

MPFR
-------
MPFR is required for Adaptb. On OpenSuse, MPFR can be installed with::

    > sudo zypper install libmpfr4 mpfr-devel

Python packages in Fedora
-------------------------
Fedora comes with python packages for numpy. You also need 
the setuptools package to be able to install the other python
packages. To install these, do::

    > sudo zypper install python-numpy \
        python-setuptools python-setuptools-devel

Python packages with easy_install
---------------------------------
The  ``nose``, ``mpi4py``, ``h5py`` and ``docutils`` can be 
installed with the ``easy_install`` command::

    > sudo easy_install nose
    > sudo easy_install mpi4py
    > sudo easy_install h5py
    > sudo easy_install docutils 
