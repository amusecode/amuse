.. _prerequisite-label:


Installation of the prerequisite software
=========================================


.. toctree::
   :maxdepth: 1
   
   install-prerequisites-ubuntu
   install-prerequisites-osx
   install-prerequisites-arch
   install-prerequisites-fedora
   install-prerequisites-centos
   install-prerequisites-suse

   

Before installing AMUSE several software packages must be installed. These 
software packages can be installed manually or with two prepared installation
scripts. The installation scripts will install python and the 
other prerequisites in a user directory. No "root" access is required.

These are the packages AMUSE needs:

* Python (version >= 2.6)
* Numpy (version >= 1.3.0)
* HDF (version 1.6.5 - 1.8.x)
* h5py (version >= 1.2.0)
* MPI (OpenMPI or MPICH2)
* mpi4py (version >= 1.0)
* nose (version >= 0.11)
* docutils (version >= 0.6)
* FFTW (version >= 3.0)
* GSL
* CMake (version >= 2.4)
* GMP (version >= 4.2.1)
* MPFR (version >= 2.3.1)

In the first two sections (compilers_ and installation_scripts_) we explain how to use the two
installation scripts to install AMUSE. In the last section (manual_) 
we have specified the required packages with the needed version for each.

.. _compilers:

Compilers
*********

To build AMUSE from source you need to have a working  build environment.
The AMUSE build system needs a C++ and fortan 90 compiler. Please check first if you
have a working build environment on your system.

In Ubuntu you can setup the environment with (as root):

.. code-block:: sh

	apt-get install build-essential curl g++ gfortran gettext zlib1g-dev



In Fedora you can setup the environment with (as root):

.. code-block:: sh

	yum groupinstall "Development Tools" "Development Libraries"

.. _installation_scripts:

Installation scripts
~~~~~~~~~~~~~~~~~~~~

We have created two installation scripts to automate the installation of
the required packages on a LINUX and OS.X system. These scripts will
install these packages in a user directory. One script downloads and
installs python while the other script downloads and installs the libraries
and python packages. As everything is installed in a user directory these
packages can be installed even if a version of the software is already 
installed on your system. 

The scripts will download and install the software in a user directory. This
user directory must be specified with the ``PREFIX`` environment variable. Before
running the installation scripts you must set the ``PREFIX`` environment 
variable and update the path and library path. For shell (bash) you need to do:

.. code-block:: sh

	export PREFIX=~/amuse/prerequisites
	export PATH=${PREFIX}/bin:${PATH}
  	export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}


One script will download, build and install python on your system. The other 
script is written in Python and will download and install the other packages. 
Both scripts can be found in the ``doc/install`` directory. 

To start the installation do:

.. code-block:: sh

	# 1. Open a shell and go to the <doc/install> directory
	>

	# 2. Set the PREFIX, PATH and LD_LIBRARY_PATH environment variables:
  	> export PREFIX=~/amuse/prerequisites
  	> export PATH=${PREFIX}/bin:${PATH}
  	> export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

	# 3. Start the installation script for python
	> ./install-python.sh

	# 4. Start the installation script for the prerequisite packages
	> ./install.py download
  	> ./install.py install

	# 5. Update your PATH variable in your profile. 
	# Make sure the `${PREFIX}/bin` directory is the first entry in the PATH!

You should now be able to install AMUSE.

Using the installation scripts on OS X
--------------------------------------

For OS.X you need to install XCode and a gfortran compiler first.
The XCode development package is available on the 
`Apple developers site <http://developer.apple.com/devcenter/mac>`_
or for Lion on the Apple Store application.

The standard XCode release does not come with a gfortran compiler. 
Go to the `HPC Mac OS X site <http://hpc.sourceforge.net/index.php>`_ 
for a recent gfortran compiler, compatible with the XCode tools.

After installing XCode and gfortan, follow the steps described in the
previous paragraph.

.. _manual:

Manually installing the prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python
------
Python is probably already installed on your system. To check the version of python do:

.. code-block:: sh

	> python --version
	Python 2.6.2

You can download python from http://www.python.org. 

Numpy
-----
To check if numpy is installed on your system do:

.. code-block:: sh

	> python -c 'import numpy; print numpy.version.version'
	1.3.0

If this fails with an error or a version before 1.3 you need to install numpy.
You can download numpy from http://www.scipy.org/NumPy. 

HDF5 library
------------
HDF5 is a data format specification. The HDF group provides a C library 
to write and access HDF files.

To check if the HDF library is installed on your system do:

.. code-block:: sh

	> h5ls -V
	h5ls: Version 1.8.3

If this fails with an error or a version before 1.6.5 you 
need to install the HDF library. 
You can download HDF from http://www.hdfgroup.org/. 

h5py
----
To access HDF5 files from python we use the ``h5py`` library.

To check if the h5py library is installed on your system do:

.. code-block:: sh

	> python -c 'import h5py; print h5py.version.version'
	1.2.0

If this fails with an error or a version before 1.2.0 you need to install h5py.
You can download h5py from http://code.google.com/p/h5py/. 

docutils
--------
To check if the python docutils are installed on your system do:

.. code-block:: sh

        > python -c 'import docutils; print docutils.__version__'
        0.6

If this fails with an error or a version before 0.6 you need to install docutils.
You can download docutils from http://docutils.sourceforge.net/

MPI
---
The installed MPI framework must be MPI 2 compatible. AMUSE will work with 
MPICH2 or OpenMPI

MPICH2
^^^^^^
MPICH2 is a portable implementation of the MPI 2 standard.

To check if MPICH2 is installed on your system do:

.. code-block:: sh
    
    > mpdhelp
    
    The following mpd commands are available.  For usage of any specific one,
    invoke it with the single argument --help .

    mpd           start an mpd daemon
    mpdtrace      show all mpd's in ring
    mpdboot       start a ring of daemons all at once
    mpdringtest   test how long it takes 
    ...
    
If this fails with an error you need to install MPICH2 or check for OpenMPI
support.
You can download MPICH2 from http://www.mcs.anl.gov/research/projects/mpich2/.

OpenMPI
^^^^^^^
OpenMPI is another portable implementation of the MPI 2 standard

To check if OpenMPI is installed on your system do:

.. code-block:: sh
    
    > mpicxx -v 
    

If this fails with an error you need to install MPICH2 or OpenMPI
support. Most examples in the dopcumentation assume OpenMPI.
You can download OpenMPI from http://www.open-mpi.org/.

MPI4PY
------
To access MPI from python we use the ``mpi4py`` software.
To check if the mpi4py library is installed on your system do:

.. code-block:: sh

	> python -c 'import mpi4py; print mpi4py.__version__'
	1.0.0

If this fails with an error or a version before 1.0 you need to install mpi4py.
You can download mpi4py from http://code.google.com/p/mpi4py/. 

Nose
----
Nose is an extension of the python testing framework. It is used for all
unit testing in AMUSE.


To check if Nose is installed on your system do:

.. code-block:: sh
    
    > nosetests --version
    nosetests version 0.11.1
    ...
    
If this fails with an error or a version before 0.11 you need to install nose.
You can download nose from http://somethingaboutorange.com/mrl/projects/nose/. 

FFTW
----
FFTW is a C subroutine library for computing discrete Fourier transforms. To 
check for the availability of fftw on your system, you can use ``fftw-wisdom``:

.. code-block:: sh

   > fftw-wisdom --version
   fftw-wisdom tool for FFTW version 3.2.1.


You can download the FFTW library from http://www.fftw.org. 

GSL
-------
The GNU Scientific Library (GSL) is a numerical library for C and C++ 
programmers. It is free software under the GNU General Public License.
To check for the availability of GSL on your system, you can use ``gsl-config``:

.. code-block:: sh

   > gsl-config --version
   1.14


You can download GSL from http://www.gnu.org/software/gsl/. 

CMake
-------
CMake is a cross-platform, open-source build system. CMake is used to control 
the software compilation process using simple platform and compiler independent
configuration files. CMake generates native makefiles and workspaces that can 
be used in the compiler environment of your choice.
CMake is used to build EVTwin. 
To check whether you have CMake installed on your system:

.. code-block:: sh

   > cmake --version
   cmake version 2.8.2


You can download CMake from http://www.cmake.org/cmake/resources/software.html. 

GMP
-------
GNU MP is a library for arbitrary precision arithmetic (ie, a bignum package). 
It can operate on signed integer, rational, and floating point numeric types.
GMP is required for Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda 
Boekholt). 
The best way to check whether you have the right version of GMP installed on your 
system depends on the package manager you use, but this should always work (note 
that the library numbers do not match the release version):

.. code-block:: sh

   > locate libgmp
   /usr/lib64/libgmp.so
   /usr/lib64/libgmp.so.10
   /usr/lib64/libgmp.so.10.0.3
   
   > locate gmp.h
   /usr/include/gmp.h
   
   > grep GNU_MP_VERSION /usr/include/gmp.h
   #define __GNU_MP_VERSION 5
   #define __GNU_MP_VERSION_MINOR 0
   #define __GNU_MP_VERSION_PATCHLEVEL 3


You can download GMP from http://gmplib.org. 

MPFR
-------
The MPFR library is a C library for multiple-precision floating-point 
computations with correct rounding.
MPFR is required for Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda 
Boekholt). 
The best way to check whether you have the right version of MPFR installed on your 
system depends on the package manager you use, but this should always work (note 
that the library numbers do not match the release version):

.. code-block:: sh

   > locate libmpfr
   /usr/lib64/libmpfr.so
   /usr/lib64/libmpfr.so.4
   /usr/lib64/libmpfr.so.4.1.0
   
   > locate mpfr.h
   /usr/include/mpfr.h
   
   > grep MPFR_VERSION /usr/include/mpfr.h
   #define MPFR_VERSION_MAJOR 3
   #define MPFR_VERSION_MINOR 1
   #define MPFR_VERSION_PATCHLEVEL 0

You can download MPFR from http://www.mpfr.org. 





    
    





    
    

    
author: Arjen van Elteren (vanelteren@strw.leidenuniv.nl)
date: 2010/09/22





