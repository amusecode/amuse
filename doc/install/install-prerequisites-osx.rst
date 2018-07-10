Installing on macOS
*******************

Installing on macOS with MacPorts 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a clean MacPorts installation. The MacPorts build
system will build most packages from source so installation may take a while.
The packages in MacPorts support different *variants*, each *variant* is built
differently. Below, with all installation commands we will specify the variant
where needed. AMUSE is tested with gcc versions 4.3 up to 7. Below, we will use
gcc 5.

.. note::
    
    If you want to use a different fortran compiler (ifort), you are better 
    off using the **install.py** script in the **doc/install** directory.

.. note::

    Make sure you have a recent MacPorts installation. Do this by running 
    
    > sudo port selfupdate

    before installing.
    
.. note::
    
    If you are unsure of your installation you can uninstall and clear the 
    packages with::
    
        port uninstall py27-docutils py27-nose py27-mpi4py py27-h5py py27-numpy hdf5 fftw-3 gsl openmpi python27
    
    To make a clean install of MacPorts, please remove the MacPorts directory
    and read the guide at:
    https://guide.macports.org/
    
All in one
----------

You can install all the packages described below in one go with::

    > sudo port install gcc5
    
    > sudo port install python27
    > sudo port install openmpi-gcc5
    > sudo port install fftw-3 +gcc5
    > sudo port install hdf5 gsl cmake gmp mpfr
    > sudo port install py27-numpy py27-h5py py27-nose py27-docutils 
    > sudo port install py27-mpi4py +openmpi
    > sudo port install py27-matplotlib

To make sure the right MacPorts compilers and python are set as default, do the
following::

    > sudo port select --set mpi openmpi-gcc5-fortran
    > sudo port select --set gcc mp-gcc5
    > sudo port select --set python python27
    > sudo port select --set nosetests nosetests27
    
After installing you will need to configure the code with the following line::

    > ./configure --with-fftw=/opt/local

.. note::

    The ``--with-fftw`` option will ensure that fftw is found in 
    ``/opt/local`` and not in any other location on macOS. Sometimes, incompatible versions of
    fftw may be installed in ``/usr/include`` or ``/usr/local/inlude``. These versions may
    have the wrong processor type (32 vs 64bits) or not contain a fortran API. For both cases
    compiling ``Fi`` will fail.      
    In case the configure script does not pick the wanted fftw directories, you
    can edit the ``config.mk`` file to point to the right version.
    

.. note::

    The ``PREFIX`` variable will make sure some support libraries for 
    community codes (gsl, gmp and mpfr) are found in ``/opt/local``.

.. note::

    Please, make sure you have no other compile variables specified
    (like CC or CXX or CFLAGS), unless you have customized MacPorts in
    any way. Some variable settings will likely give compile errors
    for the community codes. 
    
    For example, BHTree is compiled with openmpicxx and $CXX. 
    The command in the CXX variable must be compatible 
    with openmpicxx (you can do ``openmpicxx --show`` to get 
    the command openmpicxx is using)


GCC
---
By default MacPorts uses the XCode compilers, these compilers have no support
for fortran, a MacPorts gcc compiler set needs to be installed. We suggest
installing gcc 5 as the most reliable option:

.. code-block:: sh
    
    > sudo port install gcc5
    
.. note::
    
    If you have installed a different version of gcc, you need to select
    a different variant of the packages below. To select a different variant
    replace **+gcc5** with **+gcc7**, **+gcc6** or any other version
    matching your gcc installation. Note, apple-gcc versions will not work,
    these do not support fortran.

Python
------
MacPorts supports several python versions in different variants, we will install
the python27 versions

.. code-block:: sh

    > sudo port install python27
   
MPI2
----
MacPorts provides packages for mpich and openmpi. Although you can
probably install both, this is not recommended. We suggest you install
openmpi.

To install openmpi, do::

     > sudo port install openmpi +gcc5

HDF5
----
Amuse can work with HDF5 versions 1.6.*, 1.8.3 and higher. MacPorts comes
with HDF5 version 1.8.* and 1.10.* To install the most recent version, do::

    > sudo port install hdf5 

FFTW-3
------
MacPorts comes with a FFTW and FFTW-3 package, for AMUSE we need FFTW-3.
FFTW-3 can be installed with::

    > sudo port install fftw-3 +gcc5

GSL
---
GSL is used to build Gadget2, GSL can be installed with::

    > sudo port install gsl

CMake
-----
CMake is used to build EVTwin, CMake can be installed with::

    > sudo port install cmake

GMP
-------
GMP is required for Adaptb. With MacPorts, GMP can be installed with::

    > sudo port install gmp

MPFR
-------
MPFR is required for Adaptb. With MacPorts, MPFR can be installed with::

    > sudo port install mpfr


Python packages
---------------
By this point all libraries and frameworks are installed. We can now
install python packages (some depend on the installed libraries)::

    > sudo port install py27-numpy py27-h5py py27-nose py27-docutils

If you installed openmpi in the MPI2 step you need to set the
"openmpi" variant for "py27-mpi4py"::

    > sudo port install py27-mpi4py +openmpi


Matplotlib
----------
Matplotlib is not required but is highly recommended for creating graphics, 
you can install it with::

    > sudo port install py27-matplotlib
    

.. note::

    Macports will install the compilers under non standard names.
    To make sure the right MacPorts compilers and python are set as default, do
    the following::

    > sudo port select --set mpi openmpi-gcc5-fortran
    > sudo port select --set gcc mp-gcc5
    > sudo port select --set python python27
    > sudo port select --set nosetests nosetests27

    Alternatively, to use the right compilers you can specify these during the
    configure stage of AMUSE.
    
    See the output for ``configure --help`` for a list of all 
    environment variables you can set.
    
    If you installed openmpi you need to specify the mpi compilers 
    like so (replacing OPENMPICXX OPENMPICC and OPENMPIF90 with your installed compiler)::
    
        ./configure MPICXX=OPENMPICXX MPICC=OPENMPICC MPIFC=OPENMPIF90

