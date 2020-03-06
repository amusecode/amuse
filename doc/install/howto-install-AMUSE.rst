Installing the prerequisites
============================

For a full AMUSE installation, you will need to install some further dependencies that can be installed via your package manager - e.g. apt or yum on Linux; macports or homebrew on macOS.

Ubuntu
******

You can choose between openmpi and mpich as desired, both work with AMUSE. Please do not install both!
In the examples below we choose GCC-7 as the compiler, but more recent versions of GCC will also work.

* For openmpi:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      libopenmpi-dev openmpi-bin \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      git

* For mpich:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      mpich libmpich-dev \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      git


macOS
*****


On macOS, you will first need to install Xcode. You can do so via the app store.

In this section we assume a default macOS installation (up to Catalina) with MacPorts, but other methods (such as Homebrew) will also work.

You can choose between openmpi and mpich as desired, both work with AMUSE. 
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.
In the examples below we choose GCC 9 as the compiler, but other versions of GCC should also work.

* For openmpi:

.. code-block:: sh

    sudo port install gcc9 openmpi-gcc9 hdf5 gsl cmake gmp mpfr fftw-3 +gcc9
    sudo port install python38
    sudo port select --set mpi openmpi-gcc9-fortran
    sudo port select --set gcc mp-gcc9
    sudo port select --set python3 python38

* For mpich:

.. code-block:: sh

    sudo port install gcc9 mpich-gcc9 hdf5 gsl cmake gmp mpfr fftw-3 +gcc9
    sudo port install python38
    sudo port select --set mpi mpich-gcc9
    sudo port select --set gcc mp-gcc9
    sudo port select --set python3 python38



Installing AMUSE
================


After installing the prerequisites, you can install AMUSE.
Optionally, first create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don’t need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory):

.. code-block:: sh

    python3 -m venv Amuse-env

When the environment is created, you can activate it with:

.. code-block:: sh

    . Amuse-env/bin/activate

You may want to make an alias for this, e.g.:

.. code-block:: sh

    alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'

From this point, your prompt will have ‘Amuse-env’ in front of it, so you will always know when you’re in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE:

.. code-block:: sh

    pip install numpy nose docutils mpi4py h5py wheel

Probably, you’ll want to install these Python modules too:

.. code-block:: sh

    pip install scipy astropy jupyter pandas seaborn

Now we can finally install AMUSE itself.
This is done easiest via pip:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse

If you only require a subset of AMUSE, you can install any of the individual packages as such:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse-$(community_code_name)

