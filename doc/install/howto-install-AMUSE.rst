Installing the prerequisites
============================

For a full AMUSE installation, you will need to install some further dependencies that can be installed via your package manager - e.g. apt or yum on Linux; macports or homebrew on macOS.

Ubuntu
******

You can choose between openmpi and mpich as desired, both work with AMUSE. Please do not install both!

* For openmpi:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      libopenmpi-dev openmpi-bin \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      libblas-dev liblapack-dev \
      python3-venv python3-pip git

* For mpich:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      mpich libmpich-dev \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      libblas-dev liblapack-dev \
      python3-venv python3-pip git


macOS
*****


On macOS, you will first need to install Xcode. You can do so via the app store.
In the examples below we choose GCC-11 as the compiler. Older versions may not work on recent versions of macOS.
In macOS Big Sur and later, you may have to add the following line to your .bashrc or .zshrc profile:

.. code-block:: sh

    export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

In this section we assume a default macOS installation (up to Monterey) with MacPorts, but other methods (such as Homebrew) will also work.
**Please do not install packages using more than one package manager (MacPorts, Homebrew and/or Conda), as this will almost certainly lead to problems!**
On a Mac with an Apple Silicon chip (M1 and later) please use gcc12 (or later), earlier versions of gcc may not always work.

You can choose between openmpi and mpich as desired, both work with AMUSE. 
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.

* For openmpi:

.. code-block:: sh

    sudo port install gcc12 openmpi-gcc12 hdf5 gsl cmake gmp mpfr fftw-3 +gcc12 openblas lapack
    sudo port install python39
    sudo port select --set mpi openmpi-gcc12-fortran
    sudo port select --set gcc mp-gcc12
    sudo port select --set python3 python39

* For mpich:

.. code-block:: sh

    sudo port install gcc12 mpich-gcc12 hdf5 gsl cmake gmp mpfr fftw-3 +gcc12 openblas lapack
    sudo port install python39
    sudo port select --set mpi mpich-gcc12
    sudo port select --set gcc mp-gcc12
    sudo port select --set python3 python39



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

    pip install --upgrade pip

    pip install numpy docutils mpi4py h5py wheel

Probably, you’ll want to install these Python modules too:

.. code-block:: sh

    pip install scipy astropy jupyter pandas seaborn matplotlib

Now we can finally install AMUSE itself.
This is done easiest via pip:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse

If you only require a subset of AMUSE, you can install any of the individual packages as such:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse-$(community_code_name)



Re-installation notes and troubleshooting pip installs
******************************************************

The packages installed with pip are distributed as source packages that must be compiled against the libraries
installed on your local machine. After compilation pip saves a binary package version in its cache.
In case of problems with the AMUSE installation using pip or if the environment changes it may be necessary to clean the pip cache (e.g. at ```~/.cache/pip```). In addition, the cache can be disabled using the ```--no-cache-dir``` option. the ```--no-build-isolation``` may also be tried in case the virtualenv has all the prerequisites, but the build still fails.
The ```--no-clean``` pip install option preserves the build directory for debugging purposes (The actual directory is reported 
in verbose mode ```-v```). 



Development build
*****************

Alternatively, you can install amuse as a development build, which allows you to modify the source code. It is potentially also more convenient when encountering issues with installation of specific codes as the build.log file in the root directory of the repository contains the error logs of the installation process.

Installation can also be handled through pip by executing (in the root of a clone of the repository)

.. code-block:: sh

    pip install -e .

after this the codes need to be build:

.. code-block:: sh

    python setup.py develop_build

individual codes can be build with:

.. code-block:: sh

    make {code}.code

with {code} the name of the code in lower case. 
