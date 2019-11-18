Installing on macOS High Sierra with Python 2.7
===============================================

In this section we assume a default macOS installation with macports installed.

Installing prerequisites
------------------------
These prerequisites are essential for building AMUSE and the community codes.
They can be installed system-wide with the commands below.
You can choose between openmpi and mpich as desired, both work with AMUSE.
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.

For openmpi do::

  > sudo port install gcc7 openmpi-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python27 py27-virtualenv
  > sudo port select --set mpi openmpi-gcc7-fortran
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python2 python27
  > sudo port select --set virtualenv virtualenv27
  
For mpich do::
	
  > sudo port install gcc7 mpich-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python27 py27-virtualenv
  > sudo port select --set mpi mpich-gcc7
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python2 python27
  > sudo port select --set virtualenv virtualenv27

.. note:
  Please make sure not to install mpich and openmpi together. 
  When both are installed strange errors will occur and AMUSE will not work.
  If you have both installed please first remove both and then install one.
  
Installing AMUSE
----------------

First, create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don't need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory)::

  > virtualenv Amuse-env
  
When the environment is created, you can activate it with::

  > . Amuse-env/bin/activate

You may want to make an alias for this, e.g.::

  > alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'
  
From this point, your prompt will have 'Amuse-env' in front of it, so you will always know when you're in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE::

  > pip install numpy nose docutils mpi4py h5py
  
Probably, you'll want to install these Python modules too::

  > pip install scipy astropy jupyter pandas seaborn
  
Now we can finally install AMUSE itself.
First, download AMUSE or preferably make a git clone (in a desired directory)::

  > git clone https://github.com/amusecode/amuse.git

Then, change to the AMUSE directory and run configure, enabling optional GPU features if present/required::

  > cd amuse
  > ./configure [--enable-cuda] [--enable-sapporo2]

Finally, build and install AMUSE, with optionally downloaded codes if desired::

  > [export DOWNLOAD_CODES=1]
  > python setup.py install
  
Optionally, to test if your setup was successful, run (this will take a long time)::

  > python setup.py test
