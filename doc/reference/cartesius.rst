=============================================
Running AMUSE on The Cartesius supercomputer
=============================================

The Cartesius is the Dutch national supercomputer. For more information on e.g. obtaining an account, see the SURFSara website: https://userinfo.surfsara.nl/systems/cartesius.

Using AMUSE on the Cartesius is relatively straightforward. Below is a tested method for installing and using AMUSE using the prerequisites, though other options (e.g. using only software pre-installed on the machine) should also be possible. For a generic description of the installation of prerequisites, see the :ref:`prerequisite-label` section.

Obtaining AMUSE
---------------

We assume a copy of AMUSE has been downloaded to an `amuse` folder in the users home. For instance using git:

.. code-block:: sh

	> cd /home/USERNAME
	
	> git clone https://github.com/amusecode/amuse.git

Environment settings
--------------------

Since we will be using the AMUSE prerequisites software, we need to set some environment variables. Make sure the following lines are present in your .bashrc file:

.. code-block:: sh

	#file: ~/.bashrc

	#load java 1.7 needed for some codes
	module load java/oracle

	export PREFIX=~/amuse/prerequisites
	export PATH=${PREFIX}/bin:${PATH}
	export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

	#use gfortran  	
	export FC=gfortran
	export F77=gfortran
  	
Also make sure that .bashrc is loaded in your .bash_profile. This way, the enviroment is always set correctly, both in interactive and non-interactive mode.
 
.. code-block:: sh

	#file: ~/.bash_profile

	if [ -f ~/.bashrc ]; then
	    . ~/.bashrc
	fi

**Note: be sure to re-connect to the machine for these changes to take effect.**

Install AMUSE prerequisites
---------------------------

Next, we will install all prerequisites of amuse using the AMUSE supplied scripts.

.. code-block:: sh

	# create a directory for the prerequisites
	> mkdir ~/amuse/prerequisites

	# go to the <doc/install> directory
	> cd ~/amuse/doc/install
	
	# Start the installation script for Python.
	> ./install-python.sh

	# Download the prerequisite packages.
	> ./install.py download
	
	# Install prerequisites. Use hydra as the default MPI process manager. May take a while...
  	> ./install.py --hydra install
  	
  	# Optionally also install matplotlib
  	> ./install.py --matplotlib install
 
Configure and build AMUSE
-------------------------

Configuring and building amuse is now as normal.

.. code-block:: sh

	# go to the amuse directory
	> cd ~/amuse

	# configure amuse
	> ./configure MPIEXEC=mpiexec.hydra
	
	# build amuse
	> make

        # optionally also install codes requiring downloading files
	> make DOWNLOAD_CODES=1
	
	
Test the installation
---------------------

To test your AMUSE installation, run nosetests.

**Note: do not run simulations on the frontend of the cartesius. This is not allowed!**

.. code-block:: sh

	# go to the amuse directory
	> cd ~/amuse
	
	> mpiexec.hydra -n 1 nosetests -v tests

	
Running on a Cartesius node
---------------------------

Running on the Cartesius is typically done by submitting a slurm script. See the surfsare site for more info:

https://userinfo.surfsara.nl/systems/cartesius/usage/batch-usage
 
Below is a simple example script for running amuse on Cartesius.

.. code-block:: sh
	
	#!/bin/bash
	#SBATCH -N 1
	#SBATCH -n 1
	#SBATCH -p short
	#SBATCH -t 10

	cd ~/amuse

	mpiexec.hydra -n 1 ./amuse.sh examples/syllabus/gravity_simple.py
	
Submit using sbatch, get status using squeue, cancel using scancel
	
.. code-block:: sh

	#submit a script
	> sbatch example-script
	
	#list jobs of current user
	> squeue -u $USER

	#cancel job 505224
	scancel 505224

	#cancel all jobs of the current user
	> scancel -u $USER


Using multiple nodes should also work, by specifying this to slurm. MPI will automatically pickup on this and spread workers over all nodes.

.. code-block:: sh
	
	#!/bin/bash
	#SBATCH -N 2
	#SBATCH -n 10
	#SBATCH -p short
	#SBATCH -t 10

	cd ~/amuse

	#note that this simple example uses only a single worker
	mpiexec.hydra -n 1 ./amuse.sh examples/syllabus/gravity_simple.py
	
