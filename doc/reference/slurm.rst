=============================================
Running AMUSE through slurm
=============================================

Many supercomputers use slurm (see
https://slurm.schedmd.com/documentation.html) as a job-submission
environment.  The Cartesius supercomputer, for example, uses slurm,
but also the Leiden university supercomputer ALICE
https://wiki.alice.universiteitleiden.nl/index.php?title=ALICE_User_Documentation_Wiki&_ga=2.254956621.431957536.1583247665-1261213469.1568893554

The parallel AMUSE script
-------------------------

Slurm operates via batch scripts which have to be written and
subsequently submitted via slurm to the job scheduler.
We use a simple amuse example code as an example.
.. code-block:: sh

	> cd ${AMUSE_DIR}/examples/tutorial/parallel_amuse_script.py

The slurm batch script
----------------------
	
The slurm script to run this code with 6 processors can be callded run_example.sh

.. code-block:: sh

        #!/bin/bash
	sbatch <<EOT
	#!/bin/sh
	#SBATCH --mem=1000
	#SBATCH -p cpu-short
	#SBATCH -N 1
	#SBATCH -n 6

	export OMP_NUM_THREADS=6
	export OMPI_MCA_rmaps_base_oversubscribe=yes 
	export OMPI_MCA_mpi_warn_on_fork=0
	export OMPI_MCA_rmaps_base_oversubscribe=yes

	module load AMUSE/12.0.0-foss-2018a-Python-2.7.14
	module load openmpi/gcc/64/1.10.7
	module load matplotlib/3.1.1-foss-2019b-Python-3.7.4

	mpiexec -n 6 python -u parallel_amuse_script.py
	EOT

The various ingredients of this script identifies the amount of memory
requested, the queue (here called cpu-short), the number of nodes and
processors.

In the next block the environment variable are set.

Then the various modules are loaded.

And finally, the mpiexec command starts the python scrip on 6
processors and runs the script called parallel_amuse_script.py


Running the script
------------------

A slurm script is generally started with some command line, for example:	

.. code-block:: sh

	> sbatch run_example.sh

The output can subsequently be downloaded via scp	

.. code-block:: sh

	> scp user@remote_computer_address:file .
