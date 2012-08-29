-------------------------------- HiGPUs readme ----------------------------------

This is HiGPUs readme. (If there are any problems please do not 
hesitate to contact us:
Mario Spera:  mario.spera@uniroma1.it
Davide Punzo: Davide.Punzo@roma1.infn.it
Roberto Capuzzo Dolcetta : roberto.capuzzodolcetta@uniroma1.it).

-------------------------------- Introduction -----------------------------------

HiGPUs is a parallel direct N-body code based on a 6th order Hermite 
integrator. The code uses, at the same time, MPI, OpenMP and CUDA 
libraries to fully exploit all the capabilities offered by the 
most modern hybrid supercomputers in the world. Moreover it is 
implemented using block time steps taking into account the 
properties of the stiff problems such as the gravitational N-body 
problem.

If you use HiGPUs for scientific work, we kindly ask you to 
reference the following paper: "A fully parallel, high precision, N-body code running on
hybrid computing platforms" by R. Capuzzo–Dolcetta, M. Spera, D. Punzo (2012, submitted to JCP) ,
available in arXiv: 2012arXiv1207.2367C.

Feel free to suggest modifications, imporvements, additions, etc. to the code by contacting the authors.

---------------------------------- Installing and Compiling -------------------------------

First of all you should check if CUDA and MPI libraries are 
installed on your machine. The code is known to work with CUDA 
4.0 and OpenMPI 1.5.4 using the following kind of devices:

- Tesla C1060.

- Tesla C2050.

- Tesla M2070. 

- GeForce GTX 480.

- Geforce GTX 580.

In any case the GPU must support double precision floating point 
operations. If you try our code on other hardware please inform 
us about the results. If you need to install CUDA go to 
http://developer.nvidia.com/category/zone/cuda-zone and follow the 
instructions in the “Getting Started Guide” in the documentation 
section. To get OpenMPI, you can download the package at 
http://www.open-mpi.org/software/ompi/v1.5/ and follow the 
instructions contained in the readme file.

To use the code untar the downloaded archive "HiGPUs.tar" with the 
command line

> tar -zxvf HiGPUs.tgz

This will create the directory 'HiGPUs' which contains the source 
and all other files needed for running the code. The content of 
'HiGPUs' should look like:

> docs 

> exec 

> lib

> Makefile

> src

At this point, before compiling the code, you have to tell the
Makefile where to find CUDA and MPI libraries.
To do this, you have to modify the following two lines of the
Makefile setting appropriately the variables:

> CUDA_INSTALL_PATH := /usr/local/cuda

> MPI_INSTALL_PATH := /usr/openmpi

where we have assumed that the installation paths are :
/usr/local/cuda for CUDA, /opt/AMDAPP for the AMD OpenCL implementation
and /usr/openmpi for OpenMPI.
Before compiling the code you can add the following flags to the
variables MYOPTS:

-DCHECK_ERRORS : some useful functions which may help to find
possible errors activated. If this option is enabled, the program
terminates its esecution telling at what line the error has been
detected and giving further information about its type.

-DCHECK_TIMES : some functions which measure the time spent in
executing the main sections of the code are activated. The
results are showed in the output file 'times.dat'.

-DPLUMMER : (STILL EXPERIMENTAL) it enables an external Plummer-like
potential field which adds an analytical contribution to the mutual
interaction between the stars. To set the related parameters (core
radius and total mass) see below.

Once you have setup the Makefile you can proceed typing:

> make clean 

> make (note: in the AMUSE release the command to compile without 
the interface is 'make exec')

This creates the executable file 'HiGPUs.x' in the folder 'exec'. 

------------------------------------- Parameter and input files -------------------------------------------------

The file containing the parameters is called 'input_param.txt' 
and it looks like:

> 262144 // number of particles ! 

> 2 // number of gpus to use (obviously we reccomend all gpu 
which are avaible on singole node) !

> 128 // gpus threads per block (from our benchmark we obtain the 
best performance with 128) ! 

> 1.0 // integration time !

> -3.0 // exponent which defines the maximum time step allowed 
for particles (2^exponent) ! 

> -30.0 // exponent which defines the minimum time step allowed 
for particles (2^exponent) !

> 0.45 // eta parameter for determining particles time steps 
(generalized Aarseth criterion for the Hermite 6th order method) 
! 

> 0.01 // eta parameter for initializing blocks (Aarseth 
criterion) !

> 0.125 // time step for snapshots ! 

> 1000000 // maximum number of snapshots ! 

> 0 // rescale to the center of mass 1 = true, 0 = false ! 

> 0 // rescale to the center of velocity 1 = true, 0 = false !

> input_N_262144.dat // input file for positions, velocities and 
masses (it must be created in THIS order) ! 

> Tesla C2050# // GPU to use (the name must be terminated with 
the character #) .

This file must be included in the same folder you want to run the 
simulation. 

The input data file, whose name is specified in the parameters 
file, must contain the N-body data for the simulation in the 
following format :

> x y z vx vy vz mass softening

An example of input data file is in the folder 'exec'.

NOTE: if you change parameters or input file, it is not necessary 
to compile again the code .

--------------------------------- Running---------------------------------------

Now you are ready to run the program. If you are using a machine 
with one node you can simply type:

> mpirun -np 1 ./HiGPUs.x

Otherwise you have to write your own hostfile which, in case of 
two nodes is something like:

>node01 slots=1 

>node02 slots=1

where the first string of each line is the name of the node and 
the variable 'slots' indicates how many MPI processes you want to 
run on the related node. Since in our program we established a 
one to one correspondence between MPI processes and computational 
nodes, it is strongly reccomended to run only one MPI process per 
node. Each MPI process can manage all the GPUs avaible on each 
node without problems. At this point, you can use the following 
command to start the simulation:

> mpiexec -hostfile myhostfile./HiGPUs.x

If you run the program with '-h' option you can obtain further 
useful information:

> ------------------------------------------------------- 

> ----------------------Help of HiGPUs---------------

> -------------------------------------------------------

> Usage : ./HiGPUs [options] 

> Options : 

> -f [file] (input_param.txt) : it specifies the file containing 
the simulation parameters. If not specified, the default is 
'input_param.txt'.

> -h : it shows this help screen 

> -r [file] : it restarts the simulation from the specified file 
(H6Blog.dat is necessary, see below). 

> -p b=[value]:M=[value] : if you have enabled the -DPLUMMER 
option, it is necessary to specify the parameters b (core radius) 
and M (total mass).

> 

> 

>NOTE: Plummer potential \varphi(r)=\frac{G*M}{\sqrt{r^{2}+b^{2}}}

>

------------------------------ Output ----------------------------------------------

The following file will be generated or updated every time a new 
snapshot is produced:

- (MAX + i).dat : it contains positions, velocities and masses 
(in this order) at the i-th snapshot. The value of MAX is 
specified in the line 'maximum number of snapshots !' of the 
'input_param.txt'. 

- cpu_memory.dat : this file reports, every time the function 
'CPU_memcheck' is called, the ram memory used by the program at 
that point.

- gpu_memory.dat : this file reports, every time the function 
'GPU_memcheck' is called, the GPU on-board memory used by the 
program at that point.

- Blocks.dat : it reports the histogram of the block time steps 
scheme. That is in the first column the base-2 logarithm of the 
time step is reported while the second column is the number of 
particles having the same time step.

- energy.dat : it shows the simulation global time versus the 
ratio (E-E0)/E0 where E is the total energy at the current time 
and E0 is the total initial energy.

- HiGPUs.dat : it summarizes the simulation information. It could 
be something like: 

> ==================================================

> Read parameters file : input_param.txt

> N : 262144 

> Gpus : 2

> Threads per block : 128 

> Time of integration : 1.0 

> Max time step : 0.125 

> Min time step : 9.31323e-10

> eta 6th order : 0.45 

> eta 4th order : 0.01 

> time for printing : 1.0

> Max output files : 1000000

> CDM scale : 0

> CDV scale : 0

> Input file : input_N_262144 

> ==================================================

> Available : Tesla C2050 as device : 0 

> Available : Quadro as device : 1 

> Available : Tesla C2050 as device : 2 

> ============================================= 

> Using : Tesla C2050 (device 0)

> Using : Tesla C2050 (device 2)

> =============================================== 

> #Initial Total Energy : #-3.00234523435865e-3#

- times.dat : it is created only if you enable -DCHECK_TIMES 
option. It has a maximum of N data lines and reports :

- N : number of particles updated. 

- NEXT : time required to choose the particles that have to be 
updated. These particles are stored in the array named 'next'.

- PRED : time required to copy the array 'next' from the Host to 
the GPU and to perform the predictor step. 

- EVAL : time required to calculate the accelerations of the 
particles contained in the array 'next'. 

- REDU : time required to reduce forces. 

- REPOS : time required to collect accelerations. 

- CPY_ACC : time required to copy accelerations from the GPU to 
the Host. 

- MPI : time required to perform the MPI_Allreduce() function to 
reduce and collect forces from all the computational nodes. 

- CORR : time required to perform the corrector. 

- RECON : time required to reorder data on the GPUs after the 
corrector step. 


