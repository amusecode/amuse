Hierarchically split-Up Approximately sYmplectic N-body sOlver 
(HUAYNO)

Inti Pelupessy - january 2011

short description
-----------------

  HUAYNO is a code to solve the astrophysical N-body problem. It uses
  recursive Hamiltonian splitting to generate multiple-timestep integrators
  which conserve momentum to machine precision. A number of different 
  integrators are available. The code has been developed within the 
  AMUSE environment. It can make use of GPUs - for this an OpenCL 
  version can be compiled.

Use
---

  Use of the code is the same as any gravity the code within AMUSE. There are
  three parameters for the code: Smoothing length (squared) eps2, timestep
  parameter eta and a parameter to select the integrator:

  timestep_parameter( eta ): recommended values eta=0.0001 - 0.05
  epsilon_squared (eps2) : eps2 can be zero or non-zero
  inttype_parameter( inttype ): possible values for inttype are described below. 

  Miscellaneous:

  - The code assumes G=1, 
  - Collisions are not implemented (needs rewrite),
  - workermp option can be used for OpenMP parallelization,
  - The floating point precision of the calculations can be changed by setting
    the FLOAT and DOUBLE definitions in evolve.h. FLOAT sets the precision of
    the calculations, DOUBLE the precision of the position and velocity
    reductions. They can be set to e.g. float, double, long double or __float128
    It is advantageous to choose set DOUBLE at a higher precision than FLOAT.
    recommended is the combination: double/ long double  
  - the AMUSE interface uses double precision
  - It is unlikely the integer types in evolve.h would need to be changed
    (INT should be able to hold particle number, LONG should be able to hold
    interaction counts) 

  OpenCL operation:

  By compiling worker_cl it is possible to offload the force and timestep
  loops to the GPU.  The implementation is based on the STDCL library
  (www.browndeertechnology.com) so this library should be compiled first. 
  In the Makefile the corresponding directories should point to the
  installation directory of STDCL. The define in evolce_cl.h should be set
  appropiately for the OpenCL configuration: CLCONTEXT:  stdgpu or stdcpu
  NTHREAD:  64 for GPU, 2 for CPU (actually 64 will also work for CPU just
  fine) BLOCKSIZE: number of particles stored in local memory (64 good
  starting guess) evolve_kern.cl contains the OpenCL kernels. Precision of
  the calculation is controlled by FLOAT(4) defines in evolve_kern.cl and
  CLFLOAT(4) in evolve.h. They should agree with each other (i.e. float and
  cl_float or  double and cl_double)
