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

  (will be expanded)..
