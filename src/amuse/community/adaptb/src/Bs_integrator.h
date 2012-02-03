#include "Cluster.h"

#ifndef __BS_INTEGRATOR_H
#define __BS_INTEGRATOR_H

class Bs_integrator
{
  mpreal epsilon;
  int n_max;
  int k_max; 
 
  int flag;

  public:

  // Constructors
  Bs_integrator();
  Bs_integrator(mpreal e);
  Bs_integrator(mpreal e, int n, int k);

  // Set
  void set_epsilon(mpreal e);
  void set_n_max(int n);
  void set_k_max(int k);

  // Get
  mpreal get_epsilon();
  int get_n_max();
  int get_k_max();

  // Functions
  void integrate(Cluster &Cluster, mpreal &dt);
  void step(Cluster &Cluster, mpreal dt);
  void extrapol( Cluster &cl_exp, vector<mpreal> &h, vector<Cluster> &cl );
  mpreal extrapolate(vector<mpreal> h, vector<mpreal> A, mpreal x);
  int error_control(Cluster &cl0, Cluster &cl);
  int converged();
};

#endif


