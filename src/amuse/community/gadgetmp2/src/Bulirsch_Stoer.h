/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Brutus/Bulirsch_Stoer in the version used in Amuse             ### */
/* ###   Author: Brutus and Amuse contributors                                    ### */
/* ###                                                                            ### */
/* ###   Modified: September 2020                                                 ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */

//#include "Star.h"
//#include "Cluster.h"
#include "allvars.hpp"

#ifndef __Bulirsch_Stoer_h
#define __Bulirsch_Stoer_h

class Bulirsch_Stoer {
  my_float tolerance;
  int n_max, k_max;

  public:

  Bulirsch_Stoer();
  Bulirsch_Stoer(my_float tolerance);
  Bulirsch_Stoer(my_float tolerance, int n_max, int k_max);

  inline void set_tolerance(my_float tolerance);
  inline void set_n_max(int n_max);
  inline void set_k_max(int k_max);

  inline my_float get_tolerance();
  inline int get_n_max();
  inline int get_k_max();

  bool integrate(Cluster &cl, my_float &dt);
  bool step(Cluster &cl, my_float &dt);
  void extrapol(Cluster &cl_exp, vector<my_float> &dt, vector<Cluster> &c);
  my_float extrapolate(vector<my_float> x, vector<my_float> y, my_float x0);
  bool error_control(Cluster &c1, Cluster &c2);
};

#endif


