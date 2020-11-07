#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric>
#include <cstdlib>

#include "Star.h"

#ifndef __Cluster_h
#define __Cluster_h

enum{
    arg_m,
    arg_r_0,
    arg_r_1,
    arg_r_2,
    arg_v_0,
    arg_v_1,
    arg_v_2,
    #ifdef use_additional_acc
    arg_a_step_0,
    arg_a_step_1,
    arg_a_step_2,
    #endif // use_additional_acc
    arg_cnt,
};

class Cluster : public Star {
  public:

  vector<Star> s;
  mpreal eps2;
//  mpreal time, dt, dt_last;
  mpreal  dt;

  Cluster() : Star() {}

  Cluster(vector<double> data);
  Cluster(vector<mpreal> data);

//  vector<double> get_data_double();
  vector<mpreal> get_data();

  void calcAcceleration_dt();
  void calcAcceleration();
  #ifdef use_additional_acc
  void remove_step_Acceleration();
  #endif // use_additional_acc


  void updatePositions(mpreal dt);
  void updateVelocities(mpreal dt);

  void step(mpreal &dt);

//  vector<mpreal> energies();

  friend ostream & operator << (ostream &so, Cluster &cl) {
    for (vector<Star>::iterator si = cl.s.begin(); si != cl.s.end(); ++si) {
      so << *si;
    }
    return so;
  }

};

#endif


