#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>

#include "Star.h"

#ifndef __Cluster_h
#define __Cluster_h

class Cluster : public Star { 
  public:

  vector<Star> s;
  mpreal eps2;
  mpreal time, dt, dt_last;

  Cluster() : Star() {}

  Cluster(vector<double> data);
  Cluster(vector<mpreal> data);

  vector<double> get_data_double();
  vector<mpreal> get_data();    

  void calcAcceleration_dt();
  void calcAcceleration();

  void updatePositions(mpreal dt);
  void updateVelocities(mpreal dt);

  void step(mpreal &dt);
  
  vector<mpreal> energies();

  friend ostream & operator << (ostream &so, Cluster &cl) {
    for (vector<Star>::iterator si = cl.s.begin(); si != cl.s.end(); ++si) {
      so << *si;
    }
    return so;
  }
  
};

#endif


