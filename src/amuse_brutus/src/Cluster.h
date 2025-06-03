#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>

#include "Star.h"

#ifndef __Cluster_h
#define __Cluster_h

class Cluster : public Star { 
  public:

  std::vector<Star> s;
  mpreal eps2;
  mpreal time, dt, dt_last;

  Cluster() : Star() {}

  Cluster(std::vector<double> data);
  Cluster(std::vector<mpreal> data);

  std::vector<double> get_data_double();
  std::vector<mpreal> get_data();    

  void calcAcceleration_dt();
  void calcAcceleration();

  void updatePositions(mpreal dt);
  void updateVelocities(mpreal dt);

  void step(mpreal &dt);
  
  std::vector<mpreal> energies();

  friend std::ostream & operator << (std::ostream &so, Cluster &cl) {
    for (std::vector<Star>::iterator si = cl.s.begin(); si != cl.s.end(); ++si) {
      so << *si;
    }
    return so;
  }
  
};

#endif


