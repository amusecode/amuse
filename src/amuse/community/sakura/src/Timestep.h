using namespace std;
#include <vector>
#include <cmath>
#include <algorithm>

#include "Particle.h"
#include "Particles.h"

#ifndef __TIMESTEP_H
#define __TIMESTEP_H

class Timestep {

  int mode;
  double dt_max;
  double dt_param;

  public:

  Timestep();
  Timestep(double dt);
  Timestep(int mode);
  Timestep(int mode, double dt_max, double dt_param);

  void set_mode(int mode);
  void set_dt_max(double dt_max);
  void set_dt_param(double dt_param);
  void set_dt(double dt);

  int get_mode();
  double get_dt_max();
  double get_dt_param();
  double get_dt();

  double get_dt(vector<Particle> &particle);

  double get_constant_dt();
};

#endif


