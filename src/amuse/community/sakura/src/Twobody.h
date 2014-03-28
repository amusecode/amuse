#include <iostream>
using namespace std;

#include <cmath>
#include <cstdlib>

#ifndef __TWOBODY_H
#define __TWOBODY_H

class Twobody {

  double tolerance;

  public:

  Twobody();
  Twobody(double tolerance);

  void set_tolerance(double tolerance);
  double get_tolerance();

  bool solve(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt);
  bool solve_by_leapfrog(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt);
 
  bool try_solve(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt);

  bool calc_universal_anomaly_laguerre(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x);
  bool calc_universal_anomaly_newton(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x);
  bool calc_universal_anomaly_halley(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x);
  bool calc_universal_anomaly_chebyshev(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x);

  double calc_C(double &z);
  double calc_S(double &z);

  void normalize(double &mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double &dt, double &Cm, double &Cr, double &Cv);
  void denormalize(double &mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double &Cm, double &Cr, double &Cv);
  bool is_valid_number(double &x);
};

#endif


