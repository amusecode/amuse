using namespace std;
#include <cmath>
#include <vector>

#include "mpreal.h"
using namespace mpfr;

#ifndef __Star_H
#define __Star_H

class Star
{
  public:

  // Variables
  int id;
  mpreal radius;
  mpreal m, x, y, z, vx, vy, vz;
  mpreal ax, ay, az, ax0, ay0, az0;
  mpreal r2_mag, v2_mag, a2_mag;
  
  mpreal dt;

  // Constructors
  Star();
  Star(mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ);
  Star(int id, mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ);
  Star(int id, mpreal M, mpreal R, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ);

  // Reset
  void reset_a();
  void reset_dt();

  // Get
  mpreal get_a2mag();
  mpreal get_v2mag();
  mpreal get_r2mag();

  // Add
  void add_x(mpreal X);
  void add_y(mpreal Y);
  void add_z(mpreal Z);
  void add_vx(mpreal VX);
  void add_vy(mpreal VY);
  void add_vz(mpreal VZ);
  void add_ax(mpreal AX);
  void add_ay(mpreal AY);
  void add_az(mpreal AZ);

  // Calculate
  void calc_r2mag();
  void calc_v2mag();
  void calc_a2mag();
};

#endif


