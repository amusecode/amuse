#include <cmath>

#ifndef __PARTICLE_H
#define __PARTICLE_H

class Particle {
  int index;

  double m, R;

  double x, y, z;
  double vx, vy, vz;
  
  double ax, ay, az;
  double jx, jy, jz;
  double sx, sy, sz;
  double cx, cy, cz;
 
  public:

  Particle();
  Particle(double m, double x, double y, double z, double vx, double vy, double vz);
  Particle(int index, double m, double x, double y, double z, double vx, double vy, double vz);
  Particle(int index, double m, double R, double x, double y, double z, double vx, double vy, double vz);

  void set_index(int index);
  void set_mass(double m);
  void set_radius(double R);
  void set_x(double x);
  void set_y(double y);
  void set_z(double z);
  void set_vx(double vx);
  void set_vy(double vy);
  void set_vz(double vz);
  void set_ax(double ax);
  void set_ay(double ay);
  void set_az(double az);
  void set_jx(double jx);
  void set_jy(double jy);
  void set_jz(double jz);
  void set_sx(double sx);
  void set_sy(double sy);
  void set_sz(double sz);
  void set_cx(double cx);
  void set_cy(double cy);
  void set_cz(double cz);

  int get_index();
  double get_mass();
  double get_radius();
  double get_x();
  double get_y();
  double get_z();
  double get_vx();
  double get_vy();
  double get_vz();
  double get_ax();
  double get_ay();
  double get_az();
  double get_jx();
  double get_jy();
  double get_jz();
  double get_sx();
  double get_sy();
  double get_sz();
  double get_cx();
  double get_cy();
  double get_cz();
  double get_r2();
  double get_v2();
  double get_a2();
  double get_j2();
  double get_s2();
  double get_c2();
  double get_r();
  double get_v();
  double get_a();
  double get_j();
  double get_s();
  double get_c();

  void add_to_mass(double dm);
  void add_to_radius(double dR);
  void add_to_x(double dx);
  void add_to_y(double dy);
  void add_to_z(double dz);
  void add_to_vx(double dvx);
  void add_to_vy(double dvy);
  void add_to_vz(double dvz);
  void add_to_ax(double dax);
  void add_to_ay(double day);
  void add_to_az(double daz);
  void add_to_jx(double djx);
  void add_to_jy(double djy);
  void add_to_jz(double djz);
  void add_to_sx(double dsx);
  void add_to_sy(double dsy);
  void add_to_sz(double dsz);
  void add_to_cx(double dcx);
  void add_to_cy(double dcy);
  void add_to_cz(double dcz);
};

#endif


