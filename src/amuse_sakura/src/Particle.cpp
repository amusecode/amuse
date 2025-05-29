#include "Particle.h"

Particle::Particle() {
  index = -1;
  m = 0;
  R = 0;
  x = 0;
  y = 0;
  z = 0;
  vx = 0;
  vy = 0;
  vz = 0;
  ax = 0;
  ay = 0;
  az = 0;
  jx = 0;
  jy = 0;
  jz = 0;
  sx = 0;
  sy = 0;
  sz = 0;
  cx = 0;
  cy = 0;
  cz = 0;
}
Particle::Particle(double m, double x, double y, double z, double vx, double vy, double vz) {
  index = -1;
  this->m = m;
  R = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->vx = vx;
  this->vy = vy;
  this->vz = vz;
  ax = 0;
  ay = 0;
  az = 0;
  jx = 0;
  jy = 0;
  jz = 0;
  sx = 0;
  sy = 0;
  sz = 0;
  cx = 0;
  cy = 0;
  cz = 0;
}
Particle::Particle(int index, double m, double x, double y, double z, double vx, double vy, double vz) {
  this->index = index;
  this->m = m;
  R = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->vx = vx;
  this->vy = vy;
  this->vz = vz;
  ax = 0;
  ay = 0;
  az = 0;
  jx = 0;
  jy = 0;
  jz = 0;
  sx = 0;
  sy = 0;
  sz = 0;
  cx = 0;
  cy = 0;
  cz = 0;
}
Particle::Particle(int index, double m, double R, double x, double y, double z, double vx, double vy, double vz) {
  this->index = index;
  this->m = m;
  this->R = R;
  this->x = x;
  this->y = y;
  this->z = z;
  this->vx = vx;
  this->vy = vy;
  this->vz = vz;
  ax = 0;
  ay = 0;
  az = 0;
  jx = 0;
  jy = 0;
  jz = 0;
  sx = 0;
  sy = 0;
  sz = 0;
  cx = 0;
  cy = 0;
  cz = 0;
}

void Particle::set_index(int index) {
  this->index = index;
}
void Particle::set_mass(double m) {
  this->m = m;
}
void Particle::set_radius(double R) {
  this->R = R;
}
void Particle::set_x(double x) {
  this->x = x;
}
void Particle::set_y(double y) {
  this->y = y;
}
void Particle::set_z(double z) {
  this->z = z;
}
void Particle::set_vx(double vx) {
  this->vx = vx;
}
void Particle::set_vy(double vy) {
  this->vy = vy;
}
void Particle::set_vz(double vz) {
  this->vz = vz;
}
void Particle::set_ax(double ax) {
  this->ax = ax;
}
void Particle::set_ay(double ay) {
  this->ay = ay;
}
void Particle::set_az(double az) {
  this->az = az;
}
void Particle::set_jx(double jx) {
  this->jx = jx;
}
void Particle::set_jy(double jy) {
  this->jy = jy;
}
void Particle::set_jz(double jz) {
  this->jz = jz;
}
void Particle::set_sx(double sx) {
  this->sx = sx;
}
void Particle::set_sy(double sy) {
  this->sy = sy;
}
void Particle::set_sz(double sz) {
  this->sz = sz;
}
void Particle::set_cx(double cx) {
  this->cx = cx;
}
void Particle::set_cy(double cy) {
  this->cy = cy;
}
void Particle::set_cz(double cz) {
  this->cz = cz;
}

int Particle::get_index() {
  return index;
}
double Particle::get_mass() {
  return m;
}
double Particle::get_radius() {
  return R;
}
double Particle::get_x() {
  return x;
}
double Particle::get_y() {
  return y;
}
double Particle::get_z() {
  return z;
}
double Particle::get_vx() {
  return vx;
}
double Particle::get_vy() {
  return vy;
}
double Particle::get_vz() {
  return vz;
}
double Particle::get_ax() {
  return ax;
}
double Particle::get_ay() {
  return ay;
}
double Particle::get_az() {
  return az;
}
double Particle::get_jx() {
  return jx;
}
double Particle::get_jy() {
  return jy;
}
double Particle::get_jz() {
  return jz;
}
double Particle::get_sx() {
  return sx;
}
double Particle::get_sy() {
  return sy;
}
double Particle::get_sz() {
  return sz;
}
double Particle::get_cx() {
  return cx;
}
double Particle::get_cy() {
  return cy;
}
double Particle::get_cz() {
  return cz;
}
double Particle::get_r2() {
  return x*x + y*y + z*z;
}
double Particle::get_v2() {
  return vx*vx + vy*vy + vz*vz;
}
double Particle::get_a2() {
  return ax*ax + ay*ay + az*az;
}
double Particle::get_j2() {
  return jx*jx + jy*jy + jz*jz;
}
double Particle::get_s2() {
  return sx*sx + sy*sy + sz*sz;
}
double Particle::get_c2() {
  return cx*cx + cy*cy + cz*cz;
}
double Particle::get_r() {
  return sqrt( get_r2() );
}
double Particle::get_v() {
  return sqrt( get_v2() );
}
double Particle::get_a() {
  return sqrt( get_a2() );
}
double Particle::get_j() {
  return sqrt( get_j2() );
}
double Particle::get_s() {
  return sqrt( get_s2() );
}
double Particle::get_c() {
  return sqrt( get_c2() );
}
  
void Particle::add_to_mass(double dm) {
  m += dm;
}
void Particle::add_to_radius(double dR) {
  R += dR;
}
void Particle::add_to_x(double dx) {
  x += dx;
}
void Particle::add_to_y(double dy) {
  y += dy;
}
void Particle::add_to_z(double dz) {
  z += dz;
}
void Particle::add_to_vx(double dvx) {
  vx += dvx;
}
void Particle::add_to_vy(double dvy) {
  vy += dvy;
}
void Particle::add_to_vz(double dvz) {
  vz += dvz;
}
void Particle::add_to_ax(double dax) {
  ax += dax;
}
void Particle::add_to_ay(double day) {
  ay += day;
}
void Particle::add_to_az(double daz) {
  az += daz;
}
void Particle::add_to_jx(double djx) {
  jx += djx;
}
void Particle::add_to_jy(double djy) {
  jy += djy;
}
void Particle::add_to_jz(double djz) {
  jz += djz;
}
void Particle::add_to_sx(double dsx) {
  sx += dsx;
}
void Particle::add_to_sy(double dsy) {
  sy += dsy;
}
void Particle::add_to_sz(double dsz) {
  sz += dsz;
}
void Particle::add_to_cx(double dcx) {
  cx += dcx;
}
void Particle::add_to_cy(double dcy) {
  cy += dcy;
}
void Particle::add_to_cz(double dcz) {
  cz += dcz;
}




