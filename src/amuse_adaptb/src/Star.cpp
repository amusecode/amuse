#include "Star.h"

//////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////
Star::Star()
{
  id = 0;
  radius = "0";
  m = "0";
  x = "0";
  y = "0";
  z = "0";
  vx = "0";
  vy = "0";
  vz = "0";
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag= "0";
  a2_mag = "0";
  dt = "0";
} 
Star::Star(mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ)
{
  id = 0;
  radius = "0";
  m = M;
  x = X;
  y = Y;
  z = Z;
  vx = VX;
  vy = VY;
  vz = VZ;
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag = "0";
  a2_mag = "0";
  dt = "0";
}
Star::Star(int id, mpreal M, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ)
{
  this->id = id;
  radius = "0";
  m = M;
  x = X;
  y = Y;
  z = Z;
  vx = VX;
  vy = VY;
  vz = VZ;
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag = "0";
  a2_mag = "0";
  dt = "0";
}
Star::Star(int id, mpreal M, mpreal R, mpreal X, mpreal Y, mpreal Z, mpreal VX, mpreal VY, mpreal VZ)
{
  this->id = id;
  radius = R;
  m = M;
  x = X;
  y = Y;
  z = Z;
  vx = VX;
  vy = VY;
  vz = VZ;
  ax = "0";
  ay = "0";
  az = "0";
  ax0 = "0";
  ay0 = "0";
  az0 = "0";
  r2_mag = "0";
  v2_mag = "0";
  a2_mag = "0";
  dt = "0";
}
//////////////////////////////////////////////
// Reset
//////////////////////////////////////////////
void Star::reset_a()
{
  ax0 = ax;
  ay0 = ay;
  az0 = az;
  ax = "0";
  ay = "0";
  az = "0";
}
void Star::reset_dt() {
  dt = "1e100";
}
//////////////////////////////////////////////
// Get
//////////////////////////////////////////////
mpreal Star::get_a2mag() {
  return a2_mag;
}
mpreal Star::get_v2mag() {
  return v2_mag;
}
mpreal Star::get_r2mag() {
  return r2_mag;
}
//////////////////////////////////////////////
// Add
//////////////////////////////////////////////
void Star::add_x(mpreal X)
{
  x += X;
}
void Star::add_y(mpreal Y)
{
  y += Y;
}
void Star::add_z(mpreal Z)
{
  z += Z;
}
void Star::add_vx(mpreal VX)
{
  vx += VX;
}
void Star::add_vy(mpreal VY)
{
  vy += VY;
}
void Star::add_vz(mpreal VZ)
{
  vz += VZ;
}
void Star::add_ax(mpreal AX)
{
  ax += AX;
}
void Star::add_ay(mpreal AY)
{
  ay += AY;
}
void Star::add_az(mpreal AZ)
{
  az += AZ;
}
//////////////////////////////////////////////
// Calc
//////////////////////////////////////////////
void Star::calc_r2mag()
{
  r2_mag = x*x+y*y+z*z;
}
void Star::calc_v2mag()
{
  v2_mag = vx*vx+vy*vy+vz*vz;
}
void Star::calc_a2mag()
{
  a2_mag = ax*ax+ay*ay+az*az;
}


