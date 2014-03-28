#include "Timestep.h"

Timestep::Timestep() {
  mode = 1;
  dt_max = 1e-1;
  dt_param = 1e-3;
}
Timestep::Timestep(double dt) {
  dt_param = 1.0;
  dt_max = dt;
}
Timestep::Timestep(int mode) {
  this->mode = mode;
  dt_max = 1e-1;
  dt_param = 1e-3;
}
Timestep::Timestep(int mode, double dt_max, double dt_param) {
  this->mode = mode;
  this->dt_max = dt_max;
  this->dt_param = dt_param;
}

void Timestep::set_mode(int mode) {
  this->mode = mode;
}
void Timestep::set_dt_max(double dt_max) {
  this->dt_max = dt_max;
}
void Timestep::set_dt_param(double dt_param) {
  this->dt_param = dt_param;
}
void Timestep::set_dt(double dt) {
  dt_param = 1;
  dt_max = dt;
}

int Timestep::get_mode() {
  return mode; 
}
double Timestep::get_dt_max() {
  return dt_max;
}
double Timestep::get_dt_param() {
  return dt_param;
}
double Timestep::get_dt() {
  return dt_max*dt_param;
}

double Timestep::get_dt(vector<Particle> &particle) {
  return get_constant_dt();
}

double Timestep::get_constant_dt() {
  return dt_param*dt_max;
}



