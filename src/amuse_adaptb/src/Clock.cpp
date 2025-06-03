#include "Clock.h"

#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////
Clock::Clock() {
  t = "0";
  t_end = "1";

  dt_print = "1e-2";
  t_print = "1e-2";

  dt_max = "1e-1";
  dt_factor = "1e-2";

  t_cpu = "0";
  dt = dt_max;
  t_progress = "0";
}
Clock::Clock(mpreal t_begin, mpreal t_sim, mpreal dt_print, mpreal dt_max, mpreal dt_factor) {
  t = t_begin;
  t_end = t_begin + t_sim;

  this->dt_print = dt_print;
  t_print = t_begin + dt_print;

  this->dt_max = dt_max;
  this->dt_factor = dt_factor;

  t_cpu = "0";
  dt = dt_max;
  t_progress = "0";
}
///////////////////////////////////////////////////////////////
// Initialization
///////////////////////////////////////////////////////////////
void Clock::initialize(mpreal t_begin, mpreal dt_print, mpreal dt_max, mpreal dt_factor) {
  t = t_begin;
//  t_end = t_begin + t_sim;

  this->dt_print = dt_print;
  t_print = t_begin + dt_print;

  this->dt_max = dt_max;
  this->dt_factor = dt_factor;

  t_cpu = "0";
  dt = dt_max;
  t_progress = "0";
}
///////////////////////////////////////////////////////////////
// Setters
///////////////////////////////////////////////////////////////
void Clock::set_dt(mpreal dt) {
  this->dt = dt;
}
void Clock::set_t_begin_and_end(mpreal t_begin, mpreal t_sim) {
  t = t_begin;
  t_end = t_sim;
  t_print = t_begin + dt_print;
}
void Clock::set_t(mpreal t) {
  this->t = t;
}
///////////////////////////////////////////////////////////////
// Get
///////////////////////////////////////////////////////////////
mpreal Clock::get_t() {
  return t;
}
mpreal Clock::get_dt() {
  return dt;
}
mpreal Clock::get_t_print() {
  return t_print0;
}
mpreal Clock::get_t_end() {
  return t_end;
}
///////////////////////////////////////////////////////////////
// runtime Clock
///////////////////////////////////////////////////////////////
void Clock::Start_timer() {
  gettimeofday(&Tvalue,&dummy);
  timerStarted = true;
}
void Clock::stop_timer() {
  if(timerStarted == false) t_cpu = "0";
  struct timeval Tvalue2;
  struct timezone dummy2;
  gettimeofday(&Tvalue2,&dummy2);
  mpreal StartTime =  ((mpreal) Tvalue.tv_sec +"1.e-6"*((mpreal) Tvalue.tv_usec));
  mpreal endTime =  ((mpreal) Tvalue2.tv_sec +"1.e-6"*((mpreal) Tvalue2.tv_usec));
  timerStarted = false;
  t_cpu = endTime-StartTime;      
} 
mpreal Clock::get_timer() {
  return t_cpu;
}
mpreal Clock::get_progress() {
  return t/t_end*"100.0";
}
mpreal Clock::read() {
  struct timeval Tvalue2;
  struct timezone dummy2;
  gettimeofday(&Tvalue2,&dummy2);
  mpreal StartTime =  ((mpreal) Tvalue.tv_sec  +"1.e-6"*((mpreal) Tvalue.tv_usec));
  mpreal endTime =  ((mpreal) Tvalue2.tv_sec + "1.e-6"*((mpreal) Tvalue2.tv_usec));
  mpreal t_elap = endTime-StartTime;  
  return t_elap;
}
///////////////////////////////////////////////////////////////
// nbody Clock
///////////////////////////////////////////////////////////////
int Clock::alarm() {
  if(t_print <= t_end+"0.5"*dt_print) return 0;
  else return 1;
}
int Clock::to_print() {
  if(t > t_print)
  {
    t_print0 = t_print;
    t_print += dt_print;
    return 1;
  }
  else return 0;
}
void Clock::abort() {
  t_print = t_end + dt_print;
}
void Clock::tick() {
  t += dt;
}
void Clock::calc_dt( mpreal dt_est ) {
  dt = dt_factor * dt_est;
  if(dt > dt_max) dt = dt_max;
}

