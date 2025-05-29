#include <sys/time.h>

#include "mpreal.h"
using namespace mpfr;

#ifndef __CLOCK_H
#define __CLOCK_H

class Clock
{
  mpreal t, t_end;

  mpreal dt, dt_max, dt_factor;

  mpreal t_print, dt_print, t_print0;

  mpreal t_cpu, t_progress;
  struct timeval Tvalue;
  struct timezone dummy;
  bool timerStarted;

  public:

  // constructors
  Clock();
  Clock(mpreal t_begin, mpreal t_sim, mpreal dt_print, mpreal dt_max, mpreal dt_factor);

  // Setters
  void initialize(mpreal t_begin, mpreal dt_print, mpreal dt_max, mpreal dt_factor);
  void set_t_begin_and_end(mpreal t_begin, mpreal t_sim);
  void set_dt(mpreal dt);
  void set_t(mpreal t);

  // get
  mpreal get_t();
  mpreal get_dt();
  mpreal get_t_print();
  mpreal get_t_end();

  // runtime Clock
  void Start_timer();
  void stop_timer();
  mpreal get_timer();
  mpreal get_progress();
  mpreal read();

  // nbody Clock
  int alarm();
  int to_print();
  void abort();
  void tick();

  void calc_dt( mpreal dt_est );

};

#endif


