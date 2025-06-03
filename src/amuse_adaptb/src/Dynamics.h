#include "Star.h"

#ifndef __DYNAMICS_H
#define __DYNAMICS_H

class Dynamics
{
  public:

  void update_r(Star* s, mpreal dt);
  void update_v(Star* s, mpreal dt);
};

#endif
