#include "Star.h"

#ifndef __FORCE_H
#define __FORCE_H

class Force
{
  public:
  mpreal softening_sq;

  // constructors
  Force();
  Force(mpreal soft);

  // Forces
  void gravity(Star* s1, Star* s2);
  void gravity_dt(Star* s1, Star* s2);
};

#endif


