#include "Star.h"

Star::Star() {
  x.assign(3,"0");
  v.assign(3,"0");
  a.assign(3,"0");
  a0.assign(3,"0");
#ifdef use_additional_acc
  a_step.assign(3,"0");;
#endif // use_additional_acc
}
Star::Star(mpreal m, vector<mpreal> x, vector<mpreal> v
#ifdef use_additional_acc
  , vector<mpreal> a_step
#endif // use_additional_acc
)
{
  this->m = m;
  this->x = x;
  this->v = v;
  #ifdef use_additional_acc
  this->a_step = a_step;
  a.assign(3, "0");
  #endif // use_additional_acc
}



