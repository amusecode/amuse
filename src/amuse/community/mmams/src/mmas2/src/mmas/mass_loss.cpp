#include "mmas.h"


/* This routine computes mass loss */
/* age is given in units of MS time */
real mmas::mass_loss() {
  real m_a = model_a->star_mass;
  real m_b = model_b->star_mass;
  real q = min(m_a, m_b)/max(m_a, m_b);
  real f_ml = 0;
  real age = max(model_a->star_age, model_b->star_age);
//   if (age < 0.75) {
//     f_ml = 6.93 * pow(q, -2.11) * pow(2*q/(1+q), 4.16);
//   } else {
  f_ml = 8.36 * pow(q, -2.58) * pow(2*q/(1+q), 4.28);
//   }		
  PRL(f_ml);
  return f_ml;
}
