#include "Force.h"

//////////////////////////////////////////
// Constructors
//////////////////////////////////////////
Force::Force()
{
  softening_sq = "0";
}
Force::Force(mpreal soft)
{
  softening_sq = soft;
}
//////////////////////////////////////////
// Forces
//////////////////////////////////////////
void Force::gravity(Star* s1, Star* s2)
{
  mpreal dx = s2->x - s1->x;
  mpreal dy = s2->y - s1->y;
  mpreal dz = s2->z - s1->z;

  mpreal dr2 = dx*dx + dy*dy + dz*dz + softening_sq;
  mpreal dr1 = sqrt(dr2);  
  mpreal dr3 = dr2*dr1;

  mpreal a1x = s2->m/dr3*dx;
  mpreal a1y = s2->m/dr3*dy;
  mpreal a1z = s2->m/dr3*dz;

  mpreal a2x = s1->m/dr3*-dx;
  mpreal a2y = s1->m/dr3*-dy;
  mpreal a2z = s1->m/dr3*-dz;

  s1->add_ax( a1x );
  s1->add_ay( a1y );
  s1->add_az( a1z );
  s2->add_ax( a2x );
  s2->add_ay( a2y );
  s2->add_az( a2z );  
}
void Force::gravity_dt(Star* s1, Star* s2)
{
  mpreal dx = s2->x - s1->x;
  mpreal dy = s2->y - s1->y;
  mpreal dz = s2->z - s1->z;

  mpreal dr2 = dx*dx + dy*dy + dz*dz + softening_sq;
  mpreal dr1 = sqrt(dr2);  
  mpreal dr3 = dr2*dr1;

  mpreal a1x = s2->m/dr3*dx;
  mpreal a1y = s2->m/dr3*dy;
  mpreal a1z = s2->m/dr3*dz;

  mpreal a2x = s1->m/dr3*-dx;
  mpreal a2y = s1->m/dr3*-dy;
  mpreal a2z = s1->m/dr3*-dz;

  s1->add_ax( a1x );
  s1->add_ay( a1y );
  s1->add_az( a1z );
  s2->add_ax( a2x );
  s2->add_ay( a2y );
  s2->add_az( a2z );  

  mpreal da1_ij = a1x*a1x + a1y*a1y + a1z*a1z;
  mpreal da2_ij = a2x*a2x + a2y*a2y + a2z*a2z;
  mpreal dt1 = dr2 / da1_ij;
  mpreal dt2 = dr2 / da2_ij;

  if(dt1 < s1->dt) s1->dt = dt1;
  if(dt2 < s2->dt) s2->dt = dt2;
}


