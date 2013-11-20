
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"




double splint_xl_yl_D2yl(double t)
{
  double value;

  splint(xl,yl,D2yl,ZSIZE+1,t,&value);

  return value;
}
