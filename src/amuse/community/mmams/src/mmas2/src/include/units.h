#include "include/stdinc.h"

/* basic units */

#ifdef _UNITS_SI_

  const real uCM = 1.0e-2;   
  const real uGR = 1.0e-3;
  const real uSEC = 1;
  const real uERG = 1.0e-7;
  const real uKELVIN = 1;
  const real uAM = 3.3356e-10;
  const real uGS = 1.0e-4;
 
#else

  const real uCM = 1;
  const real uGR = 1;
  const real uSEC = 1;
  const real uERG = 1;
  const real uKELVIN = 1;
  const real uGS = 1;
  const real uAM = 1;

#endif


const real uMSUN       = 1.981e33 * uGR;
const real uRSUN       = 6.9599e10 * uCM;
const real uLSUN       = 3.826e33 * uERG/uSEC;
const real uC          = 2.997924800e10 * uCM/uSEC;
const real uG          = 6.67259e-8 * uCM*uCM*uCM/uGR/uSEC;
const real uE          = 4.8032068e-10 * uAM;
const real uH          = 6.6260755e-27 * uCM*uCM*uGR/uSEC;
const real uHBAR       = uH / (2.0*PI);
const real uM_E        = 9.1093897e-28 * uGR;
const real uM_P        = 1.6726231e-24 * uGR;
const real uM_N        = 1.674929e-24  * uGR;
const real uM_U        = 1.660540e-24  * uGR;
const real uALPHA      = 1.0/(uH*uC/(2*PI*uE*uE));    
const real uSIGMA_TH   = 6.6524616e-25 * uCM*uCM;
const real uK          = 1.380658e-16 * uERG/uKELVIN;
const real uN_A        = 6.0221367e+23;
const real uSIGMA_RAD  = 5.67051e-5 * uERG/uSEC/(uSEC*uSEC)/(uKELVIN*uKELVIN);
const real uA_RAD      = 4.0*uSIGMA_RAD/uC;
const real uEV         = 1.602192e-12 * uERG;
const real uKEV        = 1.0e3 * uEV;
const real uMEV        = 1.0e6 * uEV;

const real uHOUR       = 3600.0 * uSEC;
const real uDAY        = 24.0 * uHOUR;
const real uYR         = 365.242 * uDAY;
const real uLB         = 453.5924 * uGR;
const real uOZ         = uLB / 16.0;

const real uAU         = 1.49598e+13 * uCM;
const real uPC         = 3600.0*180.0/PI * uAU;
const real uLY         = uYR * uC;
const real uKPC        = 1.0e3 * uPC;
const real uMPC        = 1.0e6 * uPC;

const real uKG         = 1.0e-3 * uGR;
const real uM          = 100.0 * uCM;
const real uKM         = 1.0e3 * uM;
const real uIN         = 2.54 * uCM;
const real uFT         = 12 * uIN;
const real uMILE       = 1.609344 * uKM;
const real uMPH        = uMILE / uHOUR;
const real uKNOT       = 5.14444 * uCM / uSEC;
const real uKMH        = uKM / uHOUR;
const real uW          = 1.0e7 * uERG/uSEC;
const real uJY         = 1.0e-26 * uW / (uM*uM);
const real uMJY        = 1.0e-3 * uJY;
const real uDEG        = PI / 180.0;
const real uARCMIN     = uDEG / 60.0;
const real uARCSEC     = uARCMIN / 60.0;
