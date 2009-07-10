/* SWIG based C interface to SSE, in case f2py does not work.
 *
 * Written by Evert Glebbeek 2009.
 *
 * Modelled on the SSE_muse_interface.f written by Steve McMillan
 */

#include <stdlib.h>
#include <stdio.h>
#include "config.h"

/* We need access to some FORTRAN 77 common blocks (to set parameters)
 * We're assuming that a C int == FORTRAN integer, C double == FORTRAN
 * real*8. We're also hoping that the C compiler aligns structs in the same
 * way that the FORTRAN compiler aligns COMMON blocks.
 * Caveat emptor.
 */

/* Set up some defines to hide the FORTRAN name mangling for global
 * symbols. This is mainly to make the code more readable.
 * Note that these are not expanded recursively by the preprocessor.
 */
#define value1 GLOBAL_FORTRAN_NAME(value1)
#define value3 GLOBAL_FORTRAN_NAME(value3)
#define value4 GLOBAL_FORTRAN_NAME(value4)
#define flags  GLOBAL_FORTRAN_NAME(flags)
#define points GLOBAL_FORTRAN_NAME(points)

extern struct {
   double neta, bwind, hewind, mxns;
} value1;

extern struct {
   double sigma;
   int bhflag;
} value4;

extern struct {
   int ceflag,tflag,ifflag,nsflag,wdflag;
} flags;

extern struct {
   double pts1,pts2,pts3;
} points;

extern struct {
 int idum;
} value3;

/* We need to store some global variables as well */
static double z, zpars[20];

/* Initialisation function */
int initialize( double z_in, double neta_in, double bwind_in,
                double hewind_in, double sigma_in, int ifflag_in,
                int wdflag_in, int bhflag_in, int nsflag_in, double mxns_in,
                double pts1_in, double pts2_in, double pts3_in )
{

   /* Set input parameters passed in from MUSE */
   z = z_in;
   value1.neta = neta_in;
   value1.bwind = bwind_in;
   value1.hewind = hewind_in;
   value4.sigma = sigma_in;
   flags.ifflag = ifflag_in;
   flags.wdflag = wdflag_in;
   value4.bhflag = bhflag_in;
   flags.nsflag = nsflag_in;
   value1.mxns = mxns_in;
   points.pts1 = pts1_in;
   points.pts2 = pts2_in;
   points.pts3 = pts3_in;
      
   /* Calculate fitting parameters as a function of z */
   GLOBAL_FORTRAN_NAME(zcnsts)(&z, zpars);
   if(value3.idum > 0) value3.idum = -value3.idum;

   return 0;
}

/* Evolve function.
 * Ideally we should use the "BOTH" specifier to use the pointers as input
 * *and* output, but this doesn't seem to work reliably with SWIG any more
 * than it does using f2py (same underlying problem?), so we do it this
 * way. Stupid and ugly, but it works.
 */
void evolve(/* Input */
            int kw, double mass, double mt, double r, double lum,
            double mc, double rc, double menv, double renv,
            double ospin, double epoch, double tm, double tphys,
            double tphysf,
            /* Output */
            int *kw1, double *mass1, double *mt1, double *r1, double *lum1,
            double *mc1, double *rc1, double *menv1, double *renv1,
            double *ospin1, double *epoch1, double *tm1, double *tphys1,
            double *tphysf1)
{
   double dtp;

   /* Copying everything is probably overkill, but not wrong. */
   *kw1 = kw;
   *mass1 = mass;
   *mt1 = mt;
   *r1 = r;
   *lum1 = lum;
   *mc1 = mc;
   *rc1 = rc;
   *menv1 = menv;
   *renv1 = renv;
   *ospin1 = ospin;
   *epoch1 = epoch;
   *tm1 = tm;
   *tphys1 = tphys;
   *tphysf1 = tphysf;

   /* We don't want any output. Seriously. */
   dtp = *tphys1+1;

   /* Evolve the star to the new time, return new stellar properties */
   GLOBAL_FORTRAN_NAME(evolv1) ( kw1, mass1, mt1, r1, lum1, mc1, rc1,
         menv1, renv1, ospin1, epoch1, tm1, tphys1, tphysf1, &dtp, &z, zpars);

   return;
}

/* Timestep query function */
double get_time_step(int kw, double mass, double age, double mt, double tm, double epoch)
{
   double tscls[20], lums[10], GB[10], tn, dtm, dtr, t;

   /* Call star fuction to get stellar parameters */
   GLOBAL_FORTRAN_NAME(star)(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);

   /* Call deltat function to get next timestep */
   t = age - epoch;
   GLOBAL_FORTRAN_NAME(deltat)(&kw, &t, &tm, &tn, tscls, &dtm, &dtr);

   /* Return the smaller of the two timescales */
   return (dtr<dtm) ? dtr : dtm;
}
