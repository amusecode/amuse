#include "copyright.h"
/*============================================================================*/
/*! \file utils.c
 *  \brief A variety of useful utility functions.
 *
 * PURPOSE: A variety of useful utility functions.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - ath_strdup()     - not supplied by fancy ANSI C, but ok in C89 
 * - ath_gcd()        - computes greatest common divisor by Euler's method
 * - ath_big_endian() - run-time detection of endianism of the host cpu
 * - ath_bswap()      - fast byte swapping routine
 * - ath_error()      - fatal error routine
 * - minmax1()        - fast Min/Max for a 1d array using registers
 * - minmax2()        - fast Min/Max for a 2d array using registers
 * - minmax3()        - fast Min/Max for a 3d array using registers
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn char *ath_strdup(const char *in)
 *  \brief This is really strdup(), but strdup is not available in 
 *   ANSI  (-pendantic or -ansi will leave it undefined in gcc)
 *   much like allocate.
 */

char *ath_strdup(const char *in)
{
  char *out = (char *)malloc((1+strlen(in))*sizeof(char));
  if(out == NULL) {
    ath_perr(-1,"ath_strdup: failed to alloc %d\n",(int)(1+strlen(in)));
    return NULL; /* malloc failed */
  }
  return strcpy(out,in);
}

/*----------------------------------------------------------------------------*/
/*! \fn int ath_gcd(int a, int b)
 *  \brief Calculate the Greatest Common Divisor by Euler's method
 */

int ath_gcd(int a, int b)
{
  int c;
  if(b>a) {c=a; a=b; b=c;} 
  while((c=a%b)) {a=b; b=c;}
  return b;
}

/*----------------------------------------------------------------------------*/
/*! \fn int ath_big_endian(void)
 *  \brief Return 1 if the machine is big endian (e.g. Sun, PowerPC)
 * return 0 if not (e.g. Intel)
 */

int ath_big_endian(void)
{
  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0); /* Returns 1 on a big endian machine */
}

/*----------------------------------------------------------------------------*/
/*! \fn void ath_bswap(void *vdat, int len, int cnt)
 *  \brief Swap bytes, code stolen from NEMO  
 */
 
void ath_bswap(void *vdat, int len, int cnt)
{
  char tmp, *dat = (char *) vdat;
  int k;
 
  if (len==1)
    return;
  else if (len==2)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  else if (len==4)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  else if (len==8)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  else {  /* the general SLOOOOOOOOOW case */
    for(k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void ath_error(char *fmt, ...)
 *  \brief Terminate execution and output error message
 *  Uses variable-length argument lists provided in <stdarg.h>
 */

void ath_error(char *fmt, ...)
{
  va_list ap;
   FILE *atherr = atherr_fp();

  fprintf(atherr,"### Fatal error: ");   /* prefix */
  va_start(ap, fmt);              /* ap starts with string 'fmt' */
  vfprintf(atherr, fmt, ap);      /* print out on atherr */
  fflush(atherr);                 /* flush it NOW */
  va_end(ap);                     /* end varargs */

#ifdef MPI_PARALLEL
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/*! \fn void minmax1(Real *data, int nx1, Real *dmino, Real *dmaxo)
 *  \brief Return the min and max of a 1D array using registers
 *  Works on data of type float, not Real.
 */

void minmax1(Real *data, int nx1, Real *dmino, Real *dmaxo)
{
  int i;
  register Real dmin, dmax;

  dmin = dmax = data[0];
  for (i=0; i<nx1; i++) {
    dmin = MIN(dmin,data[i]);
    dmax = MAX(dmax,data[i]);
  }
  *dmino = dmin;
  *dmaxo = dmax;
}

/*! \fn void minmax2(Real **data, int nx2, int nx1, Real *dmino, Real *dmaxo)
 *  \brief Return the min and max of a 2D array using registers
 *  Works on data of type float, not Real.
 */
void minmax2(Real **data, int nx2, int nx1, Real *dmino, Real *dmaxo)
{
  int i,j;
  register Real dmin, dmax;

  dmin = dmax = data[0][0];
  for (j=0; j<nx2; j++) {
    for (i=0; i<nx1; i++) {
      dmin = MIN(dmin,data[j][i]);
      dmax = MAX(dmax,data[j][i]);
    }
  }
  *dmino = dmin;
  *dmaxo = dmax;
}

/*! \fn void minmax3(Real ***data, int nx3, int nx2, int nx1, Real *dmino, 
 *                   Real *dmaxo)
 *  \brief Return the min and max of a 3D array using registers
 *  Works on data of type float, not Real.
 */
void minmax3(Real ***data, int nx3, int nx2, int nx1, Real *dmino, Real *dmaxo)
{
  int i,j,k;
  register Real dmin, dmax;

  dmin = dmax = data[0][0][0];
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
	dmin = MIN(dmin,data[k][j][i]);
	dmax = MAX(dmax,data[k][j][i]);
      }
    }
  }
  *dmino = dmin;
  *dmaxo = dmax;
}

/*----------------------------------------------------------------------------*/
/*! \fn  void do_nothing_bc(GridS *pG)
 *
 *  \brief DOES ABSOLUTELY NOTHING!  THUS, WHATEVER THE BOUNDARY ARE SET TO 
 *  INITIALLY, THEY REMAIN FOR ALL TIME.
 */
void do_nothing_bc(GridS *pG)
{
}

/*============================================================================
 * ERROR-ANALYSIS FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn Real compute_div_b(GridS *pG)
 *  \brief COMPUTE THE DIVERGENCE OF THE MAGNETIC FIELD USING FACE-CENTERED 
 *  FIELDS OVER THE ENTIRE ACTIVE GRID.  RETURNS THE MAXIMUM OF |DIV B|.
 */
Real compute_div_b(GridS *pG)
{
#ifdef MHD
  int i,j,k,is,ie,js,je,ks,ke;
  Real x1,x2,x3,divB,maxdivB=0.0;
  Real lsf=1.0,rsf=1.0,dx2=pG->dx2;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef CYLINDRICAL
        rsf = (x1+0.5*pG->dx1)/x1;  lsf = (x1-0.5*pG->dx1)/x1;
        dx2 = x1*pG->dx2;
#endif
        divB = (rsf*pG->B1i[k][j][i+1] - lsf*pG->B1i[k][j][i])/pG->dx1;
        if (je > js)
          divB += (pG->B2i[k][j+1][i] - pG->B2i[k][j][i])/dx2;
        if (ke > ks)
          divB += (pG->B3i[k+1][j][i] - pG->B3i[k][j][i])/pG->dx3;

        maxdivB = MAX(maxdivB,fabs(divB));
      }
    }
  }

  return maxdivB;

#else
  fprintf(stderr,"[compute_div_b]: This only works for MHD!\n");
  exit(EXIT_FAILURE);
  return 0.0;
#endif /* MHD */
}

/*----------------------------------------------------------------------------*/
/*! \fn void compute_l1_error(const char *problem, const MeshS *pM, 
 *                            const ConsS ***RootSoln, const int errortype) 
 *  \brief COMPUTE THE L1-ERRORS IN ALL VARIABLES AT THE CURRENT 
 *  (USUALLY THE FINAL)
 *  TIMESTEP USING THE INITIAL SOLUTION.  
 *
 *  THIS MEANS THAT THE SOLUTION MUST
 *  EITHER BE STATIC (STEADY-STATE) OR MUST HAVE COMPLETED A FULL PERIOD OF
 *  ROTATION, ETC.  FOR THE ERRORTYPE FLAG, 0 MEANS ABSOLUTE ERROR, AND
 *  1 MEANS AVERAGE ERROR PER GRID CELL.
 */
void compute_l1_error(const char *problem, const MeshS *pM, 
                      const ConsS ***RootSoln, const int errortype)
{
  DomainS *pD=&(pM->Domain[0][0]);
  GridS   *pG=pM->Domain[0][0].Grid;
  int i=0,j=0,k=0;
#if (NSCALARS > 0)
   int n;
#endif
  int is,ie,js,je,ks,ke;
  Real rms_error=0.0;
  Real x1,x2,x3,dVol,totVol;
  ConsS error,total_error;
  FILE *fp;
  char *fname, fnamestr[256];
  int Nx1,Nx2,Nx3;
#if defined MPI_PARALLEL
  double err[8+NSCALARS], tot_err[8+NSCALARS];
  int mpi_err;
#endif

  /* Clear out the total_error struct */
  memset(&total_error,0.0,sizeof(ConsS));
  if (pG == NULL) return;

/* compute L1 error in each variable, and rms total error */

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      memset(&error,0.0,sizeof(ConsS));
      for (i=is; i<=ie; i++) {
        dVol = 1.0;
        if (pG->dx1 > 0.0) dVol *= pG->dx1;
        if (pG->dx2 > 0.0) dVol *= pG->dx2;
        if (pG->dx3 > 0.0) dVol *= pG->dx3;
#ifdef CYLINDRICAL
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        dVol *= x1;
#endif

        /* Sum local L1 error for each grid cell I own */
        error.d   += dVol*fabs(pG->U[k][j][i].d   - RootSoln[k][j][i].d );
        error.M1  += dVol*fabs(pG->U[k][j][i].M1  - RootSoln[k][j][i].M1);
        error.M2  += dVol*fabs(pG->U[k][j][i].M2  - RootSoln[k][j][i].M2);
        error.M3  += dVol*fabs(pG->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
#ifdef MHD
        error.B1c += dVol*fabs(pG->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
        error.B2c += dVol*fabs(pG->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
        error.B3c += dVol*fabs(pG->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
        error.E   += dVol*fabs(pG->U[k][j][i].E   - RootSoln[k][j][i].E );
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          error.s[n] += dVol*fabs(pG->U[k][j][i].s[n] - RootSoln[k][j][i].s[n]);;
#endif
      }

      /* total_error is sum of local L1 error */
      total_error.d += error.d;
      total_error.M1 += error.M1;
      total_error.M2 += error.M2;
      total_error.M3 += error.M3;
#ifdef MHD
      total_error.B1c += error.B1c;
      total_error.B2c += error.B2c;
      total_error.B3c += error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
      total_error.E += error.E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) total_error.s[n] += error.s[n];
#endif
    }
  }

#ifdef MPI_PARALLEL
/* Now we have to use an All_Reduce to get the total error over all the MPI
 * grids.  Begin by copying the error into the err[] array */

  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
#ifdef MHD
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) err[8+n] = total_error.s[n];
#endif

/* Sum up the Computed Error */
  mpi_err = MPI_Reduce(err, tot_err, (8+NSCALARS),
                       MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
  if(mpi_err)
    ath_error("[compute_l1_error]: MPI_Reduce call returned error = %d\n",
              mpi_err);

/* If I'm the parent, copy the sum back to the total_error variable */
  if(pD->DomNumber == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
#ifdef MHD
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = err[8+n];
#endif

  }
  else return; /* The child grids do not do any of the following code */
#endif /* MPI_PARALLEL */

  /* Compute total number of grid cells */
  Nx1 = pD->Nx[0];
  Nx2 = pD->Nx[1];
  Nx3 = pD->Nx[2];

  totVol = 1.0;
  if (errortype == 1) {
    if (pD->MaxX[0] > pD->MinX[0]) totVol *= pD->MaxX[0] - pD->MinX[0];
    if (pD->MaxX[1] > pD->MinX[1]) totVol *= pD->MaxX[1] - pD->MinX[1];
    if (pD->MaxX[2] > pD->MinX[2]) totVol *= pD->MaxX[2] - pD->MinX[2];
#ifdef CYLINDRICAL
    totVol *= 0.5*(pD->MinX[0] + pD->MaxX[0]);
#endif
  }


/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
#ifdef MHD
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */

  rms_error = sqrt(rms_error)/totVol;

/* Print error to file "BLAH-errors.#.dat"  */
   sprintf(fnamestr,"%s-errors",problem);
   fname = ath_fname(NULL,fnamestr,NULL,NULL,1,0,NULL,"dat");

/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[compute_l1_error]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[compute_l1_error]: Unable to open file.\n");
      free(fname);
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
#ifdef MHD
    fprintf(fp,"  B1c  B2c  B3c");
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  S[ %d ]",n);
    }
#endif
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d /totVol),
	  (total_error.M1/totVol),
	  (total_error.M2/totVol),
	  (total_error.M3/totVol));

#ifndef ISOTHERMAL
  fprintf(fp,"  %e",total_error.E/totVol);
#endif /* ISOTHERMAL */

#ifdef MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/totVol),
	  (total_error.B2c/totVol),
	  (total_error.B3c/totVol));
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  %e",total_error.s[n]/totVol);
    }
#endif

  fprintf(fp,"\n");

  fclose(fp);
  free(fname);

  return;
}

/*============================================================================
 * ROOT-FINDING FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn int sign_change(Real (*func)(const Real,const Real), const Real a0, 
 *                      const Real b0, const Real x, Real *a, Real *b)
 *  \brief SEARCH FOR A SIGN CHANGE.  
 *
 *  THIS FUNCTION PARTITIONS THE INTERVAL (a0,b0) INTO
 *  2^k EQUALLY SPACED GRID POINTS, EVALUATES THE FUNCTION f AT THOSE POINTS,
 *  AND THEN SEARCHES FOR A SIGN CHANGE IN f BETWEEN ADJACENT GRID POINTS.  THE
 *  FIRST SUCH INTERVAL FOUND, (a,b), IS RETURNED.
 */
int sign_change(Real (*func)(const Real,const Real), const Real a0, 
                             const Real b0, const Real x, Real *a, Real *b) {
  const int kmax=20;
  int k, n, i;
  Real delta, fk, fkp1;

  for (k=1; k<=kmax; k++) {
    n = pow(2,k);
    delta = (b0-a0)/(n-1);
    *a = a0;
    fk = func(x,*a);
    for (i=1; i<n; i++) {
      *b = *a + delta;
      fkp1 = func(x,*b);
      if (fkp1*fk < 0)
        return 1;
      *a = *b;
      fk = fkp1;
    }
  }
//   ath_error("[sign_change]: No sign change was detected in (%f,%f) for x=%f!\n",a0,b0,x);
  return 0;
} 


/*----------------------------------------------------------------------------*/
/*! \fn int bisection(Real (*func)(const Real,const Real), const Real a0, 
 *                    const Real b0, const Real x, Real *root) 
 *  \brief THIS FUNCTION IMPLEMENTS THE BISECTION METHOD FOR ROOT FINDING.
 */
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0,
                           const Real x, Real *root) 
{
  const Real tol = 1.0E-10;
  const int maxiter = 400;
  Real a=a0, b=b0, c, fa, fb, fc;
  int i;

  fa = func(x,a);
  fb = func(x,b);
  if (fabs(fa) < tol) {
    *root = a;
    return 1;
  }
  if (fabs(fb) < tol) {
    *root = b;
    return 1;
  }
// printf("fa = %f, fb = %f\n", fa, fb);

  for (i = 0; i < maxiter; i++) {
    c = 0.5*(a+b);
// printf("x = %f, a = %f, b = %f, c = %f\n", x,a,b,c);
#ifdef MYDEBUG
    printf("c = %f\n", c);
#endif
    if (fabs((b-a)/c) < tol) {
#ifdef MYDEBUG
      printf("Bisection converged within tolerance of %f!\n", eps);
#endif
      *root = c;
      return 1;
    }
    fc = func(x,c);
    if (fa*fc < 0) {
      b = c;
      fb = fc;
    }
    else if (fc*fb < 0) {
      a = c;
      fa = fc;
    }
    else if (fc == 0) {
      *root = c;
      return 1;
    }
    else {
      ath_error("[bisection]:  There is no single root in (%f,%f) for x = %13.10f!!\n", a, b,x);
      *root = c;
      return 0;
    }
  }

  ath_error("[bisection]:  Bisection did not converge in %d iterations for x = %13.10f!!\n", maxiter,x);
  *root = c;
  return 0;
}


/*============================================================================
 * QUADRATURE FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn Real trapzd(Real (*func)(Real), const Real a, const Real b, const int n,
 *		    const Real s)
 *  \brief THIS ROUTINE COMPUTES THE nTH STAGE OF REFINEMENT OF AN EXTENDED 
 *  TRAPEZOIDAL RULE.  
 *
 * func IS INPUT AS A POINTER TO THE FUNCTION TO BE INTEGRATED BETWEEN 
 * LIMITS a AND b, ALSO INPUT. WHEN CALLED WITH n=1, THE ROUTINE RETURNS THE 
 * CRUDEST ESTIMATE OF \int_a^b f(R) R dR.  SUBSEQUENT CALLS WITH n=2,3,... 
 * (IN THAT SEQUENTIAL ORDER) WILL IMPROVE THE ACCURACY BY ADDING 2n-2 
 * ADDITIONAL INTERIOR POINTS. 
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER 
 */
Real trapzd(Real (*func)(Real), const Real a, const Real b, const int n, 
            const Real s)
{
  Real x,tnm,sum,dx;
  int it,j;

  if (n == 1) {
    return 0.5*(b-a)*(func(a)+func(b));
  } 
  else {
    for (it=1,j=1; j<n-1; j++) it <<= 1;  /* it = 2^(n-2) */
    tnm = it;
    dx = (b-a)/tnm;  /* THIS IS THE SPACING OF THE POINTS TO BE ADDED. */
    x = a + 0.5*dx;
    for (sum=0.0,j=1; j<=it; j++,x+=dx) sum += func(x);
    return 0.5*(s+(b-a)*sum/tnm);  /* THIS REPLACES s BY ITS REFINED VALUE. */
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn Real qsimp(Real (*func)(Real), const Real a, const Real b)
 *  \brief RETURNS THE INTEGRAL OF THE FUNCTION func FROM a TO b. 
 *
 * THE PARAMETER EPS 
 * CAN BE SET TO THE DESIRED FRACTIONAL ACCURACY AND JMAX SO THAT 2^(JMAX-1) 
 * IS THE MAXIMUM ALLOWED NUMBER OF STEPS.  INTEGRATION IS PERFORMED BY 
 * SIMPSON'S RULE.
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER 
 */

#define EPS 1.0e-8
#define JMAX 20

Real qsimp(Real (*func)(Real), const Real a, const Real b) 
{
  int j;
  Real s,st,ost,os;

  ost = os = -1.0e30;
  for (j=1; j<=JMAX; j++) {
    st = trapzd(func,a,b,j,ost);
    s = (4.0*st-ost)/3.0;  /* EQUIVALENT TO SIMPSON'S RULE */
    if (j > 5)  /* AVOID SPURIOUS EARLY CONVERGENCE. */
      if (fabs(s-os) < EPS*fabs(os) || (s == 0.0 && os == 0.0)) return s;
    os=s;
    ost=st;
  }

  ath_error("[qsimp]:  Too many steps!\n");
  return 0.0;
}


/*----------------------------------------------------------------------------*/
/* FUNCTION avg1d,avg2d,avg3d
 *
 * RETURNS THE INTEGRAL OF A USER-SUPPLIED FUNCTION func OVER THE ONE-, TWO-,
 * OR THREE-DIMENSIONAL GRID CELL (i,j,k).  INTEGRATION IS PERFORMED USING qsimp.
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER 
 */
static Real xsav,ysav,zsav,xmin,xmax,ymin,ymax,zmin,zmax;
static Real (*nrfunc)(Real,Real,Real);

/*! \fn Real avg1d(Real (*func)(Real, Real, Real), const GridS *pG, 
 *                 const int i, const int j, const int k)
 *  \brief RETURNS THE INTEGRAL OF A USER-SUPPLIED FUNCTION func OVER THE ONE-
 * DIMENSIONAL GRID CELL (i,j,k).  
 *
 * INTEGRATION IS PERFORMED USING qsimp.
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER 
 */
Real avg1d(Real (*func)(Real, Real, Real), const GridS *pG, 
            const int i, const int j, const int k)
{
  Real x1,x2,x3,dvol=pG->dx1;
  Real fx(Real x);

  nrfunc=func;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;

  ysav = x2;
  zsav = x3;
#ifdef CYLINDRICAL
  dvol *= x1;
#endif

  return qsimp(fx,xmin,xmax)/dvol;
}

/*! \fn Real avg2d(Real (*func)(Real, Real, Real), const GridS *pG, 
 *                 const int i, const int j, const int k)
 *  \brief RETURNS THE INTEGRAL OF A USER-SUPPLIED FUNCTION func OVER THE TWO-
 * DIMENSIONAL GRID CELL (i,j,k).  
 *
 * INTEGRATION IS PERFORMED USING qsimp.
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER 
 */
Real avg2d(Real (*func)(Real, Real, Real), const GridS *pG, 
            const int i, const int j, const int k)
{
  Real x1,x2,x3,dvol=pG->dx1*pG->dx2;
  Real fy(Real y);

  nrfunc=func;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  ymin = x2 - 0.5*pG->dx2;  ymax = x2 + 0.5*pG->dx2;

  zsav = x3;
#ifdef CYLINDRICAL
  dvol *= x1;
#endif

  return qsimp(fy,ymin,ymax)/dvol;
}

/*! \fn Real avg3d(Real (*func)(Real, Real, Real), const GridS *pG, 
 *                 const int i, const int j, const int k)
 *  \brief RETURNS THE INTEGRAL OF A USER-SUPPLIED FUNCTION func OVER THE
 * THREE-DIMENSIONAL GRID CELL (i,j,k).  
 *
 * INTEGRATION IS PERFORMED USING qsimp.
 * ADAPTED FROM NUMERICAL RECIPES BY AARON SKINNER
 */
Real avg3d(Real (*func)(Real, Real, Real), const GridS *pG, 
            const int i, const int j, const int k)
{
  Real x1,x2,x3,dvol=pG->dx1*pG->dx2*pG->dx3;
  Real fz(Real z);

  nrfunc=func;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  ymin = x2 - 0.5*pG->dx2;  ymax = x2 + 0.5*pG->dx2;
  zmin = x3 - 0.5*pG->dx3;  zmax = x3 + 0.5*pG->dx3;

#ifdef CYLINDRICAL
  dvol *= x1;
#endif

  return qsimp(fz,zmin,zmax)/dvol;
}

Real avgXZ(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k) {
  Real x1,x2,x3;

  Real fXZ(Real z);

  nrfunc=func;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  zmin = x3 - 0.5*pG->dx3;  zmax = x3 + 0.5*pG->dx3;

  ysav = x2;
  return qsimp(fXZ,zmin,zmax)/(x1*pG->dx1*pG->dx3);

}

Real fz(Real z)
{
  Real fy(Real y);

  zsav = z;
  return qsimp(fy,ymin,ymax);
}

Real fy(Real y)
{
  Real fx(Real x);

  ysav = y;
  return qsimp(fx,xmin,xmax);
}

Real fx(Real x)
{
#ifdef CYLINDRICAL
  return x*nrfunc(x,ysav,zsav);
#else
  return nrfunc(x,ysav,zsav);
#endif
}

Real fXZ(Real z) {
        Real fx(Real x);

        zsav = z;
        return qsimp(fx,xmin,xmax);

}

/*----------------------------------------------------------------------------*/
/* FUNCTION vecpot2b1i,vecpot2b2i,vecpot2b3i
 *
 * THESE FUNCTIONS COMPUTE MAGNETIC FIELD COMPONENTS FROM COMPONENTS OF A
 * SPECIFIED VECTOR POTENTIAL USING STOKES' THEOREM AND SIMPSON'S QUADRATURE.
 * NOTE:  THIS IS ONLY GUARANTEED TO WORK IF THE POTENTIAL IS OF CLASS C^1.
 * WRITTEN BY AARON SKINNER.
 */

static Real (*a1func)(Real,Real,Real);
static Real (*a2func)(Real,Real,Real);
static Real (*a3func)(Real,Real,Real);

/*! \fn Real vecpot2b1i(Real (*A2)(Real,Real,Real), Real (*A3)(Real,Real,Real),
 *                const GridS *pG, const int i, const int j, const int k)
 *  \brief Compute B-field components from a vector potential.
 *
 * THESE FUNCTIONS COMPUTE MAGNETIC FIELD COMPONENTS FROM COMPONENTS OF A
 * SPECIFIED VECTOR POTENTIAL USING STOKES' THEOREM AND SIMPSON'S QUADRATURE.
 * NOTE:  THIS IS ONLY GUARANTEED TO WORK IF THE POTENTIAL IS OF CLASS C^1.
 * WRITTEN BY AARON SKINNER.
 */
Real vecpot2b1i(Real (*A2)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,b1i=0.0,lsf=1.0,rsf=1.0,dx2=pG->dx2;
  Real f2(Real y);
  Real f3(Real z);

  a2func = A2;
  a3func = A3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  ymin = x2 - 0.5*pG->dx2;  ymax = x2 + 0.5*pG->dx2;
  zmin = x3 - 0.5*pG->dx3;  zmax = x3 + 0.5*pG->dx3;

  xsav = xmin;
#ifdef CYLINDRICAL
  lsf = xmin;  rsf = xmin;
  dx2 = xmin*pG->dx2;
#endif

  if (A2 != NULL) {
    if (ymin == ymax)
      b1i += rsf*A2(xmin,ymin,zmin) - lsf*A2(xmin,ymin,zmax);
    else {
      zsav = zmin;
      b1i += rsf*qsimp(f2,ymin,ymax);
      zsav = zmax;
      b1i -= lsf*qsimp(f2,ymin,ymax);
    }
  }
  if (A3 != NULL) {
    if (zmin == zmax)
      b1i += A3(xmin,ymax,zmin) - A3(xmin,ymin,zmin);
    else {
      ysav = ymax;
      b1i += qsimp(f3,zmin,zmax);
      ysav = ymin;
      b1i -= qsimp(f3,zmin,zmax);
    }
  }

  if (pG->dx2 > 0.0) b1i /= dx2;
  if (pG->dx3 > 0.0) b1i /= pG->dx3;

  return b1i;
}

/*! \fn Real vecpot2b2i(Real (*A1)(Real,Real,Real), Real (*A3)(Real,Real,Real),
 *                const GridS *pG, const int i, const int j, const int k)
 *  \brief Compute B-field components from a vector potential.
 *
 * THESE FUNCTIONS COMPUTE MAGNETIC FIELD COMPONENTS FROM COMPONENTS OF A
 * SPECIFIED VECTOR POTENTIAL USING STOKES' THEOREM AND SIMPSON'S QUADRATURE.
 * NOTE:  THIS IS ONLY GUARANTEED TO WORK IF THE POTENTIAL IS OF CLASS C^1.
 * WRITTEN BY AARON SKINNER.
 */
Real vecpot2b2i(Real (*A1)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,b2i=0.0;
  Real f1(Real x);
  Real f3(Real z);

  a1func = A1;
  a3func = A3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  ymin = x2 - 0.5*pG->dx2;  ymax = x2 + 0.5*pG->dx2;
  zmin = x3 - 0.5*pG->dx3;  zmax = x3 + 0.5*pG->dx3;

  ysav = ymin;

  if (A1 != NULL) {
    if (xmin == xmax)
      b2i += A1(xmin,ymin,zmax) - A1(xmin,ymin,zmin);
    else {
      zsav = zmax;
      b2i += qsimp(f1,xmin,xmax);
      zsav = zmin;
      b2i -= qsimp(f1,xmin,xmax);
    }
  }
  if (A3 != NULL) {
    if (zmin == zmax)
      b2i += A3(xmin,ymin,zmin) - A3(xmax,ymin,zmin);
    else {
      xsav = xmin;
      b2i += qsimp(f3,zmin,zmax);
      xsav = xmax;
      b2i -= qsimp(f3,zmin,zmax);
    }
  }

  if (pG->dx1 > 0.0) b2i /= pG->dx1;
  if (pG->dx3 > 0.0) b2i /= pG->dx3;

  return b2i;
}

/*! \fn Real vecpot2b3i(Real (*A1)(Real,Real,Real), Real (*A2)(Real,Real,Real),
 *                const GridS *pG, const int i, const int j, const int k)
 *  \brief Compute B-field components from a vector potential.
 *
 * THESE FUNCTIONS COMPUTE MAGNETIC FIELD COMPONENTS FROM COMPONENTS OF A
 * SPECIFIED VECTOR POTENTIAL USING STOKES' THEOREM AND SIMPSON'S QUADRATURE.
 * NOTE:  THIS IS ONLY GUARANTEED TO WORK IF THE POTENTIAL IS OF CLASS C^1.
 * WRITTEN BY AARON SKINNER.
 */
Real vecpot2b3i(Real (*A1)(Real,Real,Real), Real (*A2)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,b3i=0.0,lsf=1.0,rsf=1.0,dx2=pG->dx2;
  Real f1(Real x);
  Real f2(Real y);

  a1func = A1;
  a2func = A2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  xmin = x1 - 0.5*pG->dx1;  xmax = x1 + 0.5*pG->dx1;
  ymin = x2 - 0.5*pG->dx2;  ymax = x2 + 0.5*pG->dx2;
  zmin = x3 - 0.5*pG->dx3;  zmax = x3 + 0.5*pG->dx3;

  zsav = zmin;
#ifdef CYLINDRICAL
  rsf = xmax;  lsf = xmin;
  dx2 = x1*pG->dx2;
#endif

  if (A1 != NULL) {
    if (xmin == xmax)
      b3i += A1(xmin,ymin,zmin) - A1(xmin,ymax,zmin);
    else {
      ysav = ymin;
      b3i += qsimp(f1,xmin,xmax);
      ysav = ymax;
      b3i -= qsimp(f1,xmin,xmax);
    }
  }
  if (A2 != NULL) {
    if (ymin == ymax)
      b3i += rsf*A2(xmax,ymin,zmin) - lsf*A2(xmin,ymin,zmin);
    else {
      xsav = xmax;
      b3i += rsf*qsimp(f2,ymin,ymax);
      xsav = xmin;
      b3i -= lsf*qsimp(f2,ymin,ymax);
    }
  }

  if (pG->dx1 > 0.0) b3i /= pG->dx1;
  if (pG->dx2 > 0.0) b3i /= dx2;

  return b3i;
}

Real f1(Real x)
{
  return a1func(x,ysav,zsav);
}

Real f2(Real y)
{
  return a2func(xsav,y,zsav);
}

Real f3(Real z)
{
  return a3func(xsav,ysav,z);
}


#if defined(PARTICLES) || defined(CHEMISTRY)
/*----------------------------------------------------------------------------*/
/*! \fn void ludcmp(Real **a, int n, int *indx, Real *d)
 *  \brief LU decomposition from Numerical Recipes
 *
 * Using Crout's method with partial pivoting
 * a is the input matrix, and is returned with LU decomposition readily made,
 * n is the matrix size, indx records the history of row permutation,
 * whereas d =1(-1) for even(odd) number of permutations.
 */
void ludcmp(Real **a, int n, int *indx, Real *d)
{
  int i,imax,j,k;
  Real big,dum,sum,temp;
  Real *rowscale;  /* the implicit scaling of each row */

  rowscale = (Real*)calloc_1d_array(n, sizeof(Real));
  *d=1.0;  /* No row interchanges yet */

  for (i=0;i<n;i++)
  { /* Loop over rows to get the implicit scaling information */
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) ath_error("[LUdecomp]:Input matrix is singular!");
    rowscale[i]=1.0/big;  /* Save the scaling */
  }

  for (j=0;j<n;j++) { /* Loop over columns of Crout's method */
    /* Calculate the upper block */
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    /* Calculate the lower block (first step) */
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      /* search for the largest pivot element */
      if ( (dum=rowscale[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    /* row interchange */
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      rowscale[imax]=rowscale[j];
    }
    indx[j]=imax; /* record row interchange history */
    /* Calculate the lower block (second step) */
    if (a[j][j] == 0.0) a[j][j]=TINY_NUMBER;
    dum=1.0/(a[j][j]);
    for (i=j+1;i<n;i++) a[i][j] *= dum;
  }
  free(rowscale);
}

/*----------------------------------------------------------------------------*/
/*! \fn void lubksb(Real **a, int n, int *indx, Real b[])
 *  \brief Backward substitution (from numerical recipies)
 *
 * a is the input matrix done with LU decomposition, n is the matrix size
 * indx id the history of row permutation
 * b is the vector on the right (AX=b), and is returned with the solution
 */
void lubksb(Real **a, int n, int *indx, Real b[])
{
  int i,ii=-1,ip,j;
  Real sum;
  /* Solve L*y=b */
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  /* Solve U*x=y */
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void InverseMatrix(Real **a, int n, Real **b)
 *  \brief Inverse matrix solver
 *
 * a: input matrix; n: matrix size, b: return matrix
 * Note: the input matrix will be DESTROYED
 */
void InverseMatrix(Real **a, int n, Real **b)
{
  int i,j,*indx;
  Real *col,d;

  indx = (int*)calloc_1d_array(n, sizeof(int));
  col = (Real*)calloc_1d_array(n, sizeof(Real));

  ludcmp(a,n,indx,&d);

  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a, n, indx, col);
    for (i=0; i<n; i++)    b[i][j] = col[i];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c)
 *  \brief Matrix multiplication: a(m*n) * b(n*l) = c(m*l) */
void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c)
{
  int i, j, k;
  for (i=0; i<m; i++)
    for (j=0; j<l; j++)
    {
      c[i][j] = 0.0;
      for (k=0; k<n; k++) c[i][j] += a[i][k] * b[k][j];
    }
}
#endif /* PARTICLES or CHEMISTRY */
