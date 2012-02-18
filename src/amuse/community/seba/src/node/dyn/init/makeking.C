
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct a King model.
////
//// Usage:  makeking [OPTIONS]
////
//// Options:
////              -b    specify Steve's rescaling parameter (< 1) [0]
////                    [models with b > 0 are just rescaled King models;
////                    models with b < 0 approach isothermal spheres
////                    as b --> -infinity]
////              -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -i    number the particles sequentially [don't number]
////              -n    specify number of particles [no default]
////                    if unspecified assumes an input snapshot with masses.
////              -o    echo value of random seed [don't echo]
////              -s    specify random seed [random from system clock]
////              -T    test options (print to cerr) [0]
////                      1: print King model, with unit central
////                      density and core radius, and exit;
////                      2: print King model, scaled to unit mass
////                      and virial radius, and exit;
////                      3: test realized velocity distribution
////              -u    leave final N-body system unscaled
////                    [scale to E=-1/4, M = 1, R = 1]
////              -w    specify King dimensionless depth [no default]
////
//// Written by Steve McMillan and Kimberly Engle.
////
//// Report bugs to starlab@sns.ias.edu.

// version 1:  June 1997    Kimberly Engle
//                          email: kim@galileo.physics.drexel.edu
//                          Physics Dept., Drexel University, Phila PA USA
//
//             Adapted from S. McMillan's FORTRAN version to the new
//             C++-based Starlab.

#include "dyn.h"

#ifdef TOOLBOX

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#define NMAX 50000

//----------------------------------------------------------------------

// Global stuff:
// ------------

real pi = M_PI;
real twopi = 2.0 * pi;
real four3pi = 4.0 * pi / 3.0;
real fourpi = 4.0 * pi; 
real zero = 0.0;
real one = 1.0;
real quarter = 0.25;
real four3 = 4.0 / 3;

// FORTRAN common blocks...

// parameter
#define  NM	2500

// profile
real rr[NM+1], d[NM+1], v2[NM+1], psi[NM+1], zm[NM+1];

// Convenient to keep these global (for now):

// Array indx is used in setpos and makeking
//		[indexes the mass distribution of the initial model],
//    "	   g  is used in gg and setvel
//		[contains the indefinite integral of y^2 exp(-y^2) dy].

#define NINDX	100
#define NG	1000
#define YMAX	4.0

int indx[NINDX+1];
real g[NG+1];

// Rescaling factors:

real beta = 0.0;
real beta_w0 = 0.0;
real scale_fac = 1.0;

//----------------------------------------------------------------------

#define FUNC(x) ((*func)(x))

real trapint(real (*func)(real), real a, real b, int n)

// Integrate func from a to b in n steps, using the trapezoid rule.

{
  static real s;
  int i, j;
  real x,  base, sum, dx;

  if (n == 1) {

      return (s = 0.5 * (b - a) * (FUNC(a) + FUNC(b)));

  } else {

      for (i = 1, j = 1; j < n - 1; j++) i <<= 1;
      base = i;   
      dx = (b - a) / base;
      x = a + 0.5 * dx;
      for (sum = 0.0, j = 1; j <= i; j++, x += dx) sum += FUNC(x);
      s = 0.5 * (s + (b - a)/base * sum);
      return s;
  }
}
#undef FUNC

#define MAX_STEPS 50

real trapzd(real (*func)(real), real a, real b, real eps)

// Trapezoid Rule for integration. Tricky code from Numerical Recipes...

{
  int i; 
  real sum, previous_sum;
  real diff, cond;

  previous_sum = -1.0e30;
  for (i = 1; i <= MAX_STEPS; i++) {

      sum = trapint(func, a, b, i);
      if (fabs(sum - previous_sum) < eps * fabs(previous_sum)) return sum;
      previous_sum = sum;
  }

  return 0.0;

}
#undef MAX_STEPS


real gaus2(real y)
{
  return y * y * exp(-y*y);
}


#define EPS  1.e-10

void initialize_global()
{
  real yo = 0;
  real dy = YMAX/NG;

  g[0] = 0;
  for (int i = 1; i <= NG; i++) {
    real y = yo + dy;
    g[i] = g[i-1] + trapzd(gaus2, yo, y, EPS);
    yo = y;
  }
}


void gg(real y, real& g2, real& g4)

// Array g contains the indefinite integral (from 0) of y^2 exp(-y^2) dy,
// up to a maximum value of YMAX, i.e.
//
//	g[i]  =  integral{0 to Y} y^2 exp(-y^2) dy,
//
// where Y  =  i*YMAX/NG (i = 0,...,NG).
//
// Function computes  g2  =  integral{0 to y} y^2 exp(-y^2) dy
//               and  g4  =  integral{0 to y} y^4 exp(-y^2) dy.

{
  real yindx = NG * y / YMAX;
  int iy = (int)yindx;

  if (iy >= NG)
      g2 = g[NG];
  else
      g2 = g[iy] + (yindx - iy) * (g[iy+1] - g[iy]);

  if (g4 > 0) {
      real y2 = y * y;
      g4 = 1.5 * g2 - 0.5 * y * y2 * exp(-y2);
  }
}

 
void get_dens_and_vel(real psi, real& dens, real& v2)

// Return density d and local velocity dispersion v2,
// given scaled potential psi (< 0).

{
  dens = 0;
  if (psi >= -beta_w0) return;

  // Distribution function is
  //
  //	f(r,v) = A exp(-psi) (exp(-v^2 / 2 sig^2) - exp(psi - beta*psi0)),
  //
  // where psi = phi/sig^2 and psi0 = -w0.
  //
  // Density is the integral of 4 pi v^2 f, velocity dispersion v2
  // is (integral of 4 pi v^4 f) / density.

  real g2, g4 = v2;

  real p_max = -psi - beta_w0;		// v_max^2 / (2 sig^2)

  real y_max = sqrt(p_max);
  real ep = exp(-psi);

  gg(y_max, g2, g4);

  dens = ep * g2 - scale_fac * y_max * p_max / 3;

  // Note that this expression for dens excludes an overall factor
  // of 8 pi sqrt(2) sig^3 A.  Scaling to the central density
  // is handled elsewhere (in rhs).

  if (v2 > 0 && dens > 0)
      v2 = 2 * (ep * g4 - 0.2 * scale_fac * y_max * p_max * p_max) / dens;

  // The unit of v2 is sig^2.

}

real dc_inverse;	// Convenient to keep this global, for now.
			// Set in poisson, used for scaling in rhs.

void rhs(real y[], real x, real ypr[])

//  Define RHS of ODE, for use by rk4.

{
  real d;

  ypr[0] = y[1];

  if (x <= 0) {

      d = 1; 

  } else {

      get_dens_and_vel(y[0]/x, d, zero);	// d is rho/rho_0
      						// zero suppresses vel
      d = d * dc_inverse;
  }

  ypr[1] = 9 * x * d;

}

//  Runge-Kutta-4 method 

void step(real y[], real& x, real dx, int N)

// Runge-Kutta-4 integrator.

{ 
    int i;

    // DEC C++ doesn't like these declarations with variable N:

    // real dydx[N], dy1[N], dy2[N], dy3[N], dy4[N],
    //               y1[N],  y2[N],  y3[N],  y4[N];

    real dydx[4], dy1[4], dy2[4], dy3[4], dy4[4],
                   y1[4],  y2[4],  y3[4],  y4[4];

    rhs(y, x, dydx);

    for (i = 0; i < N; i++) {
        dy1[i] = dx*dydx[i];
        y1[i] = y[i] + 0.5*dy1[i];
    }

    rhs(y1, x+0.5*dx, dydx);

    for (i = 0; i < N; i++) {
        dy2[i] = dx*dydx[i];
        y2[i] = y[i] + 0.5*dy2[i];
    }

    rhs(y2, x+0.5*dx, dydx); 

    for (i = 0; i < N; i++) {
        dy3[i] = dx*dydx[i];
        y3[i] = y[i] + dy3[i];
    }

    rhs(y3, x+dx, dydx);

    for (i = 0; i < N; i++) {
        dy4[i] = dx*dydx[i];
        y[i] += (dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])/6.0;
    }

    x += dx;
}

int time_to_stop(real y[], real x, real x_end, real dx)
{
    if (x < x_end - .5*dx)		// NOTE!!  Beware of rounding error!
        return 0;
    else 
        return 1;
}


void rk4(real& x, real x_next, real y[], int N, real dx)

// Integrate the y(x) from x to x_next, using an integration step of dx.
// On entry y is y at the initial x.  On exit, x is replaced by x_next
// and y is y(x).

{
    while (!time_to_stop(y, x, x_next, dx))
	step(y, x, dx, N);
}

#define RLIN 0.25
#define NLIN 105
#define RMAX 1e4
  
void poisson(real x[], int nmax, real w0, int& nprof, real& v20)

//       Self-contained 1-D (spherical) Poisson's equation solver.
//       Currently knows about normal and lowered King models.
//
//        Input:  nmax is the maximum number of points allowed
//                w0 is the dimensionless central potential
//                iout allows messages if nonzero

//        Output: x   is scaled radius (r/rc)
//                d   is scaled density (1 at center)
//                v2  is scaled velocity dispersion (1 at center)
//                psi is scaled potential (-W0 at center)
//                zm  is cumulative mass (scaling from x, d scalings)
//                nprof is the actual number of points generated
//                v20 is the central 3-D velocity dispersion (unit = sig^2)

{

  int i, iflag2; 
  real psi0, xn, xo, fac;
  real y[2];

  psi0 = - abs(w0);

  // Initialize at center of cluster.

  xn = 0;
  y[0] = 0;
  y[1] = psi0;
  x[0] = 0;
  psi[0] = psi0;
  v2[0] = 1;
  zm[0] = 0;

  // Establish density scaling factor.

  get_dens_and_vel(psi0, d[0], v2[0]);
  dc_inverse = 1./d[0];

  fac = pow(10, (log10(RMAX/RLIN) / (nmax-NLIN)));

// 	Poisson's equation is:
//
// 		(1/r^2) d/dr (r^2 dphi/dr)  =  4 pi G rho,
//
// 	where r is radius, phi is potential, and rho is density, given
// 	(for equal-mass stars) by
//
// 		rho	=  {integral (v < ve)} 4 pi v^2 f(v) dv,
//
// 	where ve is the escape velocity,
//
// 		ve^2	=  -2 phi.
//
// 	The (3-D) velocity distribution is
//
// 		f(v)	=  A (exp(-v^2 / 2 sig^2)
// 					 - exp(-ve^2 / 2 sig^2)),
//
// 	where sig^2 is a 1-D velocity dispersion (not quite the
// 	central velocity dispersion, except in the limit of infinite
// 	central potential).  In King's (1966) paper, he uses
// 	j^2 = 1 / (2 sig^2).
//
// 	Following King, we define the core radius rc by
//
// 		rc^2	=  9 sig^2 / (4 pi G rho0)
//
// 	and the dimensionless depth as
//
// 		W0	=  -phi0 / sig^2,
//
// 	where rho0 and phi0 are the central density and potential,
// 	respectively.
//
// 	We then scale as follows:
//
// 		x	=  r / rc
//
// 		d	=  rho / rho0
//
// 		psi	=  phi / sig^2,
//
// 	to obtain
//
// 		(x psi)''  =  9 x d,
//
// 	where ' = d/dx.
//
// 	We integrate this ODE from the cluster center (x = 0, d = 1,
// 	psi = -W0) to the tidal radius (d = 0) by defining
//
//		y(0)	=  (x psi)
//		y(1)	=  y(0)'
//
//	We cover the first RLIN core radii linearly with NLIN points;
//	the remaining coverage is logarithmic, out to RMAX core radii,
//	if necessary.  We stop when d <= 0.

  iflag2 = 0;

  for (i = 1; i <= nmax; i++) {

      xo = xn;
      if (i <= NLIN)
          xn = (RLIN * i) / NLIN;
      else
	  xn = fac * xo;

      real dx = 0.051*(xn-xo);

      rk4(xo, xn, y, 2, dx);

      //  N.B. Remember that y(1) is x*psi and xo is updated by step.

      xn = xo;

      x[i] = xn;
      psi[i] = y[0] / xn;

      v2[i] = 1;
      get_dens_and_vel(psi[i], d[i], v2[i]);

      if (d[i] < 0) {

 	// Density is negative, calculation is over.
 	// Interpolate to the tidal radius.

 	x[i] = x[i-1] + (x[i] - x[i-1]) / (0.1 - d[i]/d[i-1]);
 	d[i] = 0;
 	v2[i] = 0;

      }

      zm[i] = x[i] * y[1] - y[0];

      if (d[i] > 0) {

	  // Strange syntax because d = NaN (because of earlier error)
	  // will test FALSE in "if (d < 0)".

      } else {

          iflag2 = 1;
          break;

      }
  }

  if (iflag2 == 0) i = nmax;

  nprof = i;

  // Scale d and v2 to their central values.  Save v20 (unit = sig^2).

  v20 = v2[0];
  for (i = nprof; i >= 0; i--) {
      d[i] = d[i] / d[0];
      v2[i] = v2[i] / v2[0];
      zm[i] = (fourpi/9) * zm[i];
  }

}
#undef RLIN
#undef NLIN
#undef RMAX

void setpos(dyn * b, real& p)

// Obtain a random position for body b from the King profile
// and return the scaled potential at that location.

{
  //  Choose radius randomly from the mass distribution.

  real rno = randinter(0.0, 1.0);

  int i = (int)(NINDX * rno);
  int i1, iflag = 0;

  for (i1 = indx[i]; i1 <= indx[i+1]+1; i1++)
    if (zm[i1] > rno) {
      iflag = 1;
      break;
    }

  if (iflag == 0) err_exit("makeking: error in getpos");

  real rfac = (rno - zm[i1-1]) / (zm[i1] - zm[i1-1]);
  real r = rr[i1-1] + rfac * (rr[i1] - rr[i1-1]);

  p = psi[i1-1] + rfac * (psi[i1] - psi[i1-1]);

  //  Angular position random.

  real cth = 2 * randinter(0.0, 1.0) - 1;
  real sth = sqrt(1 - cth*cth);
  real ph = twopi * randinter(0.0, 1.0);
  real cph = cos(ph);
  real sph = sin(ph);

  b->set_pos(vec(r * sth * cph, r * sth * sph, r * cth));
}

void setvel(dyn * b, real p)

// Obtain a random velocity for body b from the King profile
// given its scaled potential p.

{
    static bool init_v33 = false;
    static real v33[NG+1];

    if (!init_v33) {
	for (int i = 0; i <= NG; i++)
	    v33[i] = scale_fac *  pow(((YMAX/NG) * i), 3) / 3;
	init_v33 = true;
    }

    // Array v33[] contains the second term in the integral for the density,
    // namely exp(beta*W0) * v_esc^3 / 3 (scaling v^2 by 2 sig^2, as usual).
    // As with g[], the array indices run from 0 to NG, spanning a range 0 to
    // YMAX, i.e. the cumulative distribution function for v is
    //
    //		exp(-p) * g[i] - v33[i],
    //
    // where y = i*YMAX/NG (i = 0,...,NG) and v = sqrt(2)*sig*y (sig = 1 here).

    //  Choose speed randomly from the distribution at this radius.

    real v = 0;

    if (p < -beta_w0) {

	real pfac = exp(-p);
	real ymax = sqrt(-p-beta_w0);

	// Will obtain v by bisection.  Determine maximum possible
	// range in the index i.

	int il = 0; 
	int iu = (int)((NG/YMAX) * sqrt(-p));	// Binning OK for W0 < 16,
						// *only* if beta >= 0.
	if (iu > NG) iu = NG;
      
	real rl = 0;
	real ru = pfac * g[iu] - v33[iu];

	real rno = randinter(0.0, 1.0) * ru;

	while (iu - il > 1) {
	    int  im = (il + iu) / 2;
	    real rm = pfac * g[im] - v33[im];
	    if (rm > rno) {
		iu = im;
		ru = rm;
	    } else {
		il = im;
		rl = rm;
	    }
	}

	// Maximum possible range of il here (for beta = 0) is
	//	0 to NG*sqrt(-p)/YMAX.
	// Maximum possible value of v (for beta = 0) is the local
	//      escape speed, sqrt(-2*p).

	v = (YMAX/NG) * sqrt(2.0) * (il + (rno - rl)/(ru - rl));
    }

    //  Direction is random.

    real cth = 2 * randinter(0.0,1.0) - 1;
    real sth = sqrt(1 - cth * cth);
    real ph = twopi * randinter(0.0, 1.0);
    real cph = cos(ph);
    real sph = sin(ph);

    b->set_vel(vec(v * sth * cph, v * sth * sph, v * cth));
}


local void dump_model_and_exit(int nprof, real mfac = 1, real rfac = 1)
{
    real rhofac = mfac/pow(rfac,3);
    real pfac = mfac/rfac;
    for (int i = 0; i <= nprof; i++)
	if (rr[i] > 0 && d[i] > 0)
	    cerr << i << "  "
		 << log10(rr[i]*rfac) << "  "
		 << log10(d[i]*rhofac) << "  "
		 << log10(-psi[i]*pfac) << "  "
		 << log10(v2[i]*pfac) << "  "
		 << rr[i]*rfac << "  "
		 << d[i]*rhofac << "  "
		 << -psi[i]*pfac << "  "
		 << v2[i]*pfac << "  "
		 << zm[i] << endl;		// zm is scaled to unit mass
    exit(0);
}


local void makeking(dyn * b, int n, real w0, bool n_flag, bool u_flag, int test)

// Create a King model, and optionally initialize an N-body system
// with total mass = n, core radius = 1.

{
    int i, iz, j, jcore, jhalf;
    real dz, z, rho0;
    real rhalf, zmcore;

    int nprof;
    real v20;

    if (w0 > 16) err_exit("makeking: must specify w0 < 16");

    initialize_global();

    // Compute the cluster density/velocity/potential profile
    
    poisson(rr, NM, w0, nprof, v20);

    if (test == 1)
	dump_model_and_exit(nprof);

    // Determine statistics and characteristic scales of the King model.

    rho0 = 1 / zm[nprof];		 // Central density for total mass = 1

    // Unit of velocity = sig, where rc^2 = 9 sig^2 / (4 pi G rho0)

    real sig = sqrt(four3pi * rho0 / 3); // This 3 was v20 in the f77 version...
					 // v20 is central vel. disp. / sig^2

    // Scale the zm array to unit total mass.

    for (i = 0; i <= nprof; i++)
	zm[i] = zm[i] / zm[nprof];

    // Index the mass distribution, and determine the core mass and
    // the half-mass radius.

    // By construction, rr[indx[j]] and rr[indx[j+1]] bracket the
    // radius containing a fraction j / NINDX of the total mass.

    indx[0] = 0;
    indx[NINDX] = nprof;

    dz = 1.0/NINDX;
    z = dz;
    iz = 1;
    for (j = 1; j <= nprof - 1; j++) {
	if (rr[j] < 1) jcore = j;
	if (zm[j] < 0.5) jhalf = j; 
	if (zm[j] > z) {
	    indx[iz] = j - 1;
	    z = z + dz;
	    iz = iz + 1;
	}
    }

    zmcore = zm[jcore] + (zm[jcore+1] - zm[jcore]) * (1 - rr[jcore])
		/ (rr[jcore+1] - rr[jcore]);

    rhalf = rr[jhalf] + (rr[jhalf+1] - rr[jhalf])
		* (0.5 - zm[jhalf])
		    / (zm[jhalf+1] - zm[jhalf]);

    // Compute the kinetic and potential energies, and determine the
    // virial radius and ratio.

    real kin = 0, pot =0;

    for (i = 1; i <= nprof; i++) {
	kin += (zm[i] - zm[i-1]) * (v2[i-1] + v2[i]);
	pot -= (zm[i] - zm[i-1]) * (zm[i] + zm[i-1]) / (rr[i-1] + rr[i]);
    }
    kin *= 0.25*sig*sig*v20;

    real rvirial = -0.5/pot;

    cerr << endl << "King model";
    cerr << "\n    w0 = " << w0 << "  beta = " << beta << "  nprof =" << nprof
	 <<        "  V20/sig2 = " << v20
	 <<        "  Mc/M = " << zmcore << endl
         <<   "    Rt/Rc = " << rr[nprof] << " (c = " << log10(rr[nprof])
         <<        ")  Rh/Rc = " << rhalf
         <<        "  Rvir/Rc = " << rvirial // << "  -T/U = " << -kin/pot
	 << endl
         <<   "    Rc/Rvir = " << 1/rvirial
         <<        "  Rh/Rvir = " << rhalf/rvirial
         <<        "  Rt/Rvir = " << rr[nprof]/rvirial
	 << "\n\n";

    if (test == 2) {

	// Scaling factors are for Mtotal = 1, Rvir = 1:

	dump_model_and_exit(nprof, rho0, 1/rvirial);
    }

    if (b == NULL || n < 1) return;

    // Initialize the N-body system.

    sprintf(tmp_string,
    "         King model, w0 = %.2f, Rt/Rc = %.3f, Rh/Rc = %.3f, Mc/M = %.3f",
	      w0, rr[nprof], rhalf, zmcore);
    b->log_comment(tmp_string);

    // Write essential model information to root dyn story.

    putrq(b->get_log_story(), "initial_mass", 1.0);

    // putrq(b->get_log_story(), "initial_rvirial", 0.25/kin); // assumes a lot!
							       // -- too much...

    putrq(b->get_log_story(), "initial_rtidal_over_rvirial",
	  rr[nprof] / (0.25/kin));

    // Assign positions and velocities. Note that it may actually
    // be preferable to do this in layers instead.

    for_all_daughters(dyn, b, bi) {

	if (test == 3) {

	    // Test: For a pure King model, getvel should generate a
	    //       distribution of velocities with maximum speed
	    //       sqrt(-2*p) and <v^2> = v2*v20, where p and v2
	    //       are the scaled potential and mean-square
	    //       velocity at any given radius.

	    real nsum = 10000;

	    for (int zone = 0; zone < 0.95*nprof; zone += nprof/15) {
		real v2sum = 0;
		real v2max = 0;
		for (int jj = 0; jj < nsum; jj++) {
		    setvel(bi, psi[zone]);
		    real vsq = bi->get_vel()*bi->get_vel();
		    v2sum += vsq;
		    v2max = Starlab::max(v2max, vsq);
		}

		cerr << "zone " << zone << "  r = " << rr[zone]
		     << "  v2max = " << v2max<< "  ?= "  << -2*psi[zone]
		     << "  v2mean = " << v2sum/nsum << "  ?= " << v2[zone]*v20
		     << endl;
	    }

	    exit(0);

	}

	if (n_flag)
	    bi->set_mass(1.0/n);

	real pot;
	setpos(bi, pot);
	setvel(bi, pot);

	// Unit of length = rc.
	// Unit of velocity = sig.

	bi->scale_vel(sig);
    }

    // System is in virial equilibrium in a consistent set of units
    // with G, core radius, and total mass = 1.

    // Convenient to have the "unscaled" system (-u on the command line)
    // be as close to standard units as possible, so rescale here to force
    // the virial radius to 1.  (Steve, 9/04)

    real xfac = 1/rvirial;
    real vfac = 1/sqrt(xfac);

    for_all_daughters(dyn, b, bi) {
	bi->set_pos(xfac*bi->get_pos());
	bi->set_vel(vfac*bi->get_vel());
    }

    // Transform to center-of-mass coordinates and optionally
    // scale to standard parameters.

    b->to_com();

    if (!u_flag && n > 1) {

	real kinetic, potential;

	// Note: scale_* operates on internal energies.

	get_top_level_energies(b, 0.0, potential, kinetic);
	scale_virial(b, -0.5, potential, kinetic);	// scales kinetic
	real energy = kinetic + potential;
	scale_energy(b, -0.25, energy);			// scales energy

	putrq(b->get_log_story(), "initial_total_energy", -0.25);
	putrq(b->get_log_story(), "initial_rvirial", 1.0);
	putrq(b->get_dyn_story(), "total_energy", -0.25);
    }
}

main(int argc, char ** argv)
{
    int  n;
    int  input_seed, actual_seed;

    bool c_flag = false;
    bool C_flag = false;
    bool i_flag = false;
    bool n_flag = false;
    bool o_flag = false;
    bool s_flag = false; 
    bool w_flag = false;
    bool u_flag = false;

    check_help();

    int test = 0;

    real w0;

    char  *comment;

    extern char *poptarg;
    int c;
    const char *param_string = "b:c:Cin:os:T:uw:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.15 $", _SRC_)) != -1)
        switch(c) {
            case 'b': beta = atof(poptarg);
                      break;
            case 'c': c_flag = true;
                      comment = poptarg;
                      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
                      break;
            case 'n': n_flag = true;
                      n = atoi(poptarg);
                      break;
	    case 'o': o_flag = true;
                      break;
            case 's': s_flag = true;
                      input_seed = atoi(poptarg);
                      break;
            case 'T': test = atoi(poptarg);
                      break;
	    case 'u': u_flag = true;
                      break;
            case 'w': w_flag = true;
                      w0 = atof(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
        }

    if (!w_flag) {
        cerr << "makeking: please specify the dimensionless depth";
        cerr << " with -w #\n";
        exit(1);
    }

    if (beta >= 1) {
        cerr << "makeking: beta < 1 required\n";
        exit(1);
    }

    beta_w0 = beta*w0;			// global variables!
    scale_fac = exp(beta_w0);

//    if (!n_flag && test != 1 && test != 2) {
//        cerr << "makeking: please specify the number # of";
//        cerr << " particles with -n #\n";
//        exit(1);
//    }

    dyn *b = NULL;

    if (!n_flag) {

	b = get_dyn();
	n = b->n_leaves();

	if (n < 1)
	    err_exit("makeking: n > 0 required");
 
	cerr << "makeking: read " << n << " masses from input snapshot with" 
	     << endl;
    }
    else {

	if (n < 1)
	    err_exit("makeking: n > 0 required");

	b = new dyn();
	b->set_root(b);
	dyn *by, *bo;

	bo = new dyn();
	if (i_flag)
	    bo->set_label(1);
	b->set_oldest_daughter(bo); 
	bo->set_parent(b);

	for (int i = 1; i < n; i++) {
	    by = new dyn();
	    if (i_flag)
		by->set_label(i+1);
	    by->set_parent(b);
	    bo->set_younger_sister(by);
	    by->set_elder_sister(bo);
	    bo = by;
	}

    }
    if (c_flag)
	b->log_comment(comment);		// add comment to story
    b->log_history(argc, argv);
    
    if (C_flag) b->set_col_output(true);

    if (s_flag == false) input_seed = 0;	// default
    actual_seed = srandinter(input_seed);
    
    if (o_flag)
	cerr << "makeking: random seed = " << actual_seed << endl;

    sprintf(tmp_string,
	    "       random number generator seed = %d",
	    actual_seed);
    b->log_comment(tmp_string);

    makeking(b, n, w0, n_flag, u_flag, test);

    put_dyn(b);
}

#endif

