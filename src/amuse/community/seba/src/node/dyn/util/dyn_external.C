
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  dyn_external.C: functions related to external influences on the system.
//.............................................................................
//    version 1:  Jul 2001, Steve McMillan
//    version 2:  Sep 2001, Steve McMillan
//.............................................................................
//
// Member functions:
//
//	real dyn::get_external_scale_sq()
//	vec dyn::get_external_center()	(absolute)
//
// Global functions:
//
//	void set_friction_beta
//	void set_friction_mass
//	void set_friction_vel		(CM vel, absolute)
//	void set_friction_acc		(const acc for all member stars)
//	void get_external_acc		(single node)
//	real vcirc			(at radius r)
//	real get_external_pot		(root node)
//	real get_tidal_pot		(root node)
//	real get_plummer_pot		(root node)
//	real get_external_virial	(root node)
//	void print_external
//
//.........................................................................

#include "dyn.h"

// NOTES:  1. This file should be the *only* place where tidal and other
//            external fields are specified.
//
//	   2. Must add consistent functions add_xxx(), de_xxx_pot(), and
//	      xxx_virial() for each new external field xxx introduced.
//
//	   3. pot_external uses current positions (get_pos)
//	      get_external uses positions and velocities passed as
//	      arguments (may be pos or pred_pos, depending on the
//	      application)
//
//	   4. As of 6/03, all "centers" are defined in absolute terms.
//	      Centers of mass or density now include the pos or vel of
//	      the root node.  External field centers remain absolute.

//-------------------------------------------------------------------------
// Tidal field (quadrupole):
//-------------------------------------------------------------------------

#define USE_CORIOLIS		// comment out to suppres the Coriolis terms

local inline void add_tidal(dyn *b,
			    vec pos,
			    vec vel,
			    real& pot,
			    vec& acc,
			    vec& jerk,
			    bool pot_only)
{
    // Compute the tidal components of the acceleration, jerk, and pot
    // of top-level node b.  The actual position and velocity used are
    // pos and vel, assumed relative to b.  The node pointer b is used
    // as only a means of passing global dyn data.  We *assume* that the
    // root node for the system is correctly set.

    real a1 = b->get_alpha1();
    if (a1 == 0) return;

    real a3 = b->get_alpha3();

    vec dx = pos + b->get_root()->get_pos() - b->get_tidal_center();

    if (!pot_only) {

	vec da_tidal_centrifugal = -vec(a1*dx[0], 0.0, a3*dx[2]);
	vec da_coriolis = 0;
#ifdef USE_CORIOLIS
	da_coriolis = 2 * b->get_omega() * vec(vel[1], -vel[0], 0.0);
#endif

	// Must update acc BEFORE computing dj for velocity-dependent forces!

	acc += da_tidal_centrifugal + da_coriolis;

	vec dj_tidal_centrifugal = -vec(a1*vel[0], 0.0, a3*vel[2]);
	vec dj_coriolis = 0;
#ifdef USE_CORIOLIS

	// Note that the acc here *must* be the total acc, not just
	// the tidal component.

	dj_coriolis = 2 * b->get_omega() * vec(acc[1], -acc[0], 0.0);
#endif

	jerk += dj_tidal_centrifugal + dj_coriolis;
    }

    real x = dx[0];
    real z = dx[2];

    pot += 0.5*(a1*x*x + a3*z*z);
}

local inline real tidal_pot(dyn *b,
			    void (*pot_func)(dyn *, real) = NULL)
{
    // Determine the tidal component of the potential energy
    // of root node b.

    // Add tidal and centrifugal terms for top-level nodes only.
    // (No potential term for the Coriolis force, note.)

    real a1 = b->get_alpha1();
    if (a1 == 0) return 0;

    b->set_root(b);				// safety

    real a3 = b->get_alpha3();
    vec cen = b->get_tidal_center() - b->get_pos();

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	real x = bb->get_pos()[0] - cen[0];
	real z = bb->get_pos()[2] - cen[0];
	real dp = 0.5*(a1*x*x + a3*z*z);
	dpot += bb->get_mass() * dp;
	if (pot_func) pot_func(bb, dp);
    }

    return dpot;
}

//-------------------------------------------------------------------------
// Plummer field:
//-------------------------------------------------------------------------

local inline void add_plummer(dyn *b,
			      vec pos,
			      vec vel,
			      real& pot,
			      vec& acc,
			      vec& jerk,
			      bool pot_only)
{
    // Compute the Plummer-field components of the acceleration, jerk,
    // and pot of top-level node b.  The actual position and velocity
    // used are pos and vel, assumed relative to b.  The node pointer b
    // is used only as a means of passing global dyn data.  We *assume*
    // that the root node for the system is correctly set.

    real M = b->get_p_mass();
    if (M == 0) return;

    real a2 = b->get_p_scale_sq();

    vec dx = pos + b->get_root()->get_pos() - b->get_p_center();
    real r2 = square(dx) + a2;
    real r1 = sqrt(r2);

    if (!pot_only) {
	real r3i = 1/(r1*r2);
	acc -= M*dx*r3i;
	jerk += M*(3*dx*(dx*vel)/r2 - vel)*r3i;
    }

    pot -= M/r1;
}

// Add dynamical friction term.

static bool new_rhalf = true;
static int count_warning = 0;

void set_new_rhalf(bool s)	// default = true
{
    new_rhalf = s;
    count_warning = 0;
}

static real rhalf = 0;
static real density = 0;

local inline bool get_rhalf(dyn *b)
{
    if (!new_rhalf) return true;

    dyn *root = b->get_root();
    story *s = root->get_dyn_story();

    if (find_qmatch(s, "kira_rhalf")) {

	// Get radius and density from the root dynstory.

	rhalf   = getrq(s, "kira_rhalf");
	density = getrq(s, "kira_half_density");  // assume present if rhalf is

	if (rhalf <= 0 || density <= 0) {
	    if (++count_warning < 10)
		warning("set_rhalf: rhalf or density improperly set");
	    return false;
	}

    } else {

	// Recompute radius and density.  Code follows that in sys_stats.
	// Suppressed temporarily -- reference to lagrad stuff somehow
	// forces a libsstar dependency...

	rhalf = compute_lagrangian_radii(root, 2, false, 0);    // don't print

	if (rhalf > 0) {
	    putrq(s, "kira_rhalf", rhalf);
	    density = 1.5*root->get_mass() / (4*M_PI*pow(rhalf, 3));
	    putrq(s, "kira_half_density", density);
	} else {
	    if (++count_warning < 10)
	        warning("set_rhalf: computed rhalf <= 0");
	    return false;
	}
    }

    new_rhalf = false;
    return true;
}

local inline void add_plummer_friction(dyn *b,
				       vec pos,
				       vec vel,
				       real& pot,
				       vec& acc,
				       vec& jerk,
				       bool pot_only)
{
    // Code by F-C Lin (2004).

    real M = b->get_p_mass();
    if (M == 0) return;

    if (!get_rhalf(b)) return;

    // Set parameters to reduce/cut off the friction term.
    // Still need to deal more generally with the scaling of the
    // particle masses to the mass of the background field.
   
    real mass = b->get_mass();		// scale the mass of black hole??
    real speed = sqrt(square(vel));
    real a2 = b->get_p_scale_sq();
  
    vec dx = pos + b->get_root()->get_pos() - b->get_p_center();

    real r2 = square(dx) + a2;
    real r1 = sqrt(r2);
    real sigma2 = sqrt(M/r1);
                                       
    real p_density = (0.2387324*M/pow(a2,1.5))*pow(1+square(dx)/a2,-2.5);
					// 0.2387324 = 3/(4*pi)
    real p_core_dens = 0.2387324*M/pow(a2,1.5);
  
    real core_dens = density;
    real continue_factor=1.0;
    if (core_dens > p_core_dens)
      {
	continue_factor=1.0/(1+exp(10-sqrt(square(dx))/rhalf*10));
      }

    real beta = 1.6*continue_factor;	 // calibration from N-body experiments

    real ffac = beta*40*M_PI*p_density*mass;	// assume logLamda = 10
    real X = speed/sigma2;

    vec da = 0, dj = 0;

    if (X > 0.1) {

	real erfterm = (erf(X) - 2*X*exp(-X*X)/sqrt(M_PI));
	real ffact = ffac*erfterm*pow(speed,-3);
	da = -ffact*vel;
	dj = -ffac*(erfterm*(pow(speed,-3)*acc
			     - 3*pow(speed,-5)*(vel*acc)*vel)
		    + 4*pow(X,3)/sqrt(M_PI)*exp(-X*X)
		      		*(vel*acc)*pow(speed,-5)*vel);

    } else {

	real ffact = ffac* 4 / (3*sqrt(M_PI)*pow(sigma2, 3));
	da = -ffact*vel;
	dj = -ffact*acc;

    }

    // Apply cutoff/gradual reduction here, if desired.

    acc += da;
    jerk += dj;
}

local inline real plummer_pot(dyn *b,
			      void (*pot_func)(dyn *, real) = NULL)
{
    // Determine the Plummer-field component of the potential energy
    // of root node b.

    // Add potential terms to top-level nodes only.

    real M = b->get_p_mass();
    if (M == 0) return 0;

    b->set_root(b);				// safety

    real a2 = b->get_p_scale_sq();
    vec cen = b->get_p_center() - b->get_pos();

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	vec dx = bb->get_pos() - cen;
	real r2 = square(dx) + a2;
	real dp = -M/sqrt(r2);
	dpot += bb->get_mass() * dp;
	if (pot_func) pot_func(bb, dp);
    }

    return dpot;
}

local inline real plummer_virial(dyn *b)
{
    // Determine the Plummer-field component of the virial sum
    // of root node b.

    real M = b->get_p_mass();
    if (M == 0) return 0;

    b->set_root(b);				// safety
    int debug = 0;

    real a2 = b->get_p_scale_sq();

    // Don't make any assumptions about the locations of the
    // center of mass of the center of the Plummer field...

    vec com_pos, com_vel;
    compute_com(b, com_pos, com_vel);
    if (debug) PRL(com_pos);

    vec dR = com_pos - b->get_p_center();
    vec acc_com = dR * pow(square(dR)+a2, -1.5);
    if (debug) PRL(acc_com);

    // Note that we don't actually need the acc_com term, as it should
    // sum to zero in the loop below...

    vec dcom_pos = com_pos - b->get_pos();		// com quantities
    vec dcen_pos = b->get_p_center() - b->get_pos();	// include root node
    if (debug) PRL(dcom_pos);
    if (debug) PRL(dcen_pos);

    real vir = 0;
    for_all_daughters(dyn, b, bb) {
	vec dr = bb->get_pos() - dcom_pos;		// relative to com
	if (debug > 1) PRL(dr);
	dR = bb->get_pos() - dcen_pos;			// relative to p_center
	if (debug > 1) PRL(dR);
	vec acc_ext = dR * pow(square(dR)+a2, -1.5);
	real dvir = bb->get_mass()*dr*(acc_ext - acc_com);
	if (debug > 1) PRL(dvir);
	vir += dvir;
    }
    if (debug) PRL(vir);

    return -M*vir;
}

//=========================================================================
// Power-law field, M(r) = A r^x for r >> a.  Note that the implementation
// has changed as of 2/04.  We now use a pure power law in density for r > a,
// with constant density for r < a.  (We previously, used a as a softening
// parameter, in analogy to the Plummer model).  Thus, exponent = 0 no longer
// reduces to a Plummer field with mass = A.
//
// We preserve the above asymptotic form for compatibility with previous
// work (and with the literature on the density near the Galactic center).
// However, for ease of computation, the "primary" quantity is actually
// the density rho, which we now *define* to be
//
//	rho  =  (Ax / 4pi) a^(x-3)		r < a
//		(Ax / 4pi) r^(x-3)		r >= a,
//
// whence
//
//	M(r) =  (Ax a^x /3) (r/a)^3		r < a
//		 A a^x (x/3 - 1) + A r^x	r >= a.
//
// Dynamical friction is currently implemented only for the power-law case,
// and *not* for the Plummer field -- to be fixed.
//=========================================================================

//-------------------------------------------------------------------------
//
// Notes from Steve (10/01):
//
// For now, handle here the bits and pieces related to dynamical friction.
// The expression given by Binney & Tremaine is:
//
//	Afric = -4 pi log(Lambda) beta Mfric Vcm rho
//				[erf(X) - 2Xexp(-X^2)/sqrt(pi)] / |Vcm|^3
//
// where the non-obvious terms are defined below and Lambda ~ N(<r).
// We obtain N(<r) from M(<r) assuming a mean mass of 1 Msun and using
// known scalings.  The quantity "1 Msun" is defined if physical units
// are enabled; if they aren't, then it is not clear what scaling we
// should use, or indeed what the meaning of dyamical friction is...
// Beta is a "fudge factor," = 1 according to Binney and Tremaine.
// We introduce some flexibility by allowing beta to be specified on
// the command line, and letting log Lambda default to 1 in the case
// where no physical scale is known.
//
//-------------------------------------------------------------------------

static real beta = 0;				// tunable parameter;
void set_friction_beta(real b) {beta = b;}	// BT say beta = 1

static real Mfric = 0;				// cluster effective mass
void set_friction_mass(real m) {Mfric = m;}

static vec Vcm = 0;				// cluster CM velocity
void set_friction_vel(vec v) {Vcm = v;}		// (absolute)

local inline real pl_density(dyn *b, real r)	// background density
						
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a = b->get_pl_scale();
    real x = b->get_pl_exponent();

    if (r <= a)
	return (A*x/(4*M_PI))*pow(a,x-3);
    else
	return (A*x/(4*M_PI))*pow(r,x-3);
}

local inline real pl_mass(dyn *b, real r)	// mass interior to r, the
{						// distance from pl_center
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a = b->get_pl_scale();
    real x = b->get_pl_exponent();

    if (r <= a)
	return A*x*pow(a,x)*pow(r/a,3) / 3;
    else
	return A*pow(a,x)*(x/3-1) + A*pow(r,x);
}

#define LAMBDA_FAC	1

local inline real pl_logLambda(dyn *b, real r)
{
    real mass_unit = -1;
    if (b->get_starbase())
	mass_unit = b->get_starbase()->conv_m_dyn_to_star(1);

    // Only case where this is meaningful is the power-law field.

    real LogLambda;
    if (beta <= 0 || !b->get_pl())
	LogLambda = 0;

    if (mass_unit <= 0)				// no physical mass scale
	LogLambda = 1;

    else {

	// Use M(<r), assuming <m> = 1 Msun.

	LogLambda = 6.6; // Spinnato et al 2003

	// return log(LAMBDA_FAC*pl_mass(b, r)*mass_unit);
    }

    return LogLambda;
}

local inline real pl_vcirc2(dyn *b, real r)	// circular orbit speed:
						// recall vc^2 = r d(phi)/dr
{
    if (r > 0)
	return pl_mass(b,r)/r;
    else
	return 0;
}

static vec Afric = 0;				// frictional acceleration

void set_friction_acc(dyn *b,			// root node
		      real r)			// distance from pl_center
{
    // NOTE that Dynamical friction is computed using only the "pl"
    // density (even though the cluster mass may be defined using
    // multiple external fields).

    if (beta > 0) {

	real A = b->get_pl_coeff();
	if (A <= 0) return;

	real x = b->get_pl_exponent();
	real a = b->get_pl_scale();

	// Binney & Tremaine expression needs a 1-D velocity dispersion
	// sigma for the background stars.  Note that sigma2 here is
	// sqrt(2) * sigma.

	// Define sigma2 in terms of the circular orbit speed vc.

	real vc2 = pl_vcirc2(b, r), sigma2;

	if (x < 2 && r > a)
	    sigma2 = sqrt(vc2/(2-x));		// see McM & SPZ 2003
	else
	    sigma2 = sqrt(vc2/2);

	real V = abs(Vcm);
	real X = V/sigma2;			// scaled velocity; BT p. 425

	real coeff = 4*M_PI*beta*pl_logLambda(b, r);
	real ffac = coeff * Mfric * pl_density(b, r);

	if (X > 0.1)

	    ffac *= (erf(X) - 2*X*exp(-X*X)/sqrt(M_PI)) * pow(V, -3);

	else

	    // Expand for small X:

	    ffac *= 4 / (3*sqrt(M_PI)*pow(sigma2, 3));

	if (ffac < 0.2)				// arbitrary, but we expect
	    Afric = -ffac * Vcm;		// dt ~ 1, and want Adt < v
	else {
	    cerr << "  Suppressing dynamical friction at time "
		 << b->get_system_time() << endl;
	    Afric = 0;
	}

#if 1
	cerr << "  set_friction_acc: "; PRL(Afric);
	PRI(2); PRC(A); PRC(a); PRL(beta);
	PRI(2); PRC(coeff); PRC(Mfric); PRL(pl_density(b, r));
	PRI(2); PRC(r); PRC(sigma2); PRC(V); PRL(X);
	PRI(2); PRL((erf(X) - 2*X*exp(-X*X)/sqrt(M_PI)) * pow(V, -3));
	PRI(2); PRL(4 / (3*sqrt(M_PI)*pow(sigma2, 3)));
	PRI(2); PRC(ffac); PRL(Vcm);
#endif

    }
}

//-------------------------------------------------------------------------

static bool set = false;
static real ax, ax1;			// evaluate a^x etc. once per run

local inline void set_pl(dyn *b)
{
    if (!set) {
	real a = b->get_pl_scale();
	real x = b->get_pl_exponent();
	ax = pow(a,x);
	ax1 = ax/a;
	set = true;
    }
}

real rx;				// evaluate r^x once per call

local inline void set_rx(real x, real a, real r)
{
    rx = r;
    if (r > a && x != 1) rx = pow(r,x);
}

local inline real pl_dpot(real A, real x, real a, real r)
{
    if (r <= a)
	return A*x*ax1*pow(r/a,2)/6;
    else {
	if (x == 1)
	    return A    *(x/6 + (x/3-1)*(1-a/r) + log(r/a));
	else
	    return A*ax1*(x/6 + (x/3-1)*(1-a/r) + (rx/(r*ax1)-1)/(x-1));
    }
}

local inline void add_power_law(dyn *b,
				vec pos,
				vec vel,
				real& pot,
				vec& acc,
				vec& jerk,
				bool pot_only)
{
    // Compute the power-law-field components of the acceleration, jerk,
    // and pot of top-level node b.  The actual position and velocity
    // used are pos and vel, assumed relative to b.  The node pointer b
    // is used only as a means of passing global dyn data.  We *assume*
    // that the root node for the system is correctly set.

    real A = b->get_pl_coeff();
    if (A == 0) return;

    real a = b->get_pl_scale();
    real x = b->get_pl_exponent();

    if (!set) set_pl(b);

    vec dx = pos + b->get_root()->get_pos() - b->get_pl_center();
    real r2 = square(dx);
    real r = sqrt(r2);

    set_rx(x, a, r);
    pot += pl_dpot(A, x, a, r);

    if (!pot_only) {

	real m;					// mass interior to r

	if (r <= a)
	    m = A*x*ax*pow(r/a,3) / 3;
	else
	    m = A*ax*(x/3-1) + A*rx;		// <-- rx

	real r3i = 1/(r*r2);
	real mr3 = m*r3i;

	acc -= dx*mr3;

	jerk -= vel*mr3;
	if (r > a) {				// no extra piece for r < a
	    vec vr = dx*(dx*vel);
	    jerk += A*(3-x)*vr*r3i*rx/r2;	// <-- rx
	}
    }
}

local inline real power_law_pot(dyn *b,
				void (*pot_func)(dyn *, real) = NULL)
{
    // Determine the power-law-field component of the potential energy
    // of root node b.

    // Add potential terms to top-level nodes only.

    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    b->set_root(b);				// safety

    real a = b->get_pl_scale();
    real x = b->get_pl_exponent();
    vec cen = b->get_pl_center() - b->get_pos();

    if (!set) set_pl(b);

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	real r = abs(bb->get_pos()-cen);
	set_rx(x, a, r);
	real dp = pl_dpot(A, x, a, r);
	dpot += bb->get_mass()*dp;
	if (pot_func) pot_func(bb, dp);
    }

    return dpot;
}

local inline real power_law_virial(dyn *b)
{
    // Determine the power-law-field component of the virial sum
    // of root node b.

    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    b->set_root(b);				// safety
    int debug = 0;

    real a = b->get_pl_scale();
    real x = b->get_pl_exponent();

    // Don't make any assumptions about the locations of the
    // center of mass or the center of the power-law field...

    vec com_pos, com_vel;
    compute_com(b, com_pos, com_vel);
    if (debug) PRL(com_pos);

    vec dR = com_pos - b->get_pl_center();
    real r2 = square(dR);
    real r = sqrt(r2);
    
    vec acc_com = dR * pl_mass(b,r) / (r*r2);
    if (debug) PRL(acc_com);

    // Note that we don't actually need the acc_com term, as it should
    // sum to zero in the loop below...

    vec dcom_pos = com_pos - b->get_pos();		// com quantities
    vec dcen_pos = b->get_pl_center() - b->get_pos();	// include root node
    if (debug) PRL(dcom_pos);
    if (debug) PRL(dcen_pos);

    real vir = 0;
    for_all_daughters(dyn, b, bb) {
	vec dr = bb->get_pos() - dcom_pos;		// relative to com
	if (debug > 1) PRL(dr);
	dR = bb->get_pos() - dcen_pos;			// relative to pl_center
	if (debug > 1) PRL(dR);
	r2 = square(dR);
	r = sqrt(r2);
	vec acc_ext = dR * pl_mass(b,r) / (r*r2);
	real dvir = bb->get_mass()*dr*(acc_ext - acc_com);
	if (debug > 1) PRL(dvir);
	vir += dvir;
    }
    if (debug) PRL(vir);

    return -A*vir;
}

//-------------------------------------------------------------------------
// General "external" functions:
//-------------------------------------------------------------------------

// Member functions:

real dyn::get_external_scale_sq(int bit)	// default = -1
{
    // Just enumerate the possiblilties...

    // In the case of multiple external fields, just return the first
    // scale found if bit < 0 (default); if bit >= 0, return that
    // scale, if it is defined.

    if (!get_external_field()) return 0;

    if (get_tidal_field() && (bit < 0 || bit == 0)) return 0;
    if (get_plummer() && (bit < 0 || bit == 1))     return p_scale_sq;
    if (get_pl() && (bit < 0 || bit == 2))          return pow(pl_scale,2);

    return 0;
}

vec dyn::get_external_center(int bit)		// default = -1
{
    // Just enumerate the possiblilties...

    // In the case of multiple external fields, just return the first
    // scale found if bit < 0 (default); if bit >= 0, return that
    // scale, if it is defined.

    if (!get_external_field()) return 0;

    if (get_tidal_field() && (bit < 0 || bit == 0)) return tidal_center;
    if (get_plummer() && (bit < 0 || bit == 1))     return p_center;
    if (get_pl() && (bit < 0 || bit == 2))          return pl_center;

    return 0;
}

// Other functions:

void get_external_acc(dyn *b,
		      vec pos,
		      vec vel,
		      real& pot,
		      vec& acc,
		      vec& jerk,
		      bool pot_only)	// default = false
{
    // Compute the external components of the acceleration, jerk,
    // and pot of top-level node b, using the pos and vel provided,
    // assumed relative to b.  Add the external quantities to the acc,
    // jerk, and pot provided.  The node pointer b is used only as a
    // convenient means of passing static global dyn data.  We *assume*
    // that the root node for the system is correctly set.

    // *** Note that this version no longer returns just the tidal terms
    // *** in acc, jerk, and pot, but rather updates them.  (Steve, 4/05)

    unsigned int ext = b->get_external_field();

    if (ext) {

	// Loop through the known external fields.  Must do the
	// velocity-dependent tidal field last (and note that we will
	// have to be more careful when other velocity-dependent fields
	// are added, as *all* accs should be computed before jerks are
	// updated).

	if (GETBIT(ext, 1)) {
	    add_plummer(b, pos, vel, pot, acc, jerk, pot_only);
	    if (b->get_p_friction())
		add_plummer_friction(b, pos, vel, pot, acc, jerk, pot_only);
	}

	if (GETBIT(ext, 2))
	    add_power_law(b, pos, vel, pot, acc, jerk, pot_only);

	//if (GETBIT(ext, 3))
	//     add_other(b, pos, vel, pot, acc, jerk, pot_only);

	if (GETBIT(ext, 0))
	    add_tidal(b, pos, vel, pot, acc, jerk, pot_only);

	// Add dynamical friction term to non-escapers only:

	if (getiq(b->get_dyn_story(), "esc") == 0)		// too slow?
	    acc += Afric;
    }
}

real vcirc(dyn *b, vec r)
{
    // Return the circular velocity at position r.  Node b is used only
    // as a convenient means of passing static dyn class data.  We *assume*
    // that the root node for the system is correctly set.

    vec acc = 0, jerk = 0;
    real pot = 0;

    get_external_acc(b, r, vec(0), pot, acc, jerk);

    real vc2 = -r*acc;

    if (vc2 > 0)
	return sqrt(vc2);
    else
	return -sqrt(-vc2);	    // vcirc < 0 ==> no circular orbit exists
}

// Accessors:

real get_tidal_pot(dyn *b) {return tidal_pot(b);}
real get_plummer_pot(dyn *b) {return plummer_pot(b);}
real get_power_law_pot(dyn *b) {return power_law_pot(b);}

real get_external_pot(dyn *b,
		      void (*pot_func)(dyn *, real))	// default = NULL
{
    // Determine the external component of the potential of root
    // node b, using current positions.

    real pot = 0;
    unsigned int ext = b->get_external_field();

    b->set_root(b);				// safety

    if (ext) {

	real dpot = 0;

	// Loop through the known external fields.

	if (GETBIT(ext, 0)) dpot += tidal_pot(b, pot_func);
	if (GETBIT(ext, 1)) dpot += plummer_pot(b, pot_func);
	if (GETBIT(ext, 2)) dpot += power_law_pot(b, pot_func);
	// if (GETBIT(ext, 3)) dpot += other_pot(b, pot_func);

	if (pot_func) pot_func(b, dpot);	// used by kira to set_pot()

	pot += dpot;
    }

    return pot;
}

real get_external_virial(dyn *b)
{
    // Determine the external component of the potential of root
    // node b, using current positions.

    real vir = 0;
    unsigned int ext = b->get_external_field();

    b->set_root(b);				// safety

    if (ext) {

	// Loop through the known external, non-tidal fields.

	if (GETBIT(ext, 1)) vir += plummer_virial(b);
	if (GETBIT(ext, 2)) vir += power_law_virial(b);
	// if (GETBIT(ext, 3)) vir += other_virial(b);
    }

    return vir;
}

static bool sep = false;

local void print_ext(int n)
{
    if (sep) cerr << ",";

    switch(n) {
	case 0:		cerr << "TIDAL";
			break;
	case 1:		cerr << "PLUMMER";
			break;
	case 2:		cerr << "POWER-LAW";
			break;
	default:	cerr << "?";
			break;
    }

    sep = true;
}

void print_external(unsigned int ext,
		    bool shortbits)		// default = false
{
    if (!ext) return;

    if (shortbits)

	// Just print out ext as a binary number.

	printbits(ext);

    else {

	// Want to interpret the bits in ext.  Code follows printbits.C.

	int n = 31;
	while (n >= 0) {if (GETBIT(ext, n)) break; n--;}
	while (n >= 0) {if (GETBIT(ext, n)) print_ext(n); n--;}
	sep = false;
    }
}
