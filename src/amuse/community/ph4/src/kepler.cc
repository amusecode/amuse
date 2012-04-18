
// kepler: Kepler structure functions, imported from Starlab.
//         Initialize an orbit from mass, pos, and vel, or mass,
//         semi-major axis and eccentricity, and allow the user to
//         manipulate the resulting structure.
//
// Steve McMillan, Fall 1994, Spring 1999
// SPZ and JM, Feb 1998
// Steve McMillan, Feb 2011 - imported into AMUSE

#include "hdyn.h"

#define  ITERAC	   1.0e-12	// iteration accuracy for Kepler's equation
#define  MAXITER   100		// maximum number of iterations
#define  TRIG_TOL  1.e-12       // tolerance in sine and cosine functions
				//   _______________________________________
#define  SQRT_TRIG_TOL  1.e-6   // \/ tolerance in sine and cosine functions
#define	 ECC_TOL   1.0e-12	// problems with nearly linear orbits;
				// best to treat them as linear

//----------------------------------------------------------------------------
//
// Kepler conventions and use of variables:
// ---------------------------------------
//
//			E = energy
//			h = angular momentum
//			e = eccentricity
//			a = semi-major axis
//			P = period
//			n = mean motion
//			r = separation
//
// All definitions are standard, and Kepler's equation is used in the usual
// fashion, with the following conventions:
//
//		E < 0			E = 0			E > 0
//
// h > 0	e < 1			e = 1			e > 1
//		a < _INFINITY_		a = _INFINITY_		a = _INFINITY_
//		P < _INFINITY_		P = _INFINITY_		P = _INFINITY_
//					n --> n'
//
// h = 0	e = 1			e = 1			e = 1
//		a < _INFINITY_		a = _INFINITY_		a = _INFINITY_
//		P < _INFINITY_		P = _INFINITY_		P = _INFINITY_
//					n --> n''
//		true_an = ecc_an	true_an = r^(3/2)	true_an = ecc_an
//					mean_an = true_an
//
//----------------------------------------------------------------------------

#ifndef TOOLBOX

static real cos_to_sin(real c) // special treatment near |cos| = 1
{
    if (fabs(c) >= 1) return 0;
    
    real eps = 1 - fabs(c);
    if (eps > TRIG_TOL) return sqrt(1 - c*c);
    
    return sqrt(2*eps);
}

//----------------------------------------------------------------------------

// Simple mechanism for handling trigonometric errors in the Kepler package:

static int kepler_tolerance_level = 0;
static bool print_trig_warning = false;

// Options:

//	0: quit on any trig error
//	1: quit on large trig error
//	2: force to periastron or apastron on trig error

void set_kepler_tolerance(int i)
{
    if (i < 0) i = 0;
    if (i > 3) i = 3;
    kepler_tolerance_level = i;
}

void set_kepler_print_trig_warning(bool p)
{
    print_trig_warning = p;
}

int kepler_tolerance()
{
    return kepler_tolerance_level;
}

static int check_trig_limit(kepler* k, real &c, const char *s)
{
    int forced = 0;

    if (fabs(c) > 1 + TRIG_TOL) {
	if (kepler_tolerance_level == 0) {
	    cout.precision(HIGH_PRECISION);
	    cout << s << ": c = " << c << endl << flush;
	    k->print_all(cout);
	    err_exit("cosine out of range");
	} else if (kepler_tolerance_level == 1) {
	    if (c > 1) {
		if (c > 1 + TRIG_TOL) {
		    cout.precision(HIGH_PRECISION);
		    cout << s << ": c = " << c << endl << flush;
		    k->print_all(cout);
		    err_exit("cosine > 1");
		}
		c = 1;
	    } else if (c < -1 - TRIG_TOL) {
		if (c < -1 - TRIG_TOL) {
		    cout.precision(HIGH_PRECISION);
		    cout << s << ": c = " << c << endl << flush;
		    k->print_all(cout);
		    err_exit("cosine < -1");
		}
		c = -1;
	    }
	} else {
	    if (print_trig_warning) {
		int p = cout.precision(HIGH_PRECISION);
		cout << "warning: " << s << ": c = " << c << endl << flush;
		cout.precision(p);
	    }
	    c = fmax(-1.0, fmin(1.0, c));
	}

	forced = 1;
    }

    return forced;
}

#define MAX_WARNINGS 10
static int nwarnings = 0;

static void err_or_warn(const char *s)
{
    if (kepler_tolerance_level > 0) {

	// Unlimited warnings are a problem for automated scripts...
	// Place a conservative limit here (Steve, 8/02).

	warning(s);

	if (kepler_tolerance_level == 3 && ++nwarnings > MAX_WARNINGS) {
	    cout << "err_or_warn: too many warnings: "; PRL(MAX_WARNINGS);
	    kepler_tolerance_level = 0;
	}

    } else
	err_exit(s);
}

//----------------------------------------------------------------------------

// Local "Kepler" functions:
// ------------------------

//----------------------------------------------------------------------------
//  keplerseq  --  solve Kepler's equation by iteration, for elliptic orbits
//		   litt: H. Goldstein (1980), Classical Mechanics, eq. (3-76);
//  			 S.W. McCuskey (1963), Introduction to Celestial
//					       Mechanics, Sec. 3-7.
//		   method: true_an follows from ecc_an, which is determined
//			   by inverting Kepler's equation:
//			   mean_an = ecc_an - ecc * sin(ecc_an)
//
//                 Note that ecc = 1 is OK (linear bound orbit).
//
//----------------------------------------------------------------------------
//

static int keplerseq(real mean_an,	// mean anomaly
		     real ecc,		// eccentricity
		     real &true_an,	// true anomaly (returns in (-pi,pi])
		     real &ecc_an)	// eccentric anomaly
{
    if (ecc < 0) {			// inconsistent
	cout << "keplerseq: eccentricity e = " << ecc << " < 0\n";
	exit (1);
    }

    if (ecc > 1) {			// inconsistent
	cout << "keplerseq: eccentricity e = " << ecc << " > 1\n";
	exit (1);
    }

    mean_an = sym_angle(mean_an);	// -pi < mean_an < pi
    int iter = 0;


    if (mean_an == 0)			// special case

        ecc_an = 0;

    else {

	real  delta_ecc_an;	       	// iterative increment in ecc_an
	real  function;			// function = 0  solves Kepler's eq.
	real  derivative;	       	// d function / d ecc_an

	ecc_an = mean_an;		// first guess for ecc_an
	delta_ecc_an = 1;		// just to start the while loop

	int counter = 2;

	while (counter > 0) {

	    if (fabs(delta_ecc_an) < ITERAC) counter--;

	    if (++iter > MAXITER) {
	        PRC(delta_ecc_an);PRL(ITERAC);
		err_or_warn("keplerseq: convergence too slow");
		break;
	    }

	    function = -mean_an + ecc_an - ecc * sin(ecc_an);
			       // function = 0 solves Kepler's equation
	    derivative = 1 - ecc * cos(ecc_an);
			       // d(function) / d(ecc_an)
	    delta_ecc_an = -function / derivative;
			       // use Newton's method to find roots

	    if (delta_ecc_an > 1)
	      delta_ecc_an = 1;
	    else if (delta_ecc_an < -1)
	      delta_ecc_an = -1;	// avoid large jumps

	    ecc_an += delta_ecc_an;
	}
    }

    // Note convention that true_an = ecc_an if ecc = 1.

    if (ecc < 1)
      true_an = 2 * atan(sqrt((1 + ecc)/(1 - ecc)) * tan(ecc_an/2));
    else
      true_an = ecc_an;

    return iter;
}

//----------------------------------------------------------------------------
//  keplershypeq  --  solves Kepler's equation by iteration: hyperbolic orbits
//		      litt: S.W. McCuskey (1963), Introduction to Celestial
//					          Mechanics, Sec. 3-10.
//		      method: true_an follows from ecc_an, which is determined
//			      by inverting the hyperbolic analog of
//			      Kepler's equation:
//			      mean_an = -ecc_an + ecc// sinh(ecc_an)
//
//                 Note that ecc = 1 is OK (linear unbound orbit).
//
//----------------------------------------------------------------------------

static int  keplershypeq(real mean_an,	// mean anomaly
			 real ecc,	// eccentricity
			 real &true_an,	// true anomaly
			 real &ecc_an)	// eccentric anomaly (hyperbolic analog)
{
    real  delta_ecc_an;			// iterative increment in ecc_an
    real  function;			// function = 0  solves Kepler's eq.
    real  derivative;			// d function / d ecc_an
    int  i;

    if (ecc < 1) {			// inconsistent
	cout << "keplershypeq: eccentricity e = " << ecc << " < 1\n";
	exit (1);
    }

    if (mean_an == 0) 			// special case

	ecc_an = 0;

    else {

	    ecc_an = asinh(mean_an / ecc);  // first guess for ecc_an

	    i = 0;
	    delta_ecc_an = 1;		      // to start the while loop

	    int counter = 2;

	    while (counter > 0) {

	        if (fabs(delta_ecc_an) < ITERAC) counter--;

		if (++i > MAXITER) {
		    cout << "keplershypeq:  mean_an = " << mean_an
			 << "  ecc = " << ecc << "  MAXITER = " << MAXITER
			 << endl;
		    err_or_warn("keplershypeq: convergence too slow");
		    break;
		}

		function = -mean_an - ecc_an + ecc * sinh(ecc_an);
				    // function = 0 solves Kepler's equation
		derivative = -1 + ecc * cosh(ecc_an);
				    // d(function) / d(ecc_an)
		delta_ecc_an = -function / derivative;
				    // use Newton's method to find roots

		// PRI(4); PRC(function); PRC(derivative); PRL(delta_ecc_an);

		if  (fabs(derivative) < 1) { // avoid large jumps
		    if (delta_ecc_an > 1)
		      delta_ecc_an = 1;
		    else if (delta_ecc_an < -1)
		      delta_ecc_an = -1;
		}
		ecc_an += delta_ecc_an;
	    }
	}

    // Note convention that true_an = ecc_an if ecc = 1.

    if (ecc > 1)
        true_an = 2 * atan(sqrt((ecc + 1)/(ecc - 1)) * tanh(ecc_an/2));
    else
        true_an = ecc_an;

    return i;
}

//----------------------------------------------------------------------------
//  keplerspareq  --  solves Kepler's equation for parabolic orbits
//		      litt: S.W. McCuskey (1963), Introduction to Celestial
//					          Mechanics, Sec. 3-9.
//		      method: true_an (f) is determined by inverting the
//                            parabolic analog of Kepler's equation:
//			      mean_an = tan(f/2)// (1 + tan(f/2)*tan(f/2) / 3)
//----------------------------------------------------------------------------

static real keplerspareq(real mean_an,  // mean anomaly
			real &true_an)  // true anomaly

// This should be used ONLY for non-linear zero-energy orbits.

{
    if (mean_an == 0) {

	true_an = 0;
	return 0;

    } else {

	// This analytic solution is probably expensive...

	real s = atan(2/(3*mean_an));
	real w = atan(pow(tan(0.5*fabs(s)),1./3));
	if (s < 0) w = -w;
	real t = 2/tan(2*w);
	true_an = 2*atan(t);
	
	return t;
    }
}

//----------------------------------------------------------------------------

// Basic Kepler class member functions:
// -----------------------------------

kepler::kepler()
{
    time = total_mass = energy = angular_momentum = 0;
    rel_pos = rel_vel = 0;
    normal_unit_vector = 0;
    longitudinal_unit_vector = 0;
    transverse_unit_vector = 0;
    semi_major_axis = eccentricity = 0;
    true_anomaly = mean_anomaly = 0;
    time_of_periastron_passage = 0;
    pred_time = pred_separation = pred_true_anomaly = 0;
    pred_rel_pos = pred_rel_vel = 0;
    mean_motion = period = periastron = separation = pred_mean_anomaly = 0;
    circular_binary_limit = 0;
    tidal_potential = 0;
    delta_angular_momentum = 0;
}

kepler* hdyn_to_kepler(hdyn * com,    	    // com = center of mass
		       real t)		    // default = 0
{
    kepler *k = new kepler;
    hdyn *d1, *d2;

    if (!(d1 = com->get_oldest_daughter()))
	err_exit("hdyn_to_kepler: no oldest daughter present");
    if (!(d2 = d1->get_younger_sister()))
	err_exit("hdyn_to_kepler: no second daughter present");

    k->set_time(t);
    k->set_total_mass(d1->get_mass() + d2->get_mass());
    k->set_rel_pos(d2->get_pos() - d1->get_pos());
    k->set_rel_vel(d2->get_vel() - d1->get_vel());
    k->initialize_from_pos_and_vel();

    return k;
}

void new_kepler(hdyn * com,	      	// com is the center-of-mass hdyn
		real t)			// default = 0
{
    kepler * k = hdyn_to_kepler(com, t);
    com->set_kepler(k);
}

void kepler::print_all(ostream & s)
{
    int p = s.precision(HIGH_PRECISION);

    s << "  time             = " << time << endl;
    s.precision(12);
    s << "  total_mass       = " << total_mass    << endl;
    s << "  rel_pos          = " << rel_pos    << endl;
    s << "  rel_vel          = " << rel_vel << endl;
    s << "  separation       = " << separation << endl;
    s << "  |rel_vel|        = " << abs(rel_vel) << endl;
    real t = (energy < 0 ? sym_angle(true_anomaly) : true_anomaly);
    s << "  true_anomaly     = " << t << " (" << t*180./M_PI << ")" << endl;
    real m = (energy < 0 ? sym_angle(mean_anomaly) : mean_anomaly);
    s << "  mean_anomaly     = " << m << " (" << m*180./M_PI << ")" << endl;
    s << "  semi_major_axis  = " << semi_major_axis << endl;
    s << "  eccentricity     = " << eccentricity << endl;
    s << "  circular_binary_limit = " << circular_binary_limit << endl;
    s << "  periastron       = " << periastron << endl;
    s << "  apastron         = " << (energy >= 0 ? _INFINITY_ : 
				       semi_major_axis * (1 + eccentricity))
                                    << endl;
    s << "  energy           = " << energy << endl;
    s << "  energy_check     = " << 0.5*square(rel_vel) - total_mass/separation
				    << endl;
    s << "  angular_momentum = " << angular_momentum << endl;
    s << "  period           = " << period << endl;
    s << "  mean_motion      = " << mean_motion << endl;
    s << "  normal_unit_vector       = " << normal_unit_vector << endl;
    s << "  longitudinal_unit_vector = "
				    << longitudinal_unit_vector << endl;
    s << "  transverse_unit_vector   = " << transverse_unit_vector << endl;
    s << "  time_of_periastron_passage = "
				    << time_of_periastron_passage << endl;

    s.precision(p);
}

void kepler::print_dyn(ostream & s)
{
    s << "  time             = " << time << endl;
    s << "  rel_pos          = " << rel_pos    << endl;
    s << "  rel_vel          = " << rel_vel << endl;
    s << "  separation       = " << separation << endl;
    s << "  true_anomaly     = " << (energy < 0 ? sym_angle(true_anomaly) :
					true_anomaly) << endl;
    s << "  mean_anomaly     = " << (energy < 0 ? sym_angle(mean_anomaly) :
					mean_anomaly) << endl;
    s << "  time_of_periastron_passage = "
				    << time_of_periastron_passage << endl;
}

void kepler::print_elements(ostream & s)
{
    int p = s.precision(8);

    s << "  time             = " << time << endl;
    s << "  separation       = " << separation << endl;
    s << "  semi_major_axis  = " << semi_major_axis    << endl;
    s << "  eccentricity     = " << eccentricity << endl;
    s << "  circular_binary_limit = " << circular_binary_limit << endl;
    s << "  periastron       = " << periastron << endl;
    s << "  apastron         = " << (energy >= 0 ? _INFINITY_ : 
				       semi_major_axis * (1 + eccentricity))
                                    << endl;
    s << "  energy           = " << energy << endl;
    s << "  angular_momentum = " << angular_momentum << endl;
    s << "  period           = " << period << endl;

    s.precision(p);
}

//----------------------------------------------------------------------------

// Manipulation of orbit:
// ---------------------

//----------------------------------------------------------------------------
//  kepler::pred_true_to_mean_anomaly  --  sign conventions:
//                                    elliptic orbit:
//                                      mean anomaly will return with a value
//					in (-pi, pi), independent of the
//                                      initial sign of the true anomaly.
//                                    parabolic or hyperbolic orbit:
//                                      both true and mean anomaly will return
//                                      with the same sign as the initial true
//                                      anomaly.
//                                    linear orbit:
//                                      true_anomaly is not defined, so assume
//                                      that true_anomaly = eccentric_anomaly
//                                      on entry.
//
//----------------------------------------------------------------------------

void kepler::pred_true_to_mean_anomaly()
{
    real  ecc_anomaly;

    if (energy < 0) {			    // ellipse or bound linear orbit

	if (eccentricity < 1) {
	    ecc_anomaly = acos((eccentricity + cos(pred_true_anomaly)) /
			       (1 + eccentricity*cos(pred_true_anomaly)));
	    if (pred_true_anomaly < 0) ecc_anomaly = -ecc_anomaly;
	    
	} else

	    ecc_anomaly = pred_true_anomaly;	// convention in linear case

	pred_mean_anomaly = ecc_anomaly - eccentricity * sin(ecc_anomaly);

    } else if (energy > 0) {		// hyperbola or unbound linear orbit

	if (eccentricity > 1 + ECC_TOL) {	// avoid rounding error in the
						// determination of ecc_anomaly

	    // Check that the true anomaly is legal.

	    real  maxtruean = M_PI - acos(1 / eccentricity);

	    if (pred_true_anomaly < -maxtruean) {
		cout << "pred_true_to_mean_anomaly: hyp. true anom. = "
		  << pred_true_anomaly << " < minimum = " << -maxtruean <<"\n";
		exit (1);
	    }

	    if (pred_true_anomaly > maxtruean) {
		cout << "pred_true_to_mean_anomaly: hyp. true anom. = "
		  << pred_true_anomaly << " > maximum = " << maxtruean << "\n";
		exit (1);
	    }

	    // Outside chance of an error if the denominator in the ecc_anomaly
	    // expression is zero (this may be a redundant check here).

	    real ecc_arg,  denom = 1 + eccentricity * cos(pred_true_anomaly);
	    if (fabs(denom) < ECC_TOL)
		ecc_arg = eccentricity - 1/eccentricity;
	    else
		ecc_arg = (eccentricity + cos(pred_true_anomaly)) / denom;

	    ecc_anomaly = acosh(ecc_arg);
	    if (pred_true_anomaly < 0) ecc_anomaly = -ecc_anomaly;

	} else

  	    ecc_anomaly = pred_true_anomaly;	// convention in linear case
	    
	pred_mean_anomaly = -ecc_anomaly + eccentricity * sinh(ecc_anomaly);

    } else {			    		// parabola or linear orbit

	// Note the non-standard definitions of both the true and the mean
        // anomaly in this case.

	if (angular_momentum > 0) {		// parabola

	    real tan_half_true_an = tan(0.5*pred_true_anomaly);
	    pred_mean_anomaly = tan_half_true_an
	                         *(1 + tan_half_true_an*tan_half_true_an/3);

	} else					// linear (note convention)

	    pred_mean_anomaly = pred_true_anomaly;

    }
}

//----------------------------------------------------------------------------
//  kepler::mean_anomaly_to_periastron_passage  --  sign conventions:
//                                    elliptic orbit:
//                                      the previous function
//                                      pred_true_to_mean_anomaly() will
//                                      guarantee that the mean anomaly is
//                                      already positive.
//                                      The time_of_periastron_passage is the
//                                      time of the previous periastron
//                                      passage.
//                                    hyperbolic orbit:
//                                      the previous function
//                                      pred_true_to_mean_anomaly() allows
//                                      either sign for the mean anomaly.
//                                      The time_of_periastron_passage is the
//                                      time of the one and only periastron
//                                      passage, which may be either in the
//                                      past or in the future with respect to
//                                      the present time `time'.
//----------------------------------------------------------------------------

void kepler::mean_anomaly_to_periastron_passage()
{

    // Note that the mean anomaly is defined so that this equation is
    // correct in ALL cases.

    time_of_periastron_passage = time - mean_anomaly / mean_motion;
}

void kepler::update_time_of_periastron()
    {

    if (energy >= 0) return;

    real  dt = (time - time_of_periastron_passage) / period;

    int  dt_int = (int) fabs(dt);	// conversion is machine-dependent!
    if (dt < 0.0) dt_int = -dt_int - 1;

    time_of_periastron_passage += dt_int * period;
    }

void kepler::set_real_from_pred()
{
    time = pred_time;
    rel_pos = pred_rel_pos;
    rel_vel = pred_rel_vel;
    separation = pred_separation;
    true_anomaly = pred_true_anomaly;
    mean_anomaly = pred_mean_anomaly;
    update_time_of_periastron();
}

void kepler::to_pred_rel_pos_vel(real cos_true_an, real sin_true_an)
{

    vec r_unit = cos_true_an * longitudinal_unit_vector +
                    sin_true_an * transverse_unit_vector;
    pred_rel_pos = pred_separation * r_unit;

    real rel_vel_squared = 2 * (energy + total_mass / pred_separation);
    real v_t = angular_momentum / pred_separation;
    real v_r = sqrt(fmax(0.0, rel_vel_squared - v_t * v_t));
               // it is assumed here that negative values for the argument
               // arise from round-off errors, and can be replaced by zero
    if (sin_true_an < 0) v_r = -v_r;

    pred_rel_vel = v_r * r_unit
	            + v_t * (-sin_true_an * longitudinal_unit_vector
			     + cos_true_an * transverse_unit_vector);
}

void kepler::to_pred_rel_pos_vel_linear(real true_an) // special case
{
    pred_rel_pos = -pred_separation * longitudinal_unit_vector;

    if (pred_separation > 0) {
        real v_r = sqrt(fmax(0.0, 2 * (energy + total_mass / pred_separation)));
	if (true_an < 0) v_r = -v_r;
	pred_rel_vel = -v_r * longitudinal_unit_vector;
    } else {
        warning("to_pred_rel_pos_vel_linear: r <= 0");
        pred_rel_vel = _INFINITY_*transverse_unit_vector;
    }
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// New "fast" static data and accessors:

static bool kepler_fast_flag = false;
static real kepler_fast_tol = 1.0e-6;

void set_kepler_fast_flag(bool flag)		// default = true
						{kepler_fast_flag = flag;}
void clear_kepler_fast_flag()			{kepler_fast_flag = false;}
bool get_kepler_fast_flag()			{return kepler_fast_flag;}

void set_kepler_fast_tol(real tol)		// default = 1.e-6
						{kepler_fast_tol = tol;}
void reset_kepler_fast_tol()			{kepler_fast_tol = 1.e-6;}
real get_kepler_fast_tol()			{return kepler_fast_tol;}

static inline void fast_keplerseq(real M, real e, real tol,
				 real &s, real &c)
{
    // Newton-Raphson solver, with "Danby" initial guess.
    // See test_solver.C for timing code.

    M = sym_angle(M);			// -pi < M < pi

    // Make an initial guess for E (Danby would have a coefficient of 0.85).

    real dE = 0.94*e;
    if (M < 0) dE = -dE;

    real E = M + dE;

    while (1) {

	s = sin(E);
//	c = cos(E);

	// May be marginally faster to use the trig. identity:

	c = sqrt(1 - s*s);
	if (fabs(E) > M_PI/2) c = -c;

	real func = E - e*s - M;
	if (fabs(func) <= tol) return;

	E -= func/(1 - e*c);	// NR iteration
    }
}

// New kepler member functions (Steve, 6/01):

void kepler::fast_to_pred_rel_pos_vel(real r_cos_true_an,	// r cos f
				      real r_sin_true_an)	// r sin f
{
    // Relative position vector:

    pred_rel_pos = r_cos_true_an * longitudinal_unit_vector +
                    r_sin_true_an * transverse_unit_vector;

    // Now compute the relative velocity vector.

    real ri = 1 / pred_separation;
    vec r_unit = ri * pred_rel_pos;

    real rel_vel_squared = 2 * (energy + total_mass * ri);
    real v_t = angular_momentum * ri;
    real v_r = sqrt(fmax(0.0, rel_vel_squared - v_t * v_t));
               // it is assumed here that negative values for the argument
               // arise from round-off errors, and can be replaced by zero
    if (r_sin_true_an < 0) v_r = -v_r;

    pred_rel_vel = v_r * r_unit
	            + v_t * ri * (-r_sin_true_an * longitudinal_unit_vector
				  + r_cos_true_an * transverse_unit_vector);
}

void kepler::fast_pred_mean_anomaly_to_pos_and_vel()
{
    real cos_ecc_an, sin_ecc_an;
    fast_keplerseq(pred_mean_anomaly, eccentricity, kepler_fast_tol,
		   sin_ecc_an, cos_ecc_an);

    pred_separation = semi_major_axis * (1 - eccentricity * cos_ecc_an);
    real r_cos_true_an = semi_major_axis * (cos_ecc_an - eccentricity);

    // Writing b this way avoids sqrt(1 - eccentricity*eccentricity)!

    real semi_minor_axis = angular_momentum * period / (2*M_PI*semi_major_axis);
    real r_sin_true_an = semi_minor_axis * sin_ecc_an;

    fast_to_pred_rel_pos_vel(r_cos_true_an, r_sin_true_an);
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// Modified kepler member function (Steve, 6/01):

void kepler::pred_mean_anomaly_to_pos_and_vel()
{
    if (!kepler_fast_flag || energy >= 0 || angular_momentum == 0) {

	// Standard method, as used by kira:

	real ecc_an;

	// Determine the separation and true anomaly.

	if (energy != 0) {

	    if (energy < 0) {		// ellipse or linear bound orbit

		keplerseq(pred_mean_anomaly, eccentricity,
			  pred_true_anomaly, ecc_an);

		pred_separation = semi_major_axis 
		    			* (1 - eccentricity * cos(ecc_an));

	    } else {		 	// hyperbola or linear unbound orbit


		keplershypeq(pred_mean_anomaly, eccentricity,
			     pred_true_anomaly, ecc_an);

		pred_separation = semi_major_axis 
					* (eccentricity * cosh(ecc_an) - 1);
	    }

	    // Note: Use ecc_an for better results for very eccentric orbits.
	    //
	    // pred_separation = periastron * (1 + eccentricity)
	    //                    / (1 + eccentricity * cos(pred_true_anomaly));

	} else {				// parabola or linear orbit

	    if (angular_momentum > 0) {

		real tan_half_true_an = 
		    keplerspareq(pred_mean_anomaly, pred_true_anomaly);

		pred_separation = periastron
				    * (1 + tan_half_true_an*tan_half_true_an);

	    } else {	// linear case (note convention)

		pred_true_anomaly = pred_mean_anomaly;
		pred_separation = pow(fabs(pred_true_anomaly), 2/3.0);

	    }
	}

	// Calculate position and velocity.

	if (angular_momentum > 0)
	    to_pred_rel_pos_vel(cos(pred_true_anomaly), sin(pred_true_anomaly));
	else
	    to_pred_rel_pos_vel_linear(pred_true_anomaly);

    } else {

	// New faster function for bound orbits does everything in
	// a single call.  For use by the 4tree interpolation.

	fast_pred_mean_anomaly_to_pos_and_vel();	// <--- new ---
    }
}

//----------------------------------------------------------------------------

// Initialization:
// --------------

// Align_with_axes: initialize the right-handed triad (l, t, n) with the
//		    longitudinal unit vector l pointing along the specified
// 		    axis (note that x, y, z = 1, 2, 3, Piet).

void kepler::align_with_axes(int axis)
{
    if (axis == 1) {

	longitudinal_unit_vector = vec(1, 0, 0);
	transverse_unit_vector = vec(0, 1, 0);

    } else if (axis == 2) {

	longitudinal_unit_vector = vec(0, 1, 0);
	transverse_unit_vector = vec(0, 0, 1);

    } else {

	longitudinal_unit_vector = vec(0, 0, 1);
	transverse_unit_vector = vec(1, 0, 0);
    }
    
    normal_unit_vector = longitudinal_unit_vector ^ transverse_unit_vector;

}

// initialize_from_pos_and_vel: start with time, mass, position and velocity
//				and determine all other orbital properties.

void kepler::initialize_from_pos_and_vel(bool minimal, bool verbose)
{
    // Dynamics and geometry:
    // ---------------------

    energy = 0.5 * rel_vel * rel_vel - total_mass / separation;

    normal_unit_vector = rel_pos ^ rel_vel;
    angular_momentum = abs(normal_unit_vector);

    // Attempt to take care of special cases (eccentricity = 0 or 1,
    // energy = 1, etc.) at the time of initialization, to avoid messy
    // checks elsewhere int the package...

    if (energy != 0) {

        semi_major_axis = -0.5 * total_mass / energy;

        real rdotv = rel_pos * rel_vel;
        real temp = 1 - separation / semi_major_axis;
        eccentricity = sqrt(rdotv * rdotv / (total_mass * semi_major_axis)
                            + temp * temp);

	// Deal with rounding error, and avoid problems with nearly
	// circular or nearly linear orbits:

	if (eccentricity < ECC_TOL) eccentricity = 0;

	// Do we need an energy < 0 check here also?? (Steve, 8/02)

	if ((energy > 0 && eccentricity < 1 + ECC_TOL
	     		&& eccentricity > 1 - SQRT_TRIG_TOL)
	    || (angular_momentum > 0 && eccentricity == 1)) {

	    // Force a linear orbit.

	    eccentricity = 1;
	    angular_momentum = 0;

	} else if (eccentricity > 0
		   && eccentricity < circular_binary_limit) {

	    // Force a circular orbit (SPZ 2/98).

	    // PROBLEM: this may be OK for pragmatic reasons in an
	    // N-body simulation with real stars, but it is not OK
	    // in a scattering experiment where we may be interested
	    // in small eccentricity changes.  There is no mathematical
	    // reason to expect problems near eccentricity = 0.

	    // Fix by introducing circular_binary_limit, which specifies
	    // the maximum eccentricity that should be forced to 0 here.

	    if (verbose && (print_trig_warning
			    || eccentricity > 0.1*circular_binary_limit))
		cout << "kepler: binary with ecc = " << eccentricity
		     << " forced to be circular" << endl
		     << "        in kepler::initialize_from_pos_and_vel"
		     << endl;

	    eccentricity = 0;

	    // Tidally circularized binaries pose some problem in
	    // the transformation from hdyn to kepler.
	    // This might solve the problem, but I am not sure.
	    // There may be binaries that are almost circular, and
	    // these can still give the problem.
	    //                                            SPZ 2/98
	}

        // Convention: semi_major_axis is always > 0.

        semi_major_axis = fabs(semi_major_axis);
	periastron = semi_major_axis * fabs(1 - eccentricity);
	
	mean_motion = sqrt(total_mass / pow(semi_major_axis, 3));
	period = 2*M_PI / mean_motion;

    } else {

        eccentricity = 1;

        semi_major_axis = _INFINITY_;
        periastron = 0.5 * angular_momentum * angular_momentum / total_mass;

        mean_motion = sqrt(total_mass / 
			   (2 * periastron * periastron * periastron));

        // (The real mean motion is undefined, but this is useful in
        //  determining time from radius.)

        period = _INFINITY_;
    }

    if (angular_momentum != 0) 
        normal_unit_vector /= angular_momentum;
    else {
        eccentricity = 1;
	periastron = 0;
	if (energy == 0) mean_motion = sqrt(4.5*total_mass);	// yes, really!
    }


    if (minimal) return;	// "minimal" kepler has orbital elements,
				// but no phase or orientation information


    // Phase:
    // -----

    vec r_unit = rel_pos / separation;

    if (angular_momentum == 0) {
        vec temp = vec(1,0,0);  // construct an arbitrary normal vector.
        if (fabs(r_unit[0]) > 0.5) temp = vec(0,1,0);
        normal_unit_vector = r_unit ^ temp;
        normal_unit_vector /= abs(normal_unit_vector);
    }

    vec t_unit = normal_unit_vector ^ r_unit;

    real cos_true_an = 1, sin_true_an = 0;

    // Determine the true anomaly (note conventions in special cases).

    if (eccentricity > 0) {

	if (eccentricity != 1) {

	    cos_true_an = ((periastron/separation) * (1 + eccentricity) - 1)
				/ eccentricity;

	    // Note: The order of operations above is important for very
	    //       eccentric orbits!

	    // For almost circular orbits, the above expression is
	    // likely to be very inaccurate (sep ~ peri ~ a ==> tmp =
	    // 1), so expect check_trig_limit() to complain...

	    if (fabs(cos_true_an) > 1.0) {
		cout << "  rp, rsep, ecc =  " << periastron << " " 
		     << separation << " " << eccentricity << endl;
		// print_all();
	    }

	    check_trig_limit(this, cos_true_an, "set_phase #1");

	    // If we get here, cos_true_an is always in the correct range.

	    sin_true_an = cos_to_sin(cos_true_an);
	    if (rel_pos * rel_vel < 0) sin_true_an = -sin_true_an;

	    if (1 - cos_true_an < TRIG_TOL) { // more special treatment!
		true_anomaly = 0;
		cos_true_an = 1;
		sin_true_an = 0;
	    } else
		true_anomaly = acos(cos_true_an);

	    if (sin_true_an < 0) true_anomaly = -true_anomaly;

	} else {

	    if (angular_momentum > 0) {

		// Special case: see McCuskey, p. 54.

		true_anomaly = 2*acos(sqrt(periastron/separation));
		if (rel_pos * rel_vel < 0) true_anomaly = -true_anomaly;

		cos_true_an = cos(true_anomaly);
		sin_true_an = sin(true_anomaly);

	    } else {	// linear orbit: "true anomaly" = eccentric anomaly

		if (energy < 0) {

		    cos_true_an = 1 - separation / semi_major_axis;
		    check_trig_limit(this, cos_true_an, "set_phase #2");
		    true_anomaly = acos(cos_true_an);

		} else if (energy > 0) 

		    true_anomaly = acosh(1 + separation / semi_major_axis);

		else	// special case...

		    true_anomaly = pow(separation, 1.5);

		if (rel_pos * rel_vel < 0) true_anomaly = -true_anomaly;

		cos_true_an = -1;    // (to get the unit vectors right below)
		sin_true_an = 0;
	    }
	}
    }

    longitudinal_unit_vector = cos_true_an * r_unit - sin_true_an * t_unit;
    transverse_unit_vector = sin_true_an * r_unit + cos_true_an * t_unit;

    pred_true_anomaly = true_anomaly;
    pred_true_to_mean_anomaly();
    mean_anomaly = pred_mean_anomaly;

    mean_anomaly_to_periastron_passage();
}

// initialize_from_shape_and_phase: start with time, mass, semi-major axis,
//				    eccentricity, and mean anomaly, and
//		     		    determine all other orbital properties.
//
// Does not set unit vectors, so those must be specified beforehand in order
// for rel_pos and rel_vel to be defined.

// Convention: on input, specify semi-major axis > 0 and eccentricity.  If
//             eccentricity = 1, look at periastron to distinguish a parabola
//             from a linear orbit.  If the orbit is linear, use the sign of
//             the energy to distinguish a bound from an unbound orbit.

void kepler::initialize_from_shape_and_phase()
{
    // Compute energy, angular momentum, periastron, period, and mean
    // motion...

    if (eccentricity != 1) {		// elliptical/hyperbolic motion

	energy = 0.5 * total_mass / semi_major_axis;
	if (eccentricity < 1) energy = -energy;

        periastron = semi_major_axis * fabs(1 - eccentricity);
	mean_motion = sqrt(total_mass / pow(semi_major_axis, 3));
	period = (energy < 0 ? 2*M_PI / mean_motion : _INFINITY_);
	angular_momentum = sqrt(fabs(1 - eccentricity * eccentricity)
	                        * total_mass * semi_major_axis);
	    
    } else if (periastron == 0) {	// linear orbit

	angular_momentum = 0;

	if (energy == 0) {

	    period = semi_major_axis = _INFINITY_;
	    mean_motion = sqrt(4.5*total_mass);		// yes, really!

	} else {

	    energy = 0.5 * sign(energy) * total_mass / semi_major_axis;
	    mean_motion = sqrt(total_mass / pow(semi_major_axis, 3));
	    period = (energy < 0 ? 2*M_PI / mean_motion : _INFINITY_);
	}

    } else {				// parabola

        energy = 0;
	period = semi_major_axis = _INFINITY_;
	mean_motion = sqrt(total_mass / (2*pow(periastron, 3)));
	angular_momentum = sqrt(2 * total_mass * periastron);
    }

    // ...then compute phase.

    mean_anomaly_to_periastron_passage();
    pred_time = time;
    pred_mean_anomaly = mean_anomaly;

    pred_mean_anomaly_to_pos_and_vel();
    set_real_from_pred();
}

// initialize_from_integrals_and_separation:  start with time, mass, energy,
//				    angular momentum, and separation, and
//				    determine all other orbital properties.
//
// Does not set unit vectors, so those must be specified beforehand in order
// for rel_pos and rel_vel to be defined.

void kepler::initialize_from_integrals_and_separation(bool receding)
{

    // Compute semi-major axis, eccentricity, periastron, period, and
    // mean motion...

    if (angular_momentum == 0) {	// linear orbit

        eccentricity = 1;
	periastron = 0;

	if (energy == 0) {

	    semi_major_axis = period = _INFINITY_;
	    mean_motion = sqrt(4.5*total_mass);		// yes, really!

	} else {

	    semi_major_axis = 0.5 * total_mass / fabs(energy);
	    mean_motion = sqrt(total_mass / pow(semi_major_axis, 3));
	    period = (energy < 0 ? 2*M_PI / mean_motion : _INFINITY_);
	}

    } else {

	if (energy == 0) {		// parabolic motion

	    semi_major_axis = period = _INFINITY_;
	    eccentricity = 1;
	    periastron = pow(angular_momentum, 2) / (2 * total_mass);
	    mean_motion = sqrt(total_mass / (2*pow(periastron, 3)));

	} else {			// elliptical or hyperbolic motion

	    semi_major_axis = 0.5 * total_mass / fabs(energy);
	    mean_motion = sqrt(total_mass / pow(semi_major_axis, 3));
	    period = (energy < 0 ? 2*M_PI / mean_motion : _INFINITY_);
	    real ecc_sq = 1 + sign(energy) * pow(angular_momentum, 2)
				      / (total_mass*semi_major_axis);
	    if (ecc_sq < 0) {

		// Possible rounding error, possibly wrong input.

		if (ecc_sq < -1.e-6)
		    cout << "*** warning: inconsistent energy and angular"
			 << " momentum in kepler initialization"
			 << endl << flush;
		ecc_sq = 0;		// probably incorrect
	    }
	    eccentricity = sqrt(ecc_sq);
	    periastron = semi_major_axis * fabs(1 - eccentricity);
	}
    }

    // ...then compute phase.

    pred_advance_to_radius(separation);
    if ((receding && pred_true_anomaly < 0)
	 || (!receding && pred_true_anomaly > 0)) {

	// Picked the wrong part of the orbit.  Reverse angles.

      true_anomaly = -true_anomaly;
      mean_anomaly = -mean_anomaly;
    }
    pred_time = time;
    pred_mean_anomaly_to_pos_and_vel();
    set_real_from_pred();
    mean_anomaly_to_periastron_passage();
}

//----------------------------------------------------------------------------

// Transformation functions:
// ------------------------

real  kepler::pred_advance_to_periastron()     {

    if (eccentricity >= 1 && true_anomaly > 0)
	err_or_warn("pred_advance_to_periastron: outgoing hyperbolic orbit");

    pred_time = time_of_periastron_passage;

    if (eccentricity < 1) {

	// Original formula (incorrect if mean_motion is large):
        // pred_time = time_of_periastron_passage + 2*M_PI/mean_motion;

	// Next version (slow if mean_motion very large):
	// while (pred_time - time < 0.0) pred_time += 2*M_PI/mean_motion;

	// New implementation (SPZ&JM, 2/98):

	real n_orbits = 0;
	if (pred_time - time < 0.0) {
	    
	    real n_orbits = ceil((time-pred_time)/(2*M_PI/mean_motion));
	    pred_time += n_orbits * (2*M_PI/mean_motion);
	}
	if (pred_time < time) {
	    cout << "kepler::pred_advance_to_periastron(): shouldn't happen!"
		 << endl;
	    PRC(pred_time); PRC(time); PRC(n_orbits); PRL(pred_time-time);
	}
    }

    pred_true_anomaly = 0;
    pred_mean_anomaly = 0;
    pred_separation = periastron;

    to_pred_rel_pos_vel(1, 0);

    return pred_time;
}

real  kepler::pred_return_to_periastron()
    {
    if (eccentricity >= 1 && true_anomaly < 0)
	err_or_warn("pred_return_to_periastron: incoming hyperbolic orbit");

    pred_time = time_of_periastron_passage;
    pred_true_anomaly = 0;
    pred_mean_anomaly = 0;

    pred_separation = periastron;

    to_pred_rel_pos_vel(1, 0);

    return pred_time;
    }

real  kepler::pred_advance_to_apastron()
    {
    if (eccentricity >= 1)
	err_or_warn("pred_advance_to_apastron: hyperbolic orbit");

    if (true_anomaly < 0)
	pred_time = time_of_periastron_passage + 1.5 * 2*M_PI / mean_motion;
    else
	pred_time = time_of_periastron_passage + M_PI / mean_motion;

    pred_true_anomaly = M_PI;
    pred_mean_anomaly = M_PI;
    pred_separation = semi_major_axis * (1 + eccentricity);
    to_pred_rel_pos_vel(-1, 0);

    return pred_time;
    }

real  kepler::pred_return_to_apastron()
    {
    if (eccentricity >= 1)
	err_or_warn("pred_return_to_apastron: hyperbolic orbit");

    if (true_anomaly < 0)
	pred_time = time_of_periastron_passage + M_PI / mean_motion;
    else
	pred_time = time_of_periastron_passage - M_PI / mean_motion;

    pred_true_anomaly = M_PI;
    pred_mean_anomaly = M_PI;
    pred_separation = semi_major_axis * (1 + eccentricity);
    to_pred_rel_pos_vel(-1, 0);

    return pred_time;
    }

real kepler::pred_transform_to_radius(real r, int direction)
{
    // Transform forward or backwards in time, depending on direction.

    if (energy < 0)
        true_anomaly = sym_angle(true_anomaly); // unnecessary?
    else if (separation > r && direction*true_anomaly > 0) {
	err_or_warn(
	    "pred_transform_to_radius: unbound orbit inside target radius");
	if (direction > 0)
	  cout << "Mapping to incoming branch." << endl;
	else
	  cout << "Mapping to outgoing branch." << endl;
    }

    if (r < periastron) {
	if (kepler_tolerance_level <= 1)
	    err_exit("pred_transform_to_radius: r < periastron");
	r = periastron;
    }

    pred_separation = r;
    real cos_true_an = 1, sin_true_an = 0;

    // Determine the absolute magnitude of the true anomaly in all cases.

    if (eccentricity != 1) {

	if (eccentricity > 0 && r > 0)

  	    cos_true_an = ((periastron/r) * (1 + eccentricity) - 1)
		  	    / eccentricity;
	else
	    cos_true_an = 2; 	// force the trig check below

	if (check_trig_limit(this, cos_true_an,
			     "pred_transform_to_radius #1")) {

	    // Requested r was outside the allowed range.  Map directly
	    // to periastron or apastron.

	    if (cos_true_an > 0)
 	        return (direction == 1 ? pred_advance_to_periastron()
		                       : pred_return_to_periastron());
	    else
	        return (direction == 1 ? pred_advance_to_apastron()
		                       : pred_return_to_apastron());
	}

	pred_true_anomaly = acos(cos_true_an);
	sin_true_an = cos_to_sin(cos_true_an);

    } else {

	// Treat parabolic/linear motion as a special case.

	if (angular_momentum > 0) {

	    pred_true_anomaly = 2*acos(sqrt(fmin(1.0, periastron/r)));
	    cos_true_an = cos(pred_true_anomaly);
	    sin_true_an = sin(pred_true_anomaly);

	} else {	// linear motion

	    if (energy < 0) {

		cos_true_an = 1 - r / semi_major_axis;
		check_trig_limit(this, cos_true_an,
				 "pred_transform_to_radius #2");
		pred_true_anomaly = acos(cos_true_an);

	    } else if (energy > 0) 

	        pred_true_anomaly = acosh(1 + r / semi_major_axis);

	    else

	        pred_true_anomaly = pow(r, 1.5);
	}
    }

    // Determine the sign of the true anomaly.

    if (direction * (r - separation) < 0) {
	pred_true_anomaly = -pred_true_anomaly;
	if (angular_momentum > 0) sin_true_an = -sin_true_an;
    }

    // Calculate position and velocity...

    if (angular_momentum > 0)
	to_pred_rel_pos_vel(cos_true_an, sin_true_an);
    else        
        to_pred_rel_pos_vel_linear(pred_true_anomaly);

    // ...and the mean anomaly.

    pred_true_to_mean_anomaly();

    // Ensure the smallest change in mean anomaly in the right direction.

    if (energy < 0)
      while (direction * (pred_mean_anomaly - mean_anomaly) > 2*M_PI)
	pred_mean_anomaly -= direction*2*M_PI;

    // Update the predicted time.  Looks like mean_anomaly will remain
    // in the range (-M_PI, M_PI).  Make sure the time increases or decreases
    // appropriately, depending on direction.

    pred_time = time + (pred_mean_anomaly - mean_anomaly) / mean_motion;
    while (direction*(pred_time-time) < 0) pred_time += direction*period;

    return pred_time;
}

real  kepler::pred_advance_to_radius(real r) // transform forwards in time
{
    return pred_transform_to_radius(r, +1);
}

real  kepler::pred_return_to_radius(real r)  // transform backwards in time
{
    return pred_transform_to_radius(r, -1);
}

real  kepler::pred_transform_to_time(real t)
{
    pred_mean_anomaly = mean_anomaly + mean_motion * (t - time);
//    pred_mean_anomaly = mean_anomaly + mean_motion * (fmod(t - time,period));

    pred_time = t;
    pred_mean_anomaly_to_pos_and_vel();

    return pred_separation;
}

real  kepler::advance_to_radius(real r)           // transform forwards in time
    {
    time = pred_advance_to_radius(r);
    set_real_from_pred();

    return time;
    }

real  kepler::return_to_radius(real r)           // transform backwards in time
    {
    time = pred_return_to_radius(r);
    set_real_from_pred();

    return time;
    }

real  kepler::advance_to_periastron()
    {
    time = pred_advance_to_periastron();
    set_real_from_pred();

    return time;
    }

real  kepler::return_to_periastron()
    {
    time = pred_return_to_periastron();
    set_real_from_pred();

    return time;
    }

real  kepler::advance_to_apastron()
    {
    time = pred_advance_to_apastron();
    set_real_from_pred();

    return time;
    }

real  kepler::return_to_apastron()
    {
    time = pred_return_to_apastron();
    set_real_from_pred();

    return time;
    }

real  kepler::transform_to_time(real t)
    {
    separation = pred_transform_to_time(t);
    set_real_from_pred();

    return separation;
    }

// Handy functions:

// make_standard_kepler: set standard properties for an existing
//			 kepler structure.   By default, the new
//			 orbit is in the x-y plane, with the major
//			 axis in the x direction.

void make_standard_kepler(kepler &k, real t, real mass,
			  real energy, real eccentricity,
			  real q, real mean_anomaly,
			  int align_axis)
{
    // Note: energy is the true energy of the system;
    // 	     q is the true periastron only in the case of a parabolic orbit.
    //
    // In the event of inconsistent input (e.g. energy > 0, eccentricity < 1),
    // the internals of initialize_from_shape_and_phase() will silently force
    // the energy to follow the eccentricity (i.e. if eccentricity < 1, then
    // energy will be forced < 0).

    k.set_time(t);
    k.set_total_mass(mass);
    k.set_semi_major_axis((energy == 0 ? _INFINITY_ 
			               : 0.5*mass/fabs(energy)));
    k.set_eccentricity(eccentricity);
    k.set_periastron(q);
    k.set_energy(energy);
    k.set_mean_anomaly(mean_anomaly);
    if (align_axis >= 1 && align_axis <= 3)
	k.align_with_axes(align_axis);

    k.initialize_from_shape_and_phase();	// expects a, e [, q [, E]]
}

// set_random_orientation: randomly set the orientation of a kepler
// structure.  Orbit is random in 3-D if planar = 0, planar (x-y) with
// positive circulation if planar > 0, planar with negative
// circulation if planar < 0.

void set_random_orientation(kepler &k,
			    int planar)		    // default = 0
{
    if (planar > 1) planar = 1;
    if (planar < -1) planar = -1;

    real cos_theta = (real) planar, sin_theta = 0;

    if (planar == 0) {
	cos_theta = randinter(-1, 1);
	sin_theta = sqrt(1 - cos_theta * cos_theta);
    }

    real phi = randinter(0, 2*M_PI);
    real psi = randinter(0, 2*M_PI);

    // Construct the normal vector:

    vec n = vec(sin_theta*cos(phi), sin_theta*sin(phi), cos_theta);

    // Construct unit vectors a and b perpendicular to n:

    vec temp = vec(1, 0, 0);
    if (fabs(n[0]) > 0.5) temp = vec(0, 1, 0);	// temp is not parallel to n
    if (n[2] < 0) temp = -temp;

    vec b = n ^ temp;
    b /= abs(b);
    vec a = b ^ n;
    if (n[2] < 0) a = -a;	// force (a, b) to be (x, y) for n = +/-z
    
    // Construct *random* unit vectors l and t perpendicular to each
    // other and to n (psi = 0 ==> periastron along a):

    vec l = cos(psi)*a + sin(psi)*b;
    vec t = n ^ l;

    k.set_orientation(l, t, n);
    k.initialize_from_shape_and_phase();
}

void print_orbital_elements(hdyn *bi, hdyn *bj,
			    bool verbose)	// default = true
{
    // Compute and print the orbital elements of the binary with
    // components bi and bj.

    real time = bi->get_system_time();
    real mi = bi->get_mass();
    real mj = bj->get_mass();
    real mtot = mi + mj;
    kepler k;
    k.set_time(time);
    k.set_total_mass(mtot);
    k.set_rel_pos(bj->get_pos() - bi->get_pos());
    k.set_rel_vel(bj->get_vel() - bi->get_vel());
    k.initialize_from_pos_and_vel();
    // vec cmvel = (mi*bi->get_vel() + mj*bj->get_vel())/mtot;
    if (verbose)
	k.print_all();
    else {
	int p = cout.precision(10);
	cout << time << " ";
	cout.precision(9);
	cout << k.get_energy() << " "
	     << k.get_semi_major_axis() << " "
	     << k.get_eccentricity() << " "
	     << k.get_period() << " "
//	     << 0.5*square(cmvel)
	     << endl << flush;
	cout.precision(p);
    }
}

#else

#include "kepler.h"
kepler k;

static void tt(real t)
{
    cout << "  transform_to_time " << t << endl << endl;
    k.transform_to_time(t);
    k.print_dyn();
}

static void ar(real r)
{
    cout << "  advance_to_radius " << r << endl << endl;
    k.advance_to_radius(r);
    k.print_dyn();
}

static void rr(real r)
{
    cout << "  return_to_radius " << r << endl << endl;
    k.return_to_radius(r);
    k.print_dyn();
}

int main(int argc, char **argv)
{
    bool orbit_params = true;
    real time = 0;
    real semi = 1;
    real ecc = 0;
    real mass = 1;
    real m2 = 1;
    real mean_anomaly = 0;
    vec r = vec(1,0,0);
    vec v = vec(0,1,0);

    bool anim = false;
    bool to_apo = false;
    bool to_peri = false;
    real t_end = 0;
    real dt_end = 0;		// specify in orbit periods
    real dt_step = 0;		// specify in orbit periods
    real rad_sum = 0;

    int i = 0;
    while (++i < argc)
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {

		case 'a':	semi = atof(argv[++i]);
				orbit_params = true;
				break;
		case 'A':	to_apo = 1 - to_apo;
				anim = to_apo;
				if (anim) to_peri = false;
				break;
		case 'd':	dt_step = atof(argv[++i]);
				anim = true;
				break;
		case 'e':	ecc = atof(argv[++i]);
				break;
		case 'm':	if (argv[i][2] == '2') {
				    m2 = atof(argv[++i]);
				    anim = true;
				} else
				    mass = atof(argv[++i]);
				break;
		case 'M':	mean_anomaly = atof(argv[++i]);
				break;
		case 'O':	anim = true;
				break;
		case 'P':	to_peri = 1 - to_peri;
				anim = to_peri;
				if (anim) to_apo = false;
				break;
		case 'r':
		case 'R':	rad_sum = atof(argv[++i]);
				anim = true;
				break;
		case 't':	time = atof(argv[++i]);
				break;
		case 'T':	dt_end = atof(argv[++i]);
				anim = true;
				break;
		case 'x':	for (int k = 0; k < 3; k++)
		    		    r[k] = atof(argv[++i]);
				orbit_params = false;
				break;
		case 'v':	for (int k = 0; k < 3; k++)
		    		    v[k] = atof(argv[++i]);
				break;
		default:	cerr << "unknown option \"" << argv[i]
				     << endl << flush;
				exit(1);
	    }

    set_kepler_tolerance(2);
    cout.precision(12);
    cerr.precision(12);

    real m1 = 0;
    if (anim) {
	m1 = mass;
	if (m2 <= 0) m2 = m1;
	mass = m1 + m2;
    }

    k.set_time(time);
    k.set_total_mass(mass);

    if (orbit_params) {

	k.set_semi_major_axis(semi);
	k.set_eccentricity(ecc);
	k.set_mean_anomaly(mean_anomaly);
	k.align_with_axes(1);
	k.initialize_from_shape_and_phase();

	if (!anim)
	    cout << "Initialized from shape and phase" << endl;

    } else {

	k.set_rel_pos(r);
	k.set_rel_vel(v);
	k.initialize_from_pos_and_vel();

	if (!anim)
	    cout << "Initialized from pos and vel" << endl;

    }

    if (anim) {

	if (dt_end == 0) dt_end = 1;
	t_end = time + k.get_period()*dt_end;

	real t_peri = k.get_time_of_periastron_passage();
	while (t_peri <= time) t_peri += k.get_period();

	cout << "start: " << endl;

	if (to_peri) {
	    t_end = t_peri;
	    PRI(4); PRL(t_peri);
	}

	if (to_apo) {
	    real t_apo = t_peri - 0.5*k.get_period();
	    while (t_apo < time) t_apo += k.get_period();
	    t_end = t_apo;
	    PRI(4); PRL(t_apo);
	}

	PRI(4); PRC(time); PRC(k.get_period()); PRL(t_end);

	if (dt_step == 0) dt_step = 1./64;
	dt_step *= k.get_period();

	// Make the two-body system.

	hdyn *b  = new hdyn();
	hdyn *bo = new hdyn();
	hdyn *by = new hdyn();

	b->set_oldest_daughter(bo);
	bo->set_parent(b);
	bo->set_index(1);
	bo->set_mass(m1);
	bo->set_younger_sister(by);
	by->set_parent(b);
	by->set_index(2);
	by->set_mass(m2);
	by->set_older_sister(bo);

	while (time <= t_end + 0.01*dt_step) {
	    k.transform_to_time(time);
	    if (k.get_separation() < rad_sum) break;

	    vec pos = k.get_rel_pos();
	    vec vel = k.get_rel_vel();

	    bo->set_pos(-(m2/mass)*pos);
	    bo->set_vel(-(m2/mass)*vel);
	    by->set_pos( (m1/mass)*pos);
	    by->set_vel( (m1/mass)*vel);

	    time += dt_step;
	}

	cout << "end: " << endl;
	PRI(4); PRC(k.get_time()); PRL(k.get_separation());

	delete by;
	delete bo;
	delete b;

    } else {

	cout << endl;
	k.print_all();

	// Loop over user input.

	for (;;) {

	    char opt;

	    cout << endl << "? ";
	    opt = 'q';
	    cin >> opt;
	    if (opt == 'q') exit(0);

	    if (opt == 'p') {

		cout << endl;
		k.print_all();

	    } else if (opt == 'h' || opt == '?') {

		cout << "  options:   a R    advance to radius R" << endl
		     << "             d dT   increment time by dT" << endl
		     << "             p      print" << endl
		     << "             q      quit" << endl
		     << "             r R    return to radius R" << endl
		     << "             t T    evolve to time T" << endl;

	    } else {

		real x;
		cin >> x;

		if (opt == 'd')
		    tt(k.get_time()+x);
		else if (opt == 't')
		    tt(x);
		else if (opt == 'a')
		    ar(x);
		else if (opt == 'r') 
		    rr(x);
	    }
	}

    }
    return 0;
}

#endif
