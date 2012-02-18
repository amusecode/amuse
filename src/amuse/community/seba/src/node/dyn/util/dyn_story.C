
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// dyn_story.C:  initialization of system parameters from input snapshots.
//
// Externally visible functions:
//
//	real get_initial_mass		// functions get_initial_mass and
//	real get_initial_virial_radius	// get_initial_virial radius are
//					// essentially identical
//
//	real get_initial_jacobi_radius	// similar to the above in general
//					// approach, but more complex logic
//
//	void set_tidal_params		// set tidal parameters from the
//					// input story
//	void test_tidal_params		// check consistency of disk fields
//
//	void check_set_tidal		// check and set tidal parameters
//	void check_set_plummer		// check and set Plummer parameters
//	void check_set_power_law	// check and set power-law parameters
//
//	void check_set_external		// check and set all external fields

// Created Aug/Sep 2001, Steve McMillan

#include "dyn.h"

// Require that all systems have known initial mass and initial virial
// radius.  These may be set by the make_* functions or by the function
// scale.  In general, kira or other applications should *not* compute
// these quantities -- better to determine them in advance.

// get_initial_mass:  Read the initial_mass variable from the log story.
//		      If not specified, use the input value, if any.
//		      Otherwise, use a default.
//
//		      Write the new variable to the log story in the latter
//		      two cases.

static bool mass_msg = true;

real get_initial_mass(dyn* b,
		      bool verbose,			// default = false
		      bool mass_set,			// default = false
		      real input_mass)			// default = 0
{
    real initial_mass = getrq(b->get_log_story(), "initial_mass");

    if (initial_mass > 0) {

	// The input snapshot has the information needed.  Snapshot
	// data OVERRIDES any input value.

	if (verbose && mass_msg)
	    cerr << "read initial_mass = " << initial_mass
		 << " from input snapshot" << endl;

    } else {

	// See if input_mass exists, or use a default.

	if (mass_set) {

	    if (verbose && mass_msg) {

		cerr << "setting initial mass = " << input_mass << endl;

		if ((real)b->get_system_time() > 0)
		    cerr << endl
			 << "Warning: setting initial mass with time > 0"
			 << endl;
	    }

	    initial_mass = input_mass;

	} else {

	    if ((real)b->get_system_time() <= 0)

		// Simply compute the mass if t <= 0.

		initial_mass = total_mass(b);

	    else {

		// Two options here:
		//
		//	1. Return a negative number, requiring the user
		//	   to provide the necessary data.
		//
		//	2. Adopt a default mass of 1, printing a message
		//	   so the user can use "make..." or "scale" to do
		//	   this properly.
		//
		// Choose (1) for now.

#if 1
		initial_mass = -1;

		if (verbose && mass_msg)
		    cerr << "get_initial_mass:  initial_mass unknown"
			 << " (set with scale)" << endl;
#else
		initial_mass = 1;

		if (verbose && mass_msg)
		    cerr << "adopting default initial_mass = "
			 << initial_mass << endl;
#endif
	    }
	}

	// Write the initial mass to the log story.

	if (initial_mass > 0)
	    putrq(b->get_log_story(), "initial_mass", initial_mass,
		  HIGH_PRECISION);
    }

    mass_msg = false;
    return initial_mass;
}

// get_initial_virial_radius:  Read the initial_virial_radius variable from
//			       the log story.
//			       If not specified, use the input value, if any.
//		      	       Otherwise, use a default.
//
//		      	       Write the new variable to the log story in the
//			       latter two cases.

static bool rvir_msg = true;

real get_initial_virial_radius(dyn* b,
			       bool verbose,		// default = false
			       bool r_virial_set,	// default = false
			       real input_r_virial)	// default = 0
{
    real r_virial = getrq(b->get_log_story(), "initial_rvirial");

    if (r_virial > 0) {

	// The input snapshot has the information needed.  Snapshot
	// data OVERRIDES any input value.

	if (verbose && rvir_msg)
	    cerr << "read initial r_virial = " << r_virial
		 << " from input snapshot" << endl;

    } else {

	// See if input_r_virial exists, or use a default.

	if (r_virial_set) {

	    if (verbose && rvir_msg) {

		cerr << "setting initial r_virial = " << input_r_virial << endl;

		if ((real)b->get_system_time() > 0)
		    cerr << endl
			 << "Warning: setting initial r_virial with time > 0"
			 << endl;
	    }

	    r_virial = input_r_virial;

	} else {

	    // Two options here:
	    //
	    //	1. Return a negative number, requiring the user
	    //	   to provide the necessary data.
	    //
	    //	2. Adopt a default radius of 1, printing a message
	    //	   so the user can use "make..." or "scale" to do
	    //	   this properly.
	    //
	    // Choose (1) for now.

#if 1
	    r_virial = -1;

	    if (verbose && rvir_msg)
		cerr << "warning: get_initial_virial_radius: "
		     << "initial virial radius cannot be determined" << endl
		     << "         (set before use by running scale with no"
		     << " parameters)" << endl;
#else
	    r_virial = 1;

	    if (verbose && rvir_msg)
		cerr << "adopting default initial_rvirial = "
		     << r_virial << endl;
#endif
	}

	// Write the virial radius to the log story.

	if (r_virial > 0)
	    putrq(b->get_log_story(), "initial_rvirial", r_virial);
    }

    rvir_msg = false;
    return r_virial;
}

// get_initial_jacobi_radius:  Read the initial_jacobi_radius variable from
//			       the log story.
//			       If not specified, use the input value, if any,
//			       modified in various ways.  There is no default.
//
//			       For use by kira, if specified, the variable
//			       kira_initial_jacobi_radius takes precedence
//			       over all other possibilities.
//
//		      	       Write the new variable to the log story if
//			       a value can be determined.

real get_initial_jacobi_radius(dyn* b,
			       real r_virial,
			       bool verbose,		// default = false
			       bool r_jacobi_set,	// default = false
			       real input_r_jacobi)	// default = 0
{
    real initial_r_jacobi = -1;

    // Determine the initial Jacobi radius of the system.  This radius
    // may be input as an argument, or it may be derived from the
    // information in the input snapshot.  Unlike r_virial, the input
    // value takes precedence over initial data found in the input
    // snapshot, although it is superceded by any "kira" data found
    // in that file.  The convention is such that, if a non-kira
    // version of initial_r_jacobi is found in the input snapshot, the
    // input value will *scale* that number.  A kira version is used
    // as is.  Otherwise, the input value scales the virial radius
    // (hence the two tests of r_jacobi_set below).
    //
    // (Messy logic actually makes it possible to control the
    //  Jacobi radius in a fairly natural way...)

    // See if a kira_initial_jacobi_radius exists in the input file.
    // If it does, it TAKES PRECEDENCE over any other setting.

    if (find_qmatch(b->get_log_story(), "kira_initial_jacobi_radius")) {

	real kira_initial_jacobi_radius
	    = getrq(b->get_log_story(), "kira_initial_jacobi_radius");

	if (kira_initial_jacobi_radius > 0) {

	    initial_r_jacobi = kira_initial_jacobi_radius;

	    if (verbose) {
		cerr << "using initial Jacobi radius ";
		PRL(kira_initial_jacobi_radius);
		cerr << "    from input snapshot" << endl;

		if (r_jacobi_set)
		    cerr << "ignoring input value " << input_r_jacobi << endl;
	    }

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "kira_initial_jacobi_radius "
		     << "from input snapshot"
		     << endl;
	    else
		err_exit(
	        "Error reading kira_initial_jacobi_radius from input snapshot");

	}
    }

    if (initial_r_jacobi < 0) {

	// No kira Jacobi radius known.  See if we can determine the
	// system's initial Jacobi radius from the input snapshot.

	if (find_qmatch(b->get_log_story(), "initial_rtidal_over_rvirial")) {

	    // The input snapshot contains an initial "tidal" radius.  Most
	    // likely the initial model is a King profile and this is the
	    // King radius (misnamed for historical reasons -- it may or
	    // may not have anything to do with a tidal cutoff).
	    // Alternatively, the ratio may have been set by add_tidal.

	    real r_jacobi_over_r_virial = getrq(b->get_log_story(),
						"initial_rtidal_over_rvirial");

	    if (r_jacobi_over_r_virial > 0) {

		initial_r_jacobi = r_jacobi_over_r_virial * r_virial;

		if (verbose) {
		    cerr << "got r_jacobi_over_r_virial = "
			 << r_jacobi_over_r_virial
			 << " from input snapshot"
			 << endl;
		    if (initial_r_jacobi > 0) {
			cerr << "    setting ";
			PRL(initial_r_jacobi);
		    }
		}

	    } else {

		if (verbose)
		    cerr << "Warning: error reading "
			 << "initial_rtidal_over_rvirial from input snapshot"
			 << endl;
		else
		    err_exit(
	      "Error reading initial_rtidal_over_rvirial from input snapshot");
	    }
	}

	// See if there was any input value, and resolve any conflicts.

	if (r_jacobi_set) {

	    // An input value was specified.

	    if (initial_r_jacobi > 0) {

		// The Jacobi radius has been specified both in the input file
		// and as an input argument.  If the model is an anisotropic
		// King model, we must use the snapshot data and ignore the
		// argument.  Otherwise, use the input argument to scale the
		// "known" value.

		// The variable alpha3_over_alpha1 is set *only* by
		// make_aniso_king.

		if (find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {

		    if (verbose)
			cerr << "ignoring input Jacobi radius ("
			     << input_r_jacobi << ") for anisotropic King model"
			     << endl;

		} else {

		    if (verbose)
			cerr << "input Jacobi radius ("
			     << input_r_jacobi << ") scales initial value "
			     << initial_r_jacobi
			     << endl
			     << "    from input snapshot"
			     << endl;

		    initial_r_jacobi *= input_r_jacobi;

		}

	    } else {

		// Use the input data to scale the virial radius.

		if (verbose)
		    cerr << "input Jacobi radius ("
			 << input_r_jacobi << ") scales initial virial radius "
			 << r_virial << endl;

		initial_r_jacobi = input_r_jacobi * r_virial;

	    }
	}

	// Save the kira Jacobi radius for restart.

	if (initial_r_jacobi > 0)
	    putrq(b->get_log_story(), "kira_initial_jacobi_radius",
		  initial_r_jacobi);

    }

    return initial_r_jacobi;
}

#define MATCH_TOL 0.001

local bool matches(real r, real v)
{
    return (abs(r/v - 1) < MATCH_TOL);
}

static const char* tidal_type[5] = {"none",
				    "point-mass",
				    "isothermal",
				    "disk",
				    "custom"};

// Express Galactic parameters in "Stellar" units (see Mihalas 1968):
//
//	G		=  1
//	length unit	=  1 pc
//	velocity unit	=  1 km/s
//
// ==>	time unit 	=  0.978 Myr
//	mass unit	=  232 Msun

#define OORT_A	(0.0144)	// km/s/pc
#define OORT_B	(-0.012)	// km/s/pc
#define RHO_G	(0.11/232)	// (232 Msun)/pc^3

#define Rsun_pc 2.255e-8	// R_sun/1 parsec = 6.960e+10/3.086e+18;

#define OORT_ALPHA_RATIO \
	((4*M_PI*RHO_G + 2*(OORT_A*OORT_A - OORT_B*OORT_B)) \
			 / (-4*OORT_A*(OORT_A-OORT_B)))

local void tidal_msg(int i, real alpha3_over_alpha1)
{
    cerr << "forcing tidal_field_type = " << i
	 << " for anisotropic King model "
	 << endl
	 << "    with alpha3_over_alpha1 = "
	 << alpha3_over_alpha1
	 << endl;
}

void set_tidal_params(dyn* b,
		      bool verbose,
		      real initial_r_jacobi,
		      real initial_mass,
		      int& tidal_field_type)
{
    // Set the tidal parameters for the system.
    //
    // Procedure:
    //		kira_tidal_field_type = 0 suppresses tidal fields
    //		if no type can be determined, suppress tidal field
    //		initial_mass and initial_r_jacobi set the value of alpha1
    //		tidal_field_type then sets the value of alpha3
    //
    // NOTE that tidal_field_type = 3, the field is intended to
    // model the Galactic field in the solar neighborhood, for
    // which alpha1 and alpha3 are actually determined directly by
    // the local Oort constants.  The resultant field may thus *not*
    // be consistent with the choice of radius used later in
    // establishing physical scales -- see test_tidal_params().
    //
    // Tidal field type (kira_tidal_field_type) or anisotropic King model
    // initial info (alpha3_over_alpha1) in the input snapshot OVERRIDE
    // the input data, if any.
    //
    // Thus, to decipher the tidal parameters, we need one of:
    //
    //		kira_tidal_field_type  ==> sets tidal_field_type directly and
    //					   hence a standard alpha3_over_alpha1
    //	or	alpha3_over_alpha1     ==> allows inference of tidal_field_type
    //
    // Then the tidal parameters are set from the mass and Jacobi radius.
    // To determine the Jacobi radius, we need:
    //
    //		kira_initial_jacobi_radius
    //	or	initial_r_virial  and  initial_rtidal_over_rvirial
    //
    // These come from make_*, make_king, make_aniso_king, add_tidal,
    // and scale.

    int kira_tidal_field_type = -1;

    if (find_qmatch(b->get_log_story(), "kira_tidal_field_type")) {

	kira_tidal_field_type
	    = getiq(b->get_log_story(), "kira_tidal_field_type");

	if (kira_tidal_field_type == 0)

	    return;					// take no action

	else if (kira_tidal_field_type > 0) {

	    if (verbose) {
		cerr << "tidal_field_type = " << kira_tidal_field_type
		     << " (" << tidal_type[kira_tidal_field_type]
		     << ") from input snapshot" << endl;

		if (tidal_field_type > 0)
		    cerr << "ignoring input tidal_field_type"
			 << tidal_field_type << endl;
	    }

	    tidal_field_type = kira_tidal_field_type;

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "kira_tidal_field_type "
		     << "from input snapshot"
		     << endl;
	    else
		err_exit(
	     "Error reading kira_tidal_field_type from input snapshot");

	}

    } else if (find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {

	// Try to infer tidal_field_type from alpha3_over_alpha1.

	real alpha3_over_alpha1
	    = getrq(b->get_log_story(), "alpha3_over_alpha1");

	if (alpha3_over_alpha1 < 0
	    && alpha3_over_alpha1 > -VERY_LARGE_NUMBER) {

	    // See if the input ratio matches any of the standard types.
	    // If it does, use that type; if not, flag a warning and use
	    // the input type or the default.

	    if (matches(alpha3_over_alpha1, -1.0/3)) {
		if (tidal_field_type != 1) {
		    tidal_field_type = 1;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    } else if (matches(alpha3_over_alpha1, -1.0/2)) {
		if (tidal_field_type != 2) {
		    tidal_field_type = 2;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    }  else if (matches(alpha3_over_alpha1, OORT_ALPHA_RATIO)) {
		if (tidal_field_type != 3) {
		    tidal_field_type = 3;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    } else {
		cerr << "Warning: snapshot alpha3_over_alpha1 = "
		     << alpha3_over_alpha1
		     << " does not match any standard value"
		     << endl;
	    }

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "alpha3_over_alpha1 "
		     << "from input snapshot"
		     << endl;

	}
    }

    if (tidal_field_type > 4)
	err_exit("Unknown tidal field type");

    if (verbose) {

	if (tidal_field_type <= 0)

	    cerr << "unable to determine tidal field type" << endl;

	else {

	    cerr << "using ";

	    if (tidal_field_type <= 1) {

		cerr << tidal_type[1] << " ";
		if (tidal_field_type <= 0) cerr << "(default) ";

	    } else
		cerr << tidal_type[tidal_field_type] << " ";

	    cerr << "tidal field; "
		 << "initial Jacobi radius = " << initial_r_jacobi
		 << endl;
	}
    }

    // Set up the dynamical tidal parameters.

    real alpha1, alpha3, omega;
    if (tidal_field_type > 0)
	alpha1 = -initial_mass / pow(initial_r_jacobi, 3.0);	// (definition)

    // Note that we don't use alpha3_over_alpha1, even if it is stored
    // in the snapshot.  We use the "standard" ratios, or rederive the
    // ratio from the stored Oort constants.

    if (tidal_field_type == 1) {

	// Point-mass field (see Binney & Tremaine, p. 452).

	alpha3 = -alpha1/3;
	omega = sqrt(alpha3);

    } else if (tidal_field_type == 2) {

	// Isothermal halo.

	alpha3 = -alpha1/2;
	omega = sqrt(alpha3);
	
    } else if (tidal_field_type == 3) {

	// Disk field.  Use the local Oort constants (defined above).
	// Conversion between Oort constants and alpha1/3 is taken
	// from Heggie & Ramamani (MNRAS 272, 317, 1995):
	//
	//	alpha1 = -4 A (A - B)
	//	alpha3 = 4 PI G RHO_G + 2 (A^2 - B^2)
	//
	// and recall that G = 1 by our above choice of units.

	alpha3 = alpha1 * OORT_ALPHA_RATIO;

	// Get omega from the definition of A and B:
	//
	//	A =  (1/2) (Vc/R - dVc/dR)
	//	B = -(1/2) (Vc/R + dVc/dR)
	//
	// ==>	omega = Vc/R = A - B,  dVc/dR = A + B
	//
	// so	alpha1 = -2 omega (omega - dVc/dR)
	//	       = -2 omega^2 (1 - (A + B)/(A - B))
	//	       =  4 omega^2 B / (A - B)

	omega = sqrt(0.25 * alpha1 * (OORT_A/OORT_B - 1));

    } else if (tidal_field_type == 4) {

	// Custom tidal field.  Require that the original physical
	// parameters be saved in the input snapshot.

	real Oort_A = 0;
	real Oort_B = 0;
	real rho_G = 0;

	if (find_qmatch(b->get_log_story(), "Oort_A_constant"))
	    Oort_A = getrq(b->get_log_story(), "Oort_A_constant");
	else
	    err_exit("Oort A constant not specified");

	if (find_qmatch(b->get_log_story(), "Oort_B_constant"))
	    Oort_B = getrq(b->get_log_story(), "Oort_B_constant");	
	else
	    err_exit("Oort B constant not specified");

	if (find_qmatch(b->get_log_story(), "local_mass_density"))
	    rho_G = getrq(b->get_log_story(), "local_mass_density");
	else
	    err_exit("rho_G not specified");

	cerr << "create custom tidal field from " << endl;
	PRI(4);PRC(Oort_A);PRC(Oort_B);PRL(rho_G);

	if (Oort_A != 0 && Oort_B != 0 && rho_G != 0) {

            // alpha1 = -4*Oort_A*(Oort_A-Oort_B);	// no!

	    real alpha3_over_alpha1 =
		(4*M_PI*rho_G + 2*(pow(Oort_A, 2) - pow(Oort_B, 2)))
		    / (-4*Oort_A*(Oort_A - Oort_B));

	    alpha3 = alpha1 * alpha3_over_alpha1;
	    omega = sqrt(0.25*alpha1 * (Oort_A/Oort_B - 1));
	}
	else
	  err_exit("Insufficient information to reconstruct tidal field");
    }

    if (verbose && tidal_field_type > 0) {
	PRI(4); PRC(alpha1); PRC(alpha3); PRL(omega);
    }

    if (tidal_field_type > 0) {
	b->set_tidal_field(tidal_field_type);
	b->set_omega(omega);
	b->set_alpha(alpha1, alpha3);
    }

    // Save the field type information in the root log story, if necessary.

    if (kira_tidal_field_type <= 0)
	putiq(b->get_log_story(), "kira_tidal_field_type",
	      tidal_field_type);

}

void test_tidal_params(dyn* b,
		       bool verbose,
		       real initial_r_jacobi,
		       real initial_r_virial,
		       real initial_mass)
{
    cerr << endl << "comparing initial parameters for disk tidal field:"
	 << endl;

    fprintf(stderr, "    model r_virial = %.3f, r_jacobi = %.3f",
	    initial_r_virial, initial_r_jacobi);
    fprintf(stderr, ", ratio = %.3f\n", initial_r_jacobi/initial_r_virial);

    // Convert initial mass and virial radius using conversion factors.
    // Compute Jacobi radius from Oort constants.

    initial_mass
	= b->get_starbase()->conv_m_dyn_to_star(initial_mass);	    // Msun
    initial_r_virial
	= b->get_starbase()->conv_r_dyn_to_star(initial_r_virial)
	    * Rsun_pc;						    // pc

    initial_r_jacobi = pow((initial_mass/232)
			    / (4*OORT_A*(OORT_A-OORT_B)), 1.0/3);   // pc

    fprintf(stderr, "    real  r_virial = %.3f, r_jacobi = %.3f",
	    initial_r_virial, initial_r_jacobi);
    fprintf(stderr, ", ratio = %.3f\n", initial_r_jacobi/initial_r_virial);
}

int check_set_tidal(dyn *b,
		    bool verbose)		// default = false
{
    int tidal_field_type = 0;
    b->set_tidal_field(tidal_field_type);

    real initial_mass = get_initial_mass(b, verbose);
    real initial_r_virial = get_initial_virial_radius(b, verbose);
    real initial_r_jacobi = get_initial_jacobi_radius(b, initial_r_virial,
						      verbose);
    if (initial_r_jacobi > 0)
	set_tidal_params(b, verbose,
			 initial_r_jacobi,
			 initial_mass,
			 tidal_field_type);
    else if (tidal_field_type > 0) {
	if (verbose)
	    cerr << "warning: check_set_tidal: jacobi radius < 0;"
		 << " suppressing tidal field" << endl;
	tidal_field_type = 0;
    }

    return tidal_field_type;
}

void check_set_plummer(dyn *b,
		       bool verbose)		// default = false
{
    // Check for Plummer-field data in the input snapshot, and
    // set kira parameters accordingly.  No coordination with the
    // command line because these parameters can *only* be set via
    // the input snapshot.  (Steve, 7/01)
    //
    // Expected fields in the input file are
    //
    //		kira_plummer_mass	=  M
    //		kira_plummer_scale	=  R
    //		kira_plummer_center	=  x y z

     if (find_qmatch(b->get_log_story(), "kira_plummer_mass")
	 && find_qmatch(b->get_log_story(), "kira_plummer_scale")) {

	 b->set_plummer();

	 b->set_p_mass(getrq(b->get_log_story(), "kira_plummer_mass"));
	 b->set_p_scale_sq(pow(getrq(b->get_log_story(),
				     "kira_plummer_scale"), 2));

	 vec center = 0;
	 if (find_qmatch(b->get_log_story(), "kira_plummer_center"))
	     center = getvq(b->get_log_story(), "kira_plummer_center");
	 b->set_p_center(center);

	 if (verbose) {
	     real M = b->get_p_mass();
	     real a = sqrt(b->get_p_scale_sq());
	     vec center = b->get_p_center();
	     cerr << "check_set_plummer:  "; PRC(M); PRC(a); PRL(center);
	 }
     }
}

void check_set_power_law(dyn *b,
			 bool verbose)		// default = false
{
    // Check for power-law-field data in the input snapshot, and
    // set kira parameters accordingly.  No coordination with the
    // command line because these parameters can *only* be set via
    // the input snapshot.  (Steve, 8/01)
    //
    // Expected fields in the input file are
    //
    //		kira_pl_coeff		=  A
    //		kira_pl_scale		=  R
    //		kira_pl_exponent	=  e
    //		kira_pl_center		=  x y z

     if (find_qmatch(b->get_log_story(), "kira_pl_coeff")
	 && find_qmatch(b->get_log_story(), "kira_pl_scale")) {

	 b->set_pl();

	 b->set_pl_coeff(getrq(b->get_log_story(), "kira_pl_coeff"));
	 b->set_pl_scale(getrq(b->get_log_story(), "kira_pl_scale"));

	 real exponent = 0;
	 if (find_qmatch(b->get_log_story(), "kira_pl_exponent"))
	     exponent = getrq(b->get_log_story(), "kira_pl_exponent");
	 b->set_pl_exponent(exponent);

	 vec center = 0;
	 if (find_qmatch(b->get_log_story(), "kira_pl_center"))
	     center = getvq(b->get_log_story(), "kira_pl_center");
	 b->set_pl_center(center);

	 if (verbose) {
	     real A = b->get_pl_coeff();
	     real a = b->get_pl_scale();
	     real x = b->get_pl_exponent();
	     vec center = b->get_pl_center();
	     cerr << "check_set_power_law:  ";
	     PRC(A); PRC(a); PRC(x); PRL(center);
	 }
     }
}

void check_set_external(dyn *b,
			bool verbose,		// default = false
			int fric_int)		// default = -1
{
    // Check for and set parameters for all known external fields.

    check_set_tidal(b, verbose);
    check_set_plummer(b, verbose);
    check_set_power_law(b, verbose);

    // See whether or not to apply dynamical friction.  This really only
    // applies to Plummer models, for now, but check here because the
    // setting is controlled by a command-line argument.
    // Specifying fric_int on the command line will override any settings
    // found in the root log story, so long as the option is legal.

    if (find_qmatch(b->get_log_story(), "kira_plummer_friction")) {

	// A friction flag exists in the input data.  Use it unless
	// the command line supercedes it.

	int fric_int_snap = getiq(b->get_log_story(), "kira_plummer_friction");
	if (fric_int < 0)
	    fric_int = fric_int_snap;
	else
	    cerr << "command-line flag -f " << fric_int
		 << (fric_int == fric_int_snap ? " is identical to"
		     			     : " overrides")
		 << " value in input snapshot" << endl;
    }

    if (fric_int > 0) fric_int = 1;

    if (fric_int > 0 && !b->get_plummer()) {
	cerr << "check_set_external: ignoring '-f' option"
	     << " -- no Plummer field" << endl;
	fric_int = -1;
    }

    if (fric_int >= 0) {

	// Set the log story friction flag and dyn static data.

	cerr << "turning " << (fric_int == 0 ? "off" : "on")
	     << " internal dynamical friction " << endl;
	putiq(b->get_log_story(), "kira_plummer_friction", fric_int);
	b->set_p_friction(fric_int);
    }
}

void check_set_ignore_internal(dyn *b,
			       bool verbose)		// default = false
{
    // Check for and set flag to skip internal interactions.

    bool ignore_internal = false;
    if (find_qmatch(b->get_log_story(), "ignore_internal")
	&& getiq(b->get_log_story(), "ignore_internal") == 1)
	ignore_internal = true;

    b->set_ignore_internal(ignore_internal);
}
