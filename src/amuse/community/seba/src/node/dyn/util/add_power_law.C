
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add power-law parameters to an input snapshot and write it out
//// again.  Do not change the input data.  External field is basically
//// that of a power-law mass distribution
////                 M(r)  =  A r^x
//// except that the density is constant for r < a.  Parameters are 
//// interpreted in physical units if they have been specified with 
//// add_star, are converted to N-body units if necessary, and are 
//// written to the root log story for subsequent use by set_com,
//// scale, kira, etc.  Note that -e 0 is currently implemented as a 
//// Plummer field, and is flagged as such.
////
//// Usage: add_power_law [OPTIONS] < input > output
////
//// Options:      
////               -A/M  specify coefficient [1]
////               -c    add comment [none]
////               -C    specify center [(0,0,0)]
////               -e/E/x/X  specify exponent [0]
////               -a/R  specify scale [1]
////               -G    select parameters (physical units) appropriate
////                     for the Galactic center (Mezger et al. 1999) [no]
////               -n    force interpretation of all parameters in
////                     N-body units [no]
////
//// Written by Steve McMillan.
////
//// Report bug to starlab@sns.ias.edu.


//   version 1:  Aug/Sep 2001   Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

bool get_physical_scales(dyn *b, real& mass, real& length, real& time)
{
    mass = length = time = -1;

    if (b)
	if (b->get_starbase()) {

#define Rsun_pc 2.255e-8	// R_sun/1 parsec = 6.960e+10/3.086e+18;

	    mass = b->get_starbase()->conv_m_dyn_to_star(1);
	    length = b->get_starbase()->conv_r_dyn_to_star(1) * Rsun_pc;
	    time = b->get_starbase()->conv_t_dyn_to_star(1);

	    return (mass > 0 && length > 0 && time > 0);
	}

    return false;
}

void add_power_law(dyn *b,
		   real coeff, real exponent, real scale,
		   vec center,		// default = (0, 0, 0)
		   bool n_flag,			// default = false
		   bool verbose,		// default = false
		   bool G_flag)			// default = false
{
    // Add power-law parameters for an external field.

    char id[10];
    int ind;

    if (exponent == 0) {
	strcpy(id, "plummer");
	ind = 14;
    } else {
	strcpy(id, "power_law");
	ind = 16;
    }

    // Check for physical parameters.

    real mass, length, time;
    bool phys = get_physical_scales(b, mass, length, time);

    // If phys, we are using physical units.  Mass and length are
    // the physical equivalents of 1 N-body unit, in Msun and pc.
    // Time is not used here.

    if (phys && !n_flag) {

	cerr << "add_" << id
	     << ":  converting input physical parameters to N-body units"
	     << endl;
	PRI(ind);
	if (exponent == 0)
	    cerr << "M";
	else
	    cerr << "A";
	cerr << " = " << coeff << " Msun,  a = " << scale
	     << " pc" << endl;
	PRI(ind); cerr << "center = (" << center << ") pc" << endl;
	PRI(ind); cerr << "N-body mass scale = " << mass
	               << " Msun,  length scale = " << length
		       << " pc" << endl;

	// Convert to N-body units.

	coeff *= pow(length, exponent) / mass;
	scale /= length;
	center /= length;

    } else {

	cerr << "add_" << id
	     << ":  interpreting input parameters as N-body units"
	     << endl;

	if (G_flag) {		// won't happen in Plummer case

	    // Warn if physical parameters are not set or are being ignored.

	    PRI(ind); cerr << "warning:  -G but ";
	    if (phys)
		cerr << "ignoring";
	    else
		cerr << "no";
	    cerr << " physical scaling" << endl;
	}
    }

    PRI(ind);
    if (phys && !n_flag) cerr << "N-body ";

    if (exponent != 0) {

	putrq(b->get_log_story(), "kira_pl_coeff", coeff);
	putrq(b->get_log_story(), "kira_pl_exponent", exponent);
	putrq(b->get_log_story(), "kira_pl_scale", scale);
	putvq(b->get_log_story(), "kira_pl_center", center);

	
	cerr << "A = " << coeff << ",  a = " << scale
	     << ",  x = " << exponent << endl;

    } else {

	putrq(b->get_log_story(), "kira_plummer_mass", coeff);
	putrq(b->get_log_story(), "kira_plummer_scale", scale);
	putvq(b->get_log_story(), "kira_plummer_center", center);
	putiq(b->get_log_story(), "kira_plummer_friction", 0);

	cerr << "M = " << coeff << ", a = " << scale << endl;
    }

    PRI(ind); cerr << "center = (" << center << ")" << endl;
}

#else

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    real coeff = 1, scale = 1, exponent = 0;
    vec center = 0;

    bool G_flag = false;
    bool n_flag = false;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    const char *param_string = "A:a:c:C:::e:E:GM:nR:x:X:";

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.15 $", _SRC_)) != -1) {
	switch (c) {
	    case 'A':
	    case 'M':	coeff = atof(poptarg);
			break;

	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;

	    case 'C':	center = vec(atof(poparr[0]),
					atof(poparr[1]),
					atof(poparr[2]));
			break;
	    case 'e':
	    case 'E':
	    case 'x':
	    case 'X':
			exponent = atof(poptarg);
			break;

	    case 'a':
	    case 'R':	scale = atof(poptarg);
			break;

	    case 'G':	G_flag = true;
	    		break;
	    case 'n':	n_flag = true;
	    		break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    dyn *b = get_dyn();
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);
    if (c_flag)	b->log_comment(comment);

    if (G_flag) {

	// Parameters from Mezger et al. (1999):

	coeff = 4.25e6;
	exponent = 1.2;

	scale = 0.1;
	center = vec(0,0,0);
    }

    add_power_law(b, coeff, exponent, scale, center, n_flag, true, G_flag);
    put_dyn(b);
}

#endif
