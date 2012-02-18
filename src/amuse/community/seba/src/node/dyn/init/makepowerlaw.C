
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct a system of test particles drawn from a power-law mass
//// distribution, in virial equilibrium with the field due to that
//// distribution.  The mass distribution is of the form:
////
////               M(r)  =  A r^x		(0 < x < 3)
////
//// except that r is modified to include a softening scale parameter:
////
////                     r  -->  sqrt(r^2 + a^2)
////
//// The power-law particle distribution is created between some minimum
//// and maximum radii.  The inner radius is for numerical convenience, to
//// limit the number of stars on very tight orbits.  In addition, we allow
//// the possibility that the power-law potential is cut off within some
//// radius b and that there a point mass at the center.  The default mass
//// for the central point is the power-law mass within the cutoff.  The
//// point-mass potential is softened to avoid numerical problems.  If b is
//// nonzero, a is set to zero.
////
//// We assume that x > 0.  Particle parameters will be chosen so that
//// the total mass of the N-body system is 1, independent of the actual
//// mass of the background distribution.  For now, we are interested only
//// in test particles.  The particles are currently always drawn from the
//// power-law distribution, even within the cutoff radius.
////
//// The output snapshot will have the external field already enabled, and
//// will contain a flag to disable internal interactions.
////
//// Usage:  makepowerlaw [OPTIONS]
////
//// Options:
////              -A    specify the coefficient A [1]
////              -a/R  specify scale [1]
////              -b    cutoff radius [0]
////              -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -d    make orbits disk-like within the cutoff [isotropic]
////              -e    softening parameter [b/100]
////              -i    number the particles sequentially [don't number]
////              -l    minimum radius of the N-body system [0.1]
////              -m/M  central mass [mass at cutoff]
////              -n    specify number of particles [no default]
////              -o    echo value of random seed [don't echo]
////              -u    maximum radius of the N-body system [10]
////              -s    specify random seed [random from system clock]
////              -x    specify exponent [1]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

local void  makepowerlaw(dyn * root, int n,
			 real A, real a, real x,
			 real cutoff, real mass,
			 real r_min, real r_max,
			 bool d_flag)
{
    real xi = 1/x, pmass = 1./n, armax = 0, vc;
    root->set_mass(1);

    if (a > 0) {
	armax = pow(a/r_max, x);
	vc =sqrt(3*A*pow(a, x-1)/(3-x));	// 3-D vrms at r = a
    } else
	vc = sqrt(3*A*pow(cutoff, x-1)/(3-x));	// 3-D vrms at r = cutoff

    for_all_daughters(dyn, root, bi) {

	bi->set_mass(pmass);

	// Radii are distributed between r_min and r_max, roughly uniformly
	// in r^x.  Power-law is OK if a << r_max; otherwise, better to use
	// a "core-halo" approximation to m(r):
	//
	//	m(r)  ~  (r/a)^3 (a/r_max)^x	(r < a)
	//		 (r/r_max)^x		(r > a)
	//
	// Use this distribution even inside the cutoff radius, if any.

	real radius = 0;
	while (radius < r_min) {	// lower limit for convenience only
	    if (randinter(0, 1) < armax)
		radius = a*pow(randinter(0, 1), 1./3);
	    else
		radius = r_max*pow(randinter(armax, 1), xi);
	}

	real costheta;

	// Details of orbits within the cutoff radius are to be determined.
	// For now, the orbits form a disk close to the central point mass.

	if (radius > cutoff || !d_flag)

	    // Spherically symmetric distribution.

	    costheta = randinter(-1, 1);

	else {

	    // Progressively force particles into the x-y plane as r --> 0.

	    real lim = sqrt(radius/cutoff);
	    costheta = randinter(-lim, lim);
	}

	real sintheta = 1 - costheta*costheta;
	if (sintheta > 0) sintheta = sqrt(sintheta);
	real phi = randinter(0.0, TWO_PI);

        bi->set_pos(radius*vec(sintheta * cos(phi),
				  sintheta * sin(phi),
				  costheta));

	// Choose velocities to ensure equilibrium in the external field.
	// Modify the velocity dispersion inside the cutoff to take the
	// central point mass into account.

	real vrms;

	if (a > 0) {

	    // Note:  scale > 0 ==> cutoff = 0.

	    vrms = vc;
	    if (radius > a) vrms *= pow(radius/a, (x-1)/2);

	} else {

	    if (radius > cutoff) 
		vrms = sqrt(3*A*pow(radius, x-1)/(3-x));
	    else
		// vrms = sqrt(mass/radius);
		vrms = Starlab::max(sqrt(mass/radius), vc);
	}

	// Play with orbital orientations.  At present, orbits become more
	// nearly planar and circular as we approach the center if d_flag
	// is set (see above).

	if (radius > cutoff || !d_flag) {

	    // Isotropic.

	    bi->set_vel(vrms*vec(randinter(-1,1),
				    randinter(-1,1),
				    randinter(-1,1)));
	} else {

	    // Some handy unit vectors.

	    vec rhat = bi->get_pos()/radius;
	    vec zhat = vec(0,0,1);

	    vec xyhat = rhat^zhat;		// transverse unit vector
	    xyhat /= abs(xyhat);		// in the x-y plane

	    real lim = sqrt(radius/cutoff);
	    real rfrac = randinter(-lim, lim);
	    real zfrac = randinter(-lim, lim);
	    real tfrac = sqrt(Starlab::max(0.0, 1-rfrac*rfrac-zfrac*zfrac));

	    bi->set_vel(vrms*(rfrac*rhat+tfrac*xyhat+zfrac*zhat));
	}
    }

    putrq(root->get_log_story(), "initial_mass", 1.0);
}

#define  SEED_STRING_LENGTH  60

int main(int argc, char ** argv) {
    int  i, n;
    int  input_seed, actual_seed;

    bool a_flag = false;
    bool c_flag = false;
    bool C_flag = false;
    bool d_flag = false;
    bool e_flag = false;
    bool i_flag = false;
    bool m_flag = false;
    bool n_flag = false;
    bool o_flag = false;
    bool v_flag = false;
    bool s_flag = false;

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    real coeff = 1, scale = 1, exponent = 1;
    vec center = 0;

    real cutoff = 0, mass = 0, softening = 0;

    real r_max = 10, r_min = 0.1;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "A:a:b:c:Cde:il:M:m:n:oR:u:s:x:v:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.13 $", _SRC_)) != -1)
	switch(c) {
	    case 'A': coeff = atof(poptarg);
		      break;
	    case 'a':
	    case 'R':
		      a_flag = true;
		      scale = atof(poptarg);
		      break;
	    case 'b': cutoff = atof(poptarg);
		      break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'd': d_flag = true;
		      break;
	    case 'e': e_flag = true;
	    	      softening = atof(poptarg);
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'l': r_min = atof(poptarg);
		      break;
	    case 'M':
	    case 'm': m_flag = true;
		      mass = atof(poptarg);
		      break;
	    case 'n': n_flag = true;
		      n = atoi(poptarg);
		      break;
	    case 'o': o_flag = true;
                      break;
	    case 'u': r_max = atof(poptarg);
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'x': exponent = atof(poptarg);
		      break;
	    case 'v': v_flag = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    if (!n_flag) {
        cerr << "makepowerlaw: must specify the number # of";
	cerr << " particles with -n#" << endl;
	exit(1);
    }
    
    if (n < 1) {
        cerr << "makepowerlaw: n > 0 required" << endl;
	exit(1);
    }

    if (exponent <= 0 || exponent >= 3) {
        cerr << "makepowerlaw: 0 < x < 3 required" << endl;
	exit(1);
    }

    if (cutoff > 0) {
	if (!m_flag) mass = coeff*pow(cutoff, exponent);
	if (!e_flag) softening = cutoff/1000;
	if (a_flag)
	    cerr << "makepowerlaw: scale = " << scale
		 << " inconsistent with cutoff = " << cutoff
		 << "; setting scale = 0" << endl;
	scale = 0;
    }

    // Create a linked list of (labeled) nodes.

    dyn *b, *by, *bo;

    b = new dyn();
    b->set_root(b);
    if (i_flag) b->set_label("root");

    if (C_flag) b->set_col_output(true);

    bo = new dyn();
    if (i_flag) bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    for (i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag) by->set_label(i+1);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    // Add comments, etc.

    if (c_flag) b->log_comment(comment);

    b->log_history(argc, argv);

    if (!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makepowerlaw: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    // Create dyn data.

    makepowerlaw(b, n, coeff, scale, exponent,
		 cutoff, mass, r_min, r_max, d_flag);

    // Flag actions for use by kira.

    if(v_flag) {
	putrq(b->get_log_story(), "kira_pl_coeff", coeff);
	putrq(b->get_log_story(), "kira_pl_exponent", exponent);
	putrq(b->get_log_story(), "kira_pl_scale", scale);
	putvq(b->get_log_story(), "kira_pl_center", center);

	putrq(b->get_log_story(), "kira_pl_cutoff", cutoff);
	putrq(b->get_log_story(), "kira_pl_mass", mass);
	putrq(b->get_log_story(), "kira_pl_softening", softening);

	putiq(b->get_log_story(), "ignore_internal", 1);
	putrq(b->get_log_story(), "r_reflect", r_max);
    }

    put_dyn(b);
    return 0;
}

#endif

/* end of: makepowerlaw.c */
