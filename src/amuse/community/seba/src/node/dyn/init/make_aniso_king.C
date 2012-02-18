
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct an anisotropic King model.
////
//// The resulting model fills its Roche lobe in the specified tidal
//// field.  Its length scale is set by G = M = 1 and the value of alpha1.
////
//// Note that, if we scale the model to standard units (default), then the
//// absolute value of alpha1 is irrelevant.  Only alpha3/alpha1 matters.
//// The default alpha3/alpha1 corresponds to a point-mass tidal field,
//// with alpha1 = -3 Omega^2 (see kira or Binney & Tremaine).
////
//// When scaling, the tidal potential is taken into account in determining
//// the virial radius.
////
//// Usage:  make_aniso_king [OPTIONS]
////
//// Options:
////              -a    specify alpha1 [-0.0009]
////              -A    specify "Oort A"
////              -b    specify alpha3/alpha1 [-1/3]
////              -B    specify "Oort B"
////              -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -F    specify tidal field type as in kira [1]
////                    1: point-mass;
////                    2: isothermal halo;
////                    3: Galactic disk;
////                    4: Custom (input Oort constants A and B,
////                    and local density in units of 232 Msun/pc^3)   
////              -i    number the particles sequentially [don't number]
////              -n    specify number of particles [no default]
////              -o    echo value of random seed [don't echo]
////              -s    specify random seed [random from system clock]
////              -T    test options [0]
////                    1: print the King model to cerr;
////                    2: save model and surface brightness profiles
////              -u    leave final N-body system unscaled
////                    [scale to E=-1/4, M = 1, R = 1]
////              -w    specify King dimensionless depth [no default]
////
//// Written by Douglas Heggie and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// version 1:  Jul 1998     Steve McMillan
//	       Jan 1999	    steve@zonker.drexel.edu
//                          Physics Dept., Drexel University, Phila PA, USA
//                          (Wrapper for D.C. Heggie's Fortran 77 program)
//             Jun 1999     Steve McMillan: added non-point-mass fields
//                          to this program and kira

#include "dyn.h"

#ifdef TOOLBOX

#define AKING_F77 F77_FUNC(aking, AKING)
#define RAN2_F77 F77_FUNC(ran2, RAN2)

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];


// Should really be 'extern "F77"'...

extern "C" void AKING_F77(real*, real*, real*, real*, real*, real*, real*,
		      real*, real*, int*, int*, real*, int*,
		      real*, real*, real*, real*);

extern "C" real RAN2_F77(int &seed)
{
    return randinter(0, 1);
}

local void make_aniso_king(dyn * b, int n, real w0, real alpha1, real alpha3,
			   bool u_flag, int test, int seed)

// Create an anisotropic King model using Heggie's Fortran code, and
// optionally initialize an N-body system with G = 1, total mass = 1,
// and length scale set by alpha1.

{
    // These tests are redundant: they duplicate the tests at the
    // start of aking.

    if (alpha3 != 0 && alpha1 >= 0)
	err_exit("make_aniso_king: must specify alpha1 < 0");
    if (w0 < 1)
	err_exit("make_aniso_king: must specify w0 > 1");
    if (w0 > 12)
	err_exit("make_aniso_king: must specify w0 < 12");

    // Compute the cluster density/velocity/potential profile

    real* m = new real[n];
    real* x = new real[n];
    real* y = new real[n];
    real* z = new real[n];
    real* vx = new real[n];
    real* vy = new real[n];
    real* vz = new real[n];

    // Enforce limit on seed required by internal random-number generator
    // in aking (ran2, from Fortran Numerical Recipes, 1ed).

    seed %= 150999;

    real potential, kinetic, tidal_energy;	// will be computed in aking;
    						// no point recomputing...

    real* coord = new real[3*n];		// workspace for aking...
						// note that this duplicates
						// the x, y, and z arrays

    AKING_F77(m, x, y, z, vx, vy, vz,			// DCH's Fortran
	  &alpha1, &alpha3, &seed, &n, &w0, &test,	// function
	  &potential, &kinetic, &tidal_energy, coord);

    if (b == NULL || n < 1) return;

    // Initialize the N-body system.

    sprintf(tmp_string,
    "       Anisotropic King model, w0 = %.2f, alpha1 = %f, alpha3 = %f",
	    w0, alpha1, alpha3);
    b->log_comment(tmp_string);

    // Assign positions and velocities.

    int i = 0;
    for_all_daughters(dyn, b, bi) {

	bi->set_mass(m[i]);
	bi->set_pos(vec(x[i], y[i], z[i]));
	bi->set_vel(vec(vx[i], vy[i], vz[i]));

	i++;
    }

    delete m, x, y, z, vx, vy, vz;

    // System is in virial equilibrium in a consistent set of units
    // with G and total mass = 1 and length scale determined by the
    // fact that the cluster exactly fills its Roche lobe in the
    // external field.

    // Transform to center-of-mass coordinates.
    // Note: must recompute the kinetic energy and the tidal potential.

    b->to_com();	

    kinetic = get_kinetic_energy(b);
    tidal_energy = get_tidal_pot(b);

    // Definition of r_virial assumes GM = 1 and *does not include* the
    // tidal potential.

    // real r_virial = -0.5/(potential + tidal_energy);
    real r_virial = -0.5/potential;

    // Fake "Jacobi radius" is used to transmit tidal (alpha1) information
    // to kira -- r_jacobi actually is the Jacobi radius only for a 1/r
    // cluster potential (actually, not a bad approximation).  Assumes GM = 1.

    real r_jacobi = pow(-alpha1, -1.0/3);

    // Optionally scale to standard parameters.  Scaling is equivalent
    // to using "scale -s" as a tool.

    if (!u_flag && n > 1) {

	// Scaling is OK because cluster is Roche-lobe filling.
	// Function scale_virial adjusts kinetic.

	// Note: scale_* operates on internal energies.

	scale_virial(b, -0.5, potential, kinetic);  // potential + tidal_energy

	real energy = kinetic + potential + tidal_energy;

	// Scale_energy uses the specified energy to rescale velocities,
	// adjusts energy accordingly, and returns the factor by which
	// the positions were scaled.

	real fac = scale_energy(b, -0.25, energy);

	// Note: when recomputing energies for test purposes, we must
	// remember to rescale alpha1,3 to preserve Roche-lobe filling.

	alpha1 /= pow(fac,3);
	alpha3 /= pow(fac,3);

	sprintf(tmp_string,
		"       Rescaled alpha1 = %f, alpha3 = %f", alpha1, alpha3);
	b->log_comment(tmp_string);

	// Recompute other quantities mainly for completeness.

	r_virial *= fac;		// should be 1
	kinetic /= fac;			// should be 0.25
	potential /= fac;		// should be -0.5
	tidal_energy /= fac;
	r_jacobi *= fac;

	// Story output mimics makeking where possible.

	putrq(b->get_log_story(), "initial_total_energy", -0.25);
	putrq(b->get_dyn_story(), "total_energy", -0.25);
    }

    // Write essential model information to the root dyn story.
    // Story output mimics makeking where possible.

    putrq(b->get_log_story(), "initial_mass", 1.0);
    putrq(b->get_log_story(), "initial_rvirial", r_virial);
    putrq(b->get_log_story(), "initial_rtidal_over_rvirial",
	  r_jacobi/r_virial, 10);

    cerr << " make_aniso_king: "; PRL(r_jacobi/r_virial);

    // Note that, for this model, kira *must* use this Jacobi
    // radius to be consistent.

    // Additional information (indicator of anisotropic King model):

    putrq(b->get_log_story(), "alpha3_over_alpha1", alpha3/alpha1, 10);
}

#define OMEGA_2	3.e-4

// Express Galactic parameters in "stellar" units (see Mihalas 1968):
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

#define OORT_ALPHA_RATIO \
	((4*M_PI*RHO_G + 2*(OORT_A*OORT_A - OORT_B*OORT_B)) \
			 / (-4*OORT_A*(OORT_A-OORT_B)))

main(int argc, char ** argv)
{
    real w0;
    int  n = 0;

    real Oort_A	= OORT_A; 	// km/s/pc
    real Oort_B = OORT_B;	// km/s/pc
    real rho_G  = RHO_G;	// (232 Msun)/pc^3

    // Default tidal parameters are appropriate to a cluster in a
    // a circular orbit of radius D = 1000 length units around a
    // point-mass "galaxy" of mass Mg = 3.0e5
    //
    //	    alpha1  = -3 G Mg / D^3  = -3 Omega^2  (Binney & Tremaine, p. 450)
    //	    alpha3  =    G Mg / D^3  =    Omega^2
	
    real alpha1 = -3 * OMEGA_2;
    real alpha3 = OMEGA_2;
    real alpha3_over_alpha1 = alpha3/alpha1;
    bool a_flag = false;
    bool b_flag = false;

    int  input_seed, actual_seed;

    bool A_flag = false;
    bool B_flag = false;
    bool G_flag = false;

    bool c_flag = false;
    bool C_flag = false;
    bool i_flag = false;
    bool n_flag = false;
    bool o_flag = false;
    bool s_flag = false; 
    bool w_flag = false;
    bool u_flag = false;

    int tidal_type = 1;
    bool F_flag = false;

    int test = 0;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "A:a:B:b:c:CF:G:in:os:T:uw:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.14 $", _SRC_)) != -1)
        switch(c) {
            case 'A': Oort_A = atof(poptarg);
	    	      A_flag = true;
                      break;
            case 'B': Oort_B = atof(poptarg);
	    	      B_flag = true;
                      break;
            case 'G': rho_G = atof(poptarg);
	    	      G_flag = true;
                      break;
            case 'a': alpha1 = atof(poptarg);
	    	      a_flag = true;
                      break;
            case 'b': alpha3_over_alpha1 = atof(poptarg);
	    	      b_flag = true;
		      F_flag = false;
                      break;
            case 'c': c_flag = true; 
                      comment = poptarg;
                      break;
	    case 'C': C_flag = true;
		      break;
            case 'F': tidal_type = atoi(poptarg);
		      b_flag = false;
	    	      F_flag = true;
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
        cerr << "make_aniso_king: please specify the dimensionless depth";
        cerr << " with -w #\n";
        exit(1);
    }

//    if (test == 0) {
//        cerr << "make_aniso_king: please specify the number # of";
//        cerr << " particles with -n #\n";
//        exit(1);
//    }

    if (n < 0)
	err_exit("make_aniso_king: n > 0 required");

    if (A_flag && B_flag && G_flag) {

      cerr << " make_aniso_king: Custom tidal field with Oort constants" 
	   << endl;
      PRI(4);PRC(Oort_A);PRC(Oort_B);PRL(rho_G);
      tidal_type = 4;
      F_flag = true;

    }

    if (F_flag) {

	if (tidal_type == 1)

	    alpha3_over_alpha1 = -1.0/3;

	else if (tidal_type == 2)

	    alpha3_over_alpha1 = -0.5;

	else if (tidal_type == 3) {

	    // An odd way to do things...

	    alpha1 = -4*OORT_A*(OORT_A-OORT_B);
	    alpha3_over_alpha1 = OORT_ALPHA_RATIO;
	    a_flag = true;

	} else if(tidal_type == 4) {
	  
	    alpha1 = -4*Oort_A*(Oort_A-Oort_B);
	    alpha3_over_alpha1 = ((4*M_PI*rho_G 
				   + 2*(Oort_A*Oort_A - Oort_B*Oort_B)) \
			       / (-4*Oort_A*(Oort_A-Oort_B)));

	    a_flag = true;
	}
	else
	    F_flag = false;

	if (F_flag) {
	    cerr << " make_aniso_king: tidal_type = " << tidal_type;
	    if (tidal_type >= 3) cerr << ", alpha1 = " << alpha1;
	    cerr << ", alpha3/alpha1 = " << alpha3_over_alpha1
		 << endl;
	    b_flag = true;
	}
    }

    if (!b_flag && !F_flag)
	cerr << " make_aniso_king: using default alpha3/alpha1 = "
	     << alpha3_over_alpha1 << endl;

    alpha3 = alpha3_over_alpha1 * alpha1;

    if (!a_flag) {
	cerr << " make_aniso_king: using default "; PR(alpha1);
	if (!b_flag) {
	    cerr << ", "; PR(alpha3);
	}
	cerr << endl;
    }

    dyn *b = NULL;

    if (!n_flag) {

	b = get_dyn();
	n = b->n_leaves();
    }
    else {

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

    if (C_flag) b->set_col_output(true);

    if (c_flag)
	b->log_comment(comment);		// add comment to story
    b->log_history(argc, argv);

    if (s_flag == false) input_seed = 0;	// default
    actual_seed = srandinter(input_seed);

    if (o_flag)
	cerr << "make_aniso_king: random seed = " << actual_seed << endl;

    sprintf(tmp_string,
	    "       random number generator seed = %d",
	    actual_seed);
    b->log_comment(tmp_string);

    make_aniso_king(b, n, w0, alpha1, alpha3, u_flag, test, actual_seed);

    b->set_tidal_field(alpha1, alpha3);
    putiq(b->get_log_story(), "kira_tidal_field_type", tidal_type);

    if(tidal_type == 4) {
	if (A_flag) putrq(b->get_log_story(), "Oort_A_constant", Oort_A, 10);
	if (B_flag) putrq(b->get_log_story(), "Oort_B_constant", Oort_B, 10);
	if (G_flag) putrq(b->get_log_story(), "local_mass_density", rho_G, 10);
    }

    if (n_flag)
	put_dyn(b);
}

#endif

