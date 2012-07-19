
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Complete the binary formation process by adding binary orbits to
//// an existing binary tree.  It is assumed that the binary tree
//// structure and all masses have already been set (e.g. by
//// makesecondary), and that positions and velocities of all top-level
//// nodes are already known.
////
//// If possible, use the system parameters to scale the binary parameters.
//// If the total system energy is already known (saved in the snapshot
//// dyn story), then energies are in units of kT.  Otherwise, energies
//// are in absolute units.
////
//// Usage:  makeplanetary [OPTIONS]
////
//// Options:
////            -f    function select option [3]
////
////                  1: angular momentum per unit reduced mass
////                  (L^2 = am[1-e^2]), solar units;
////
////                  2: semi-major axis or peri/apo, solar units;
////
////                  3: energy
////
////            -l    lower limit on selected binary parameter [1]
////
////            -o    specify interpretation of limits [1]
////
////                  (-f 1)
////                  1: angular momentum,
////                  2: angular momentum, detached binary
////
////                  (-f 2)
////                  1: semi-major axis,
////                  2: semi-major axis, detached,
////                  3: -l gives peri, -u gives apo, and detached
////
////                  (-f 3)
////                  1: |binary energy|,
////                  2: |binary energy| per unit reduced mass,
////                  3: |binary energy| per unit binary mass
////
////            -s    specify random seed [take from system clock]
////            -u    upper limit on selected binary parameter [1]
////
//// Written by Steve McMillan and Simon Portegies Zwart..
////
//// Report bugs to starlab@sns.ias.edu.

//   Steve McMillan, July 1996
//   Modified by SPZ 21 April 1998

#include "dyn.h"

#ifdef TOOLBOX

local void add_dynamics(dyn* cm, real ecc, real energy)
{
    dyn* primary = cm->get_oldest_daughter();
    dyn* secondary = primary->get_younger_sister();

    real m_total = cm->get_mass();
    real m1 = primary->get_mass();
    real m2 = secondary->get_mass();

    // Set internal orbital elements:

    kepler k;

    real peri = 1; // Default value (unimportant unless ecc = 1).
    if (ecc == 1) peri = 0;

    // For now, binary phase is random.

    real mean_anomaly = randinter(-PI, PI);

    // Energies here are really binding energies (i.e. > 0 for a bound
    // orbit) per unit mass; kepler package expects energy < 0.

    energy = -energy;

    make_standard_kepler(k, 0, m_total, energy, ecc,
			 peri, mean_anomaly);
    //PRC(m_total);
    //PRC(energy);
    //PRL(ecc);

    set_random_orientation(k);

    // Set positions and velocities.

    primary->set_pos(-m2 * k.get_rel_pos() / m_total);
    primary->set_vel(-m2 * k.get_rel_vel() / m_total);

    secondary->set_pos(m1 * k.get_rel_pos() / m_total);
    secondary->set_vel(m1 * k.get_rel_vel() / m_total);
}

local real minimum_semi_major_axis(dyn* b1, dyn* b2)
{
    real ms_prim = b1->get_starbase()->conv_m_dyn_to_star(b1->get_mass());
    real ms_sec  = b2->get_starbase()->conv_m_dyn_to_star( b2->get_mass());
    real rs_prim = b1->get_starbase()->conv_r_star_to_dyn(ms_prim);
    real rs_sec  = b2->get_starbase()->conv_r_star_to_dyn(ms_sec);
    real ms_tot  = ms_prim + ms_sec;
  
    real sma_prim = 2.2*rs_prim*pow(ms_tot/ms_prim, ONE_THIRD);
    real sma_sec  = 2.2*rs_sec*pow(ms_tot/ms_sec, ONE_THIRD);
    // a_min is Roche-lobe limited.

    return Starlab::min(sma_prim, sma_sec);
}

local void add_secondary(node* original, real mass_ratio)
{
    node* primary = new node;
    node* secondary = new node;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    primary->set_mass(original->get_mass());
    secondary->set_mass(mass_ratio*original->get_mass());
    original->inc_mass(secondary->get_mass());

    // Naming convention:

    if (original->get_name() == NULL)
	if (original->get_index() >= 0) {
	    char tmp[64];
	    sprintf(tmp, "%d", original->get_index());
	    original->set_name(tmp);
	}

    if (original->get_name() != NULL) {

	// Label components "a" and "b".

	primary->set_name(original->get_name());
	secondary->set_name(original->get_name());
	strcat(primary->get_name(), "a");
	strcat(secondary->get_name(), "b");

	// Make standard CM name.

	char tmp[256];
	sprintf(tmp, "(%s,%s)", primary->get_name(), secondary->get_name());
	original->set_name(tmp);
    }
}

local void makebinary(dyn* b, real lower, real upper,
		      int select, int option, real emax)
{
    // Binary parameters are drawn from a 1/x distribution,
    // between the specified limits.
    // For now, use thermal eccentricities.

    real frac = upper/lower;

    for_all_daughters(dyn, b, bi)
	if(!bi->is_parent()) {
	    add_secondary(bi, 0.0001);
	}

    for_all_daughters(dyn, b, bi)
	if (bi->get_oldest_daughter()) {

	    dyn* primary = bi->get_oldest_daughter();
	    dyn* secondary = primary->get_younger_sister();
	    real m_total = bi->get_mass();
	    real mu = primary->get_mass() * secondary->get_mass() / m_total;

	    // Function select options:
	    //
	    // Choose binary parameters:  select = 1:  angular momentum
	    //                            select = 2:  semi-major axis
	    //                            select = 3:  energy
	    //
    	    // 1. Select on angular momentum:
	    //
	    //   option = 1 ==> limits refer to angular momentum
	    //   option = 2 ==> limits refer to angular momentum + detached
	    //
    	    // 2. Select on semi-major axis:
	    //
	    //   option = 1 ==> limits refer to semi-major axis
	    //   option = 2 ==> limits refer to semi-major axis + detached
	    //   option = 3 ==> limits refer to pericenter & apocenter
	    //							+ detached
    	    //
	    // 3. Select on energy:
	    //
	    //   option = 1 ==> limits refer to binary energy.
	    //   option = 2 ==> limits refer to binary energy per unit
	    //				              		reduced mass.
	    //   option = 3 ==> limits refer to binary energy per unit
	    //							binary mass.
	    //
	    // However, add_dynamics expects energy per unit reduced mass,
	    // so scale the energy in cases 1 and 3.
	    //
	    // Selection on angular momentum and semi-major axis are
	    // determined iteratively using a do-loop.
	    // Infinite looping is prohibited by an initial check.
	    // Beware, however, such looping might be expensive/dangerous.

	    real ecc, energy = 0;

	    if (select == 1) {

		real l_const, sma, angular_momentum;

		if (option == 2) {
		    bool detached = false;
		    real a_min = minimum_semi_major_axis(primary, secondary);
		    if (pow(upper, 2) <= a_min*m_total*(1-pow(ecc, 2)))
			err_exit(
			    "makebinary: Illegal limits on angular momentum");
		    do {
			ecc = sqrt(randinter(0,emax));
			l_const = log(upper) - log(lower);
			angular_momentum = lower * pow(frac, randinter(0,1));
			sma = pow(angular_momentum, 2)
			    	/ (m_total*(1-pow(ecc, 2)));
			if (sma*(1-pow(ecc, 2)) > a_min) detached = true;
		    }
		    while(!detached);
		}
		else {
		    ecc = sqrt(randinter(0,emax));
		    l_const = log(upper) - log(lower);
		    angular_momentum = lower * pow(frac, randinter(0,1));
		    sma = pow(angular_momentum, 2) / (m_total*(1-pow(ecc, 2)));
		}
		energy = 0.5 * m_total / sma;
	    }

	    else if (select == 2) {

		real semi_major_axis, a_const;

		if (option >= 2) {
		    bool detached = false;
		    real a_min = minimum_semi_major_axis(primary, secondary);
		    if(upper<=a_min)
			err_exit(
			    "makebinary: Illegal limits on angular momentum");

		    if (option == 3) {// limits between peri- and apocenter.
			real pericenter, apocenter;
			do {
			    ecc = sqrt(randinter(0,emax));
			    pericenter = lower*(1-ecc);
			    apocenter = upper*(1+ecc);
			    a_const = log(upper) - log(lower);
			    semi_major_axis = lower*exp(randinter(0., a_const));
			    if (semi_major_axis*(1-pow(ecc, 2)) > a_min  &&
				(semi_major_axis > pericenter &&
				 semi_major_axis < apocenter))
				detached = true;
			}
			while(!detached);
		    }
		    else {	// separation limited by semi-detached.

			// Assume for now that
			//	stellar radius [Rsun] = mass [Msun].
			// This is of course not correct, but for the moment
			// good enough.  makebinary does not know much about
			// stars.
			do {
			    ecc = sqrt(randinter(0,emax));
			    a_const = log(upper) - log(lower);
			    semi_major_axis = lower*exp(randinter(0., a_const));
			    if(semi_major_axis*(1-pow(ecc, 2))>a_min)
				detached = true;
			}
			while(!detached);
		    }
		}
		else {
		    ecc = sqrt(randinter(0,emax));
		    a_const = log(upper) - log(lower);
		    semi_major_axis = lower*exp(randinter(0., a_const));
		}
		energy = 0.5 * m_total / semi_major_axis;
	    }

	    else {				// default is select = 3

		ecc = sqrt(randinter(0,emax));
		energy = lower * pow(frac, randinter(0,1));
		if (option == 1)
		    energy /= mu;
		else if (option == 3)
		    energy *= m_total / mu;
	    }

	    //PR(m_total);
	    //PR(primary->get_mass());
	    //PRL(secondary->get_mass());
	    //PRL(energy);
            //real sma = 0.5 * (primary->get_mass() + secondary->get_mass())
	    //               / energy;
	    //real L = sqrt(m_total*sma*(1-pow(ecc, 2)));
	    //PR(sma); PR(ecc); PRL(L);

	    add_dynamics(bi, ecc, energy);
	    add_dynamics(bi, ecc, 2*energy);
	    add_dynamics(bi, ecc, 4*energy);
	    add_dynamics(bi, ecc, 8*energy);
	}
}

local void scale_limits(real& e1, real& e2, int option,
			int N, real M, real E)
{
    real scale = 1.0;

PRL(E);
    if (option == 3) {

	// Limits refer to kinetic energy per unit mass.  Scale to -E/M.

	scale = -E / M;

    } else if (option != 2) {

	// Limits refer to energy.  Scale to kT = -(2/3) E/N.

	scale = -E / (1.5*N);
PRL(scale);
    }

    e1 *= scale;
    e2 *= scale;
}

int main(int argc, char ** argv)
{
    real lower = 1.0, upper = 1.0;
    real emax = 1;
    int function_select = 3;
    int option = 1;
    int random_seed = 0;
    char seedlog[64];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "f:l:s:o:u:e:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {

	    case 'f': function_select = atoi(poptarg);
		      break;
	    case 'l': lower = atof(poptarg);
		      break;
	    case 'o': option= atoi(poptarg);
		      break;
	    case 'e': emax = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
	    case 'u': upper = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    if (lower <= 0 || upper <= 0 || lower > upper)
	err_exit("makebinary: Illegal limits");

    dyn* b;
    b = get_dyn();

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d",
	    actual_seed);
    b->log_comment(seedlog);

    // Look in the dyn story to see if a total system energy is known.
    // If it is, scale lower and upper accordingly.

    // Thus, in case function_select = 3, if the energy is known,
    // lower and upper are specified in "kT" units.  If the energy
    // is not known, lower and upper are absolute limits on the
    // binary kinetic energy.

    const char* energy_string = "initial_total_energy";

    if (function_select == 1) {
	if (b->get_starbase()->get_stellar_evolution_scaling()) {
	    real scale = sqrt(b->get_starbase()->conv_m_star_to_dyn(1)
			      * b->get_starbase()->conv_r_star_to_dyn(1));
	    lower *= scale;
	    upper *= scale;
	}
	else
	    cerr << "makebinary: Using unscaled angular-momentum limits.\n";
    }
    else if (function_select == 2) {
	if (b->get_starbase()->conv_r_star_to_dyn(1)>0) {
	    lower = b->get_starbase()->conv_r_star_to_dyn(lower);
	    upper = b->get_starbase()->conv_r_star_to_dyn(upper);
	}
	else
	    cerr << "makebinary: Using unscaled semi-major axis limits.\n";
    }
    else if (function_select == 3) {
	if (find_qmatch(b->get_log_story(), energy_string))
	    scale_limits(lower, upper, option,
			 b->n_daughters(), b->get_mass(),
			 getrq(b->get_log_story(), energy_string));
	else
	    cerr << "makebinary: Using unscaled energy limits.\n";
    }

    makebinary(b, lower, upper, function_select, option, emax);

    put_dyn(b);
    rmtree(b);
    return 0;
}
#endif
