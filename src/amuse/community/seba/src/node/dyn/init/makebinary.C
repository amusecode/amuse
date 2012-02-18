
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Complete the binary formation process by adding binary orbits to an
//// existing binary tree.  It is assumed that the binary tree structure
//// and all masses have already been set (e.g. by makesecondary), and
//// that positions and velocities of all top-level nodes are already known.
////
//// If possible, use the system parameters to scale the binary parameters.
//// If the total system energy is already known (saved in the snapshot
//// log or dyn story), then energies are in units of kT.  Otherwise, energies
//// are in absolute units.
////
//// Usage:  makebinary [OPTIONS]
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
////            -e    maximum eccentricity [1]
////            -l    lower limit on selected binary parameter [1]
////
////            -o    specify interpretation of limits [1]
////
////                  (-f 1)
////                  1: angular momentum,
////                  2: angular momentum, detached binary
////
////                  (-f 2) [unit: Rsun]
////                  1: semi-major axis,
////                  2: semi-major axis, detached and hard
////                  (note: -u gives an upper limit on sma),
////                  3: -l gives peri, -u gives apo, and detached
////
////                  (-f 3) [N-body units]
////                  1: |binary energy|,
////                  2: |binary energy| per unit reduced mass,
////                  3: |binary energy| per unit binary mass
////
////            -s    specify random seed [take from system clock]
////            -u    upper limit on selected binary parameter [1]
////
//// Written by Steve McMillan and Simon Portegies Zwart.
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

real zero_age_main_sequnece_radius(const real mass) {

    real alpha, beta, gamma, delta, kappa, lambda;

    real log_mass = log10(mass);

    if (mass > 1.334) {

	alpha = 0.1509 + 0.1709*log_mass;
	beta  = 0.06656 - 0.4805*log_mass;
	gamma = 0.007395 + 0.5083*log_mass;
	delta = (0.7388*pow(mass, 1.679) - 1.968*pow(mass, 2.887))
    	      / (1.0 - 1.821*pow(mass, 2.337));
    } 
    else {
      
	alpha = 0.08353 + 0.0565*log_mass;
	beta  = 0.01291 + 0.2226*log_mass;
	gamma = 0.1151 + 0.06267*log_mass;
	delta = pow(mass, 1.25) * (0.1148 + 0.8604*mass*mass)
              / (0.04651 + mass*mass);
    }

    real radius = delta;

    if (radius < 0.1) {
 
	// assumed brown dwarf or planet.

	radius = 0.1;
    }

    return radius;
}

local real roche_radius(const real m1, const real m2) {

  real q = m1/m2;
  real q1_3 = pow(q, ONE_THIRD);
  real q2_3 = pow(q1_3, 2);   
  
  return 0.49*q2_3/(0.6*q2_3 + log(1 + q1_3));
}


local real minimum_semi_major_axis(dyn* b1, dyn* b2)
{
    real ms_prim = b1->get_starbase()->conv_m_dyn_to_star(b1->get_mass());
    real ms_sec  = b2->get_starbase()->conv_m_dyn_to_star(b2->get_mass());

    // Use stellar mass as radius indicater.
    // makebinary known little about stars
    real rs_prim = zero_age_main_sequnece_radius(ms_prim);
    real rs_sec = zero_age_main_sequnece_radius(ms_sec);

    // real rs_prim = b1->get_starbase()->conv_r_star_to_dyn(ms_prim);
    // real rs_sec  = b2->get_starbase()->conv_r_star_to_dyn(ms_sec);
  
    real sma_prim = rs_prim/roche_radius(ms_prim, ms_sec);
    real sma_sec = rs_sec/roche_radius(ms_sec, ms_prim);

    return Starlab::max(sma_prim, sma_sec);
}

local bool dyn_present(dyn* bi) {

  bool dynamics_present = true;

  if(bi->get_pos()[0]==0 && bi->get_pos()[1]==0 && bi->get_pos()[2]==0 &&
     bi->get_vel()[0]==0 && bi->get_vel()[1]==0 && bi->get_vel()[2]==0)
    dynamics_present = false;

  return dynamics_present;
}

local void makebinary(dyn* b, real lower, real upper,
		      int select, int option, real emax)
{
    // Binary parameters are drawn from a 1/x distribution,
    // between the specified limits.
    // For now, use thermal eccentricities.

    for_all_daughters(dyn, b, bi)
	if (bi->get_oldest_daughter()) {

	    dyn* primary = bi->get_oldest_daughter();
	    dyn* secondary = primary->get_younger_sister();
	    vec nul = 0;
	    if(!dyn_present(primary) || !dyn_present(secondary)) {

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
		    a_min = b->get_starbase()->conv_r_star_to_dyn(a_min);

		    if (pow(upper, 2) <= a_min*m_total*(1-pow(ecc, 2))) {
		      PRC(upper);PRL(a_min*m_total*(1-pow(ecc, 2)));
		      err_exit("makebinary: Illegal limits on angular momentum");
		    }
		    a_min = Starlab::max(lower, a_min);
		    real frac = upper/lower;

		    do {
			ecc = sqrt(randinter(0,emax));
			l_const = log(upper) - log(a_min);
			angular_momentum = a_min
			                 * pow(frac, randinter(0,1));
			sma = pow(angular_momentum, 2)
			    	/ (m_total*(1-pow(ecc, 2)));
			if (sma*(1-ecc) > a_min) 
			  detached = true;
		    }
		    while(!detached);
		}
		else {
		    real frac = upper/lower;
		    ecc = sqrt(randinter(0,emax));
		    l_const = log(upper) - log(lower);
		    angular_momentum = lower*pow(frac, randinter(0,1));
		    sma = pow(angular_momentum, 2) / (m_total*(1-pow(ecc, 2)));
		}
		energy = 0.5 * m_total / sma;
	    }

	    else if (select == 2) {

		real semi_major_axis, a_const;

		if (option >= 2) {

		    bool detached = false;
		    real a_min = 5*minimum_semi_major_axis(primary, secondary);
		    a_min = b->get_starbase()->conv_r_star_to_dyn(a_min);

		    //binding energy = 0.5 * m_total / sma;
		    real min_binding_energy = 10; // [kT] -- arbitrary!
		    real a_max = 0.5 * m_total/min_binding_energy;
		    if(upper>a_min) {
		      a_max = Starlab::min(a_max, upper);
		    }

		    if(a_max<=a_min) {
		      PRC(a_max);PRL(a_min);
		      err_exit("makebinary: Illegal limits on semi major axis");
		    }
		    a_min = Starlab::max(lower, a_min);

		    if (option == 3) {// limits between peri- and apocenter.
			real pericenter, apocenter;
			do {
			    ecc = sqrt(randinter(0, emax));
			    pericenter = a_min*(1-ecc);
			    apocenter = a_max*(1+ecc);
			    a_const = log(a_max) - log(a_min);
			    semi_major_axis = a_min
			                    * exp(randinter(0., a_const));
			    if (semi_major_axis*(1-ecc) > a_min  &&
				(semi_major_axis > pericenter &&
				 semi_major_axis < apocenter))
				detached = true;
			}
			while(!detached);
//	PRC(b->get_starbase()->conv_r_dyn_to_star(a_min));
//      PRC(b->get_starbase()->conv_r_dyn_to_star(semi_major_axis));PRL(ecc);
		    }
		    else {	// separation limited by semi-detached.

			do {
			    ecc = sqrt(randinter(0, emax));
			    a_const = log(a_max) - log(a_min);
			    semi_major_axis = a_min
			                    * exp(randinter(0., a_const));
			    // if(semi_major_axis*(1-ecc)>a_min)
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
		    semi_major_axis = b->get_starbase()
                                       ->conv_r_star_to_dyn(semi_major_axis);


		}
		energy = 0.5 * m_total / semi_major_axis;
	    }

	    else {				// default is select = 3

	        real frac = upper/lower;
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
	}
	}
}

local void scale_limits(real& e1, real& e2, int option,
			int N, real M, real K)
{
    real scale = 1.0;

    // PRL(K);

    if (option == 3) {

	// Limits refer to kinetic energy per unit mass.  Scale to K/M.

	scale = K / M;

    } else if (option != 2) {

	// Limits refer to energy.  Scale to kT = (2/3) K/N.

	scale = K / (1.5*N);
	// PRL(scale);
    }

    e1 *= scale;
    e2 *= scale;
}

int main(int argc, char ** argv)
{
  real lower = 1.0, upper = 1.0;   //Units depends on -o and -f
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
		    "$Revision: 1.26 $", _SRC_)) != -1)
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

    if (lower < 0 || upper < 0 || lower > upper)
	err_exit("makebinary: Illegal limits");

    dyn* b;
    b = get_dyn();

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d",
	    actual_seed);
    b->log_comment(seedlog);

    // In case function_select = 3, if the energy is known, then
    // lower and upper are specified in "kT" units.  If the energy
    // is not known, lower and upper are absolute limits on the
    // binary kinetic energy.

    // Look in the log or dyn story to see if a total system energy
    // is known.  If it is, scale lower and upper accordingly.

    const char* log_energy_string = "initial_total_energy";
    const char* dyn_energy_string = "total_energy";

    // Scale input to N-body units where appropriate.
    // Possible options are:
    //   1: Angular momentum   [cgs]
    //   2: Semi-major axis    [Rsun]
    //   3: Binding energy     [N-body]

    cerr << "makebinary: "; PRL(function_select);
    if (function_select == 1) {

	if (b->get_starbase()->get_stellar_evolution_scaling()) {
	    real scale = sqrt(b->get_starbase()->conv_m_star_to_dyn(1)
			      * b->get_starbase()->conv_r_star_to_dyn(1));
	    lower *= scale;
	    upper *= scale;
	}
	else
	    cerr << "makebinary: using unscaled angular-momentum limits.\n";

    } else if (function_select == 2) {

	if (b->get_starbase()->conv_r_star_to_dyn(1)>0) {
	  //	    lower = b->get_starbase()->conv_r_star_to_dyn(lower);
	  //	    upper = b->get_starbase()->conv_r_star_to_dyn(upper);
	    cerr << "makebinary: using unscaled semi-major axis limits"
		 << endl;
	}
	else
	    cerr << "makebinary: using unscaled semi-major axis limits"
		 << endl;

    } else if (function_select == 3) {

	if (find_qmatch(b->get_dyn_story(), dyn_energy_string)
	    || find_qmatch(b->get_log_story(), log_energy_string)) {

	    // Use energy_string as an indicator that some energies
	    // are known, but don't use its value as an estimator
	    // of the kinetic energy, as it may include an external
	    // field and also includes any CM motion of the system.

	    // Recompute the top-level kinetic energy in the CM frame.

	    vec com_pos, com_vel;
	    compute_com(b, com_pos, com_vel);
	    com_vel -= b->get_pos();		// energies are relative
						// to the parent node

	    real kin = get_top_level_kinetic_energy(b)
			    - 0.5*b->get_mass()*square(com_vel);

	    // PRL(get_top_level_kinetic_energy(b));
	    // PRL(0.5*b->get_mass()*square(com_vel));

	    scale_limits(lower, upper, option,
			 b->n_daughters(), b->get_mass(), kin);

	} else
	    cerr << "makebinary: Using unscaled energy limits.\n";

    }

    makebinary(b, lower, upper, function_select, option, emax);

    put_dyn(b);
    rmtree(b);
    return 0;
}
#endif
