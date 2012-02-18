
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// adddouble: Adds double star class to excisting node structure.
////
//// Options:    -A    maximum initial semi-major axis [in solar units]
////             -a    minimum initial semi-major axis [in solar units]
////             -c    comment to put in the starbase log structure.
////             -E    maximum initial eccentricity
////             -e    minimum initial eccentricity
////             -M    Mass scaling [total cluster mass in solar units].
////             -Q    Virial ratio (used if time scaling is omitted)
////             -R    Dynamical size scaling for the star
////                   [in units of the virial radius].
////             -S    use internal scaling?
////             -s    Random seed.
////             -T    Dynamical time scaling
////                   [in units of the NBODY time].
////             -t    start time [in million years].
////
//++ Notes:
//++  Online arguments over rule parameters from input node.
//++
//++  If no time unit is explicitly specified it is derived
//++  from the other units.
//++
//++  The mass of the system can be set independedly from the total
//++  mass of the input N-body system (e.g. when modeling a large
//++  cluster by a smaller N-body system).
//++
//++ Example of usage:
//++  adddouble
//++
//++ See also: addstar

//// Version 1.0: Feb 1995   Simon Portegies Zwart
//                           spz@grape.c.u-tokyo.ac.jp

#include "double_star.h"
#include "util_io.h"
#include "dstar_to_dyn.h"
#include "kepler.h"

#ifndef TOOLBOX

//#define MASS_UPDATE_LIMIT 0.0001
//#define SMA_UPDATE_LIMIT 0.0001

#define REPORT_ADD_DOUBLE        false
#define REPORT_DEL_DOUBLE        false
#define REPORT_EVOLVE_DOUBLE     false

//This is a challenge, does it improve bugginess.

local void randomize_created_binary(dyn *the_binary_node)
{
     dyn *d1 = the_binary_node->get_oldest_daughter();
     dyn *d2 = d1->get_younger_sister();
     real m1 = d1->get_mass();
     real m2 = d2->get_mass();
     real m_total = m1 + m2;
     real sma = the_binary_node->get_starbase()->get_semi();
     real ecc = the_binary_node->get_starbase()->get_eccentricity();

     if (REPORT_ADD_DOUBLE) {
       cerr << "randomize_created_binary: ";
       PRC(m1);PRC(m2);PRC(sma);PRL(ecc);
     }

     if (sma<0 || ecc<0 || ecc>1) {
        cerr << "randomize_created_binary(): This is not a binary!"<<endl;
        return;
     }

     kepler johannes;

     real peri = 1; // Default value (unimportant unless ecc = 1).

     // For now, binary phase is random.

     real mean_anomaly = randinter(-PI, PI);

     make_standard_kepler(johannes, 0, m_total, -0.5 * m_total / sma, ecc,
                          peri, mean_anomaly);
     set_random_orientation(johannes);
     johannes.print_all();
}


void update_dyn_from_binary_evolution(dyn* the_binary_node,
				      dyn* d1, real dm1_fast,
				               real dm2_fast)
{

    dyn *d2 = NULL;
    if (d1==NULL) {
	d1 = the_binary_node->get_oldest_daughter();
	d2 = d1->get_younger_sister();
    }
    else
	d2 = d1->get_binary_sister();

    real dyn_m1 = get_total_mass(d1) - dm1_fast;
    real dyn_m2 = get_total_mass(d2) - dm2_fast;
    real dyn_mtot = dyn_m1 + dyn_m2;

    // kepler and binary exists.
    // As there is a kepler it should have been initialized.

    //		Needed to think about more carefully.
    kepler *johannes = the_binary_node->get_oldest_daughter()->get_kepler();
    d1->set_pos(-dyn_m2 * johannes->get_rel_pos() / dyn_mtot);
    d1->set_vel(-dyn_m2 * johannes->get_rel_vel() / dyn_mtot);

    d2->set_pos(dyn_m1 * johannes->get_rel_pos() / dyn_mtot);
    d2->set_vel(dyn_m1 * johannes->get_rel_vel() / dyn_mtot);
    d1->set_mass(dyn_m1);
    d2->set_mass(dyn_m2);


}


/*-----------------------------------------------------------------------------
 *  adddouble  -- for all particles, add a star part using "new double()".
 *
 *  Adds a double_star to a parent node of two leaves, i.e.; the
 *  double_star is attached to the center-of-mass node of the
 *  two stars.
 *  Note that the kepler pointer is attached to the leave nodes.
 *  The function checks (too extensively) wheter or not a double_star
 *  is really needed. The function is fully self-supporting:
 *  call adddouble from any root-node and it finds itself what it needs,
 *  (unlike delete_double see dstar_to_kira.C).
 *-----------------------------------------------------------------------------
 */

local void do_the_adding(hdyn *bi, real stellar_time,	// bi is binary CM
			 binary_type type,
			 real sma, real ecc)
{
  hdyn* od = (hdyn*)bi->get_oldest_daughter();

    if (stellar_time <= 0) {
	//cerr << "adding at Time=zero"<<endl;
	// Should be the individual binaries time!!
	stellar_time = bi->get_starbase()
	    		 ->conv_t_dyn_to_star(od->get_time());
				      // 3/99 was get_system_time()
    }

    if (REPORT_ADD_DOUBLE) {
	int p = cerr.precision(HIGH_PRECISION);
	cerr << "Before adddouble: at " << stellar_time << " to: ";
	bi->pretty_print_tree(cerr);
	cerr << " a=" << sma << " e=" << ecc << endl;
	od->get_kepler()->print_all(cerr);
	cerr.precision(p);
    }

    // addstar(bi, stellar_time); // Should not be needed.

    int id;
    if((od->is_low_level_node() &&
	od->get_younger_sister()->is_low_level_node()) &&
       (od->get_elder_sister() == NULL) &&
       bi->n_leaves()==2 ) {

	if (!has_dstar(od)) {
	    //bi->get_parent()->get_starbase()->get_element_type()!=Double) {

	    //    if (bi->is_parent() && !bi->is_root())
	    story * old_story = bi->get_starbase()->get_star_story();
	    bi->get_starbase()->set_star_story(NULL);

	    id = bi->get_index();

	    // cerr << "Adding binary to "<< id << " at time = "
	    //      << stellar_time << endl;

	    double_star* new_double
		= new_double_star(bi, sma, ecc, stellar_time, id, type);
	    // synchronize_binary_components(dynamic_cast(double_star*,
	    // 						  bi->get_starbase()));

	    if (REPORT_ADD_DOUBLE) {
		cerr<<"Adddouble: id, t, a, e: "
		    <<new_double->get_identity()<<" "
		    <<new_double->get_current_time() << " "
		    <<new_double->get_semi()<<" "
		    <<new_double->get_eccentricity()<<endl;
		put_state(make_state(new_double));
	    }

	    // Give the new binary the old star_story.

	    new_double->set_star_story(old_story);

	    if (REPORT_ADD_DOUBLE) {
		int p = cerr.precision(HIGH_PRECISION);
		cerr<<"After adddouble: at "<<stellar_time<< " to: ";
		new_double->dump(cerr);
		od->get_kepler()->print_all(cerr);
		put_state(make_state(new_double));
		cerr << endl;
		cerr.precision(p);
	    }
	}
	else {
	    if (REPORT_ADD_DOUBLE)
		cerr << "No double_star to node needed in adddouble"<<endl;
	}
    }
    else {
	if (REPORT_ADD_DOUBLE)
	    cerr << "Sorry, no binary node in adddouble"<<endl;
    }
}

void adddouble(hdyn *b,					// b is binary CM
	       real dyn_time, binary_type type,
	       bool random_initialization,
	       real a_min, real a_max,
	       real e_min, real e_max)
{
    if (REPORT_ADD_DOUBLE) {
	int p = cerr.precision(HIGH_PRECISION);
	cerr<<"adddouble: "<<b<<" at T="<<dyn_time<<endl;
	cerr.precision(p);
    }

    real stellar_time = b->get_starbase()->conv_t_dyn_to_star(dyn_time);
    addstar(b, stellar_time);

    real ecc, sma;
    real binary_age=0;
//      binary_type local_type = Unknown_Binary_Type;
    real a_const = log(a_max) - log(a_min);
    int id;
//      kepler johannes;

    for_all_nodes(dyn, b, bi) {
	if (bi->is_parent() && !bi->is_root()) {

	    // Note that b's star_story is DELETED here.
	    // Is this intentional?

            story * old_story = b->get_starbase()->get_star_story();
//            extract_story_chapter(local_type, sma, ecc, *old_story);
            b->get_starbase()->set_star_story(NULL);
//            delete old_story;
	    
	    binary_age = Starlab::max(bi->get_oldest_daughter()->get_starbase()
			       ->get_current_time(),
			     bi->get_oldest_daughter()->get_binary_sister()
			       ->get_starbase()->get_current_time());
	    // cerr << "times: " << binary_age<<" "<<stellar_time<<endl;
            id = bi->get_index();

            if (random_initialization) {

               sma = a_min*exp(randinter(0., a_const));
               do {
                  ecc = sqrt(randinter(0., 1.));
               }
               while (ecc<e_min && ecc>=e_max);

            } else { //if (local_type==Unknown_Binary_Type)

	      // cerr << "Non-random creation"<<endl;
//               type = local_type;
//            else {

	        if (!bi->get_oldest_daughter()->get_kepler())
		    new_child_kepler(bi, dyn_time, MAX_CIRC_ECC);
		else {
		    bi->get_oldest_daughter()->get_kepler()->	  // Probably
			set_circular_binary_limit(MAX_CIRC_ECC);  // unnecessary
		    dyn_to_child_kepler(bi, dyn_time);
		}
		
		sma = b->get_starbase()->
		    		conv_r_dyn_to_star(bi->get_oldest_daughter()->
				get_kepler()->get_semi_major_axis());
		ecc = bi->get_oldest_daughter()->get_kepler()->
		    		get_eccentricity();
            }

	    do_the_adding((hdyn*)bi, binary_age, type, sma, ecc);

            if (random_initialization) {
		randomize_created_binary(bi);
		//            update_dynamical_part_of_binary(bi);
		// The above function no longer exist
		// I hope the following does the expected
		// thing.... (JM)
		cerr << "try update dyn here\n";
		update_dyn_from_binary_evolution(bi);
		cerr << "return update dyn here\n";
	    }
	}
    }
}

#else

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
    {
    int  c;
    bool  A_flag = FALSE;
    bool  a_flag = FALSE;
    bool  E_flag = FALSE;
    bool  e_flag = FALSE;
    bool  t_flag = FALSE;
    bool  M_flag = FALSE;
    bool  R_flag = FALSE;
    bool  Q_flag = FALSE;
    bool  T_flag = FALSE;
    bool  p_flag = FALSE;
    bool  S_flag = FALSE;
    bool  c_flag = FALSE;
    bool  random_initialization = false;
    real  m_tot = 1;
    real  r_vir = 1;
    real  t_hc = 1;
    real  Q_vir = 0.5;
    real  a_min = 1;
    real  a_max = 1.e+6;
    real  e_min = 0;
    real  e_max = 1;
    real  t_start = 0;
    binary_type type = Detached;

    int random_seed = 0;
    char seedlog[64];
    char  *comment;
    extern char *poptarg;

    const char * param_string = "A:a:E:e:M:R:Q:T:t:Ss:c:";
    check_help();

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c)
	    {
            case 'A': A_flag = true;
		      a_max = atof(poptarg);
                      break;
            case 'a': a_flag = true;
                      a_min = atof(poptarg);
                      break;
            case 'E': E_flag = true;
                      e_max = atof(poptarg);
                      break;
            case 'e': e_flag = true;
                      e_min = atof(poptarg);
                      break;
            case 'M': M_flag = true;
		      m_tot = atof(poptarg);
                      break;
            case 'R': R_flag = true;
	              r_vir = atof(poptarg);
                      break;
            case 'Q': Q_flag = true;
	              Q_vir = atof(poptarg);
	              break;
            case 'T': T_flag = true;
		      t_hc = atof(poptarg);
                      break;
	    case 't': t_start = atof(poptarg);
		      break;
            case 's': random_seed = atoi(poptarg);
                      break;
            case 'S': S_flag = true;
                      break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	    }

    int actual_seed = srandinter(random_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

       dyn* b = get_dyn();

       b->log_history(argc, argv);
       b->log_comment(seedlog);

       if (!T_flag) t_hc = 21.0 * sqrt(Q_vir/m_tot) * pow(r_vir, 1.5);

       if (S_flag)
          if (M_flag && R_flag && T_flag)
             b->get_starbase() ->set_stellar_evolution_scaling(m_tot, r_vir,
							       t_hc);
	   else {
	       cerr << "Proper scaling not set." << endl; cerr
	            << "Use the -M -R an -T options to set the scaling for "
		    << "the mass, \nthe size and the time" << endl; exit(1);
	   }

          if (A_flag || a_flag || E_flag || e_flag)
           random_initialization=true;
       adddouble(b, t_start, type, random_initialization,
	         a_min, a_max, e_min, e_max);
       put_dyn(b);
       delete b;
    }

#endif

/* endof: adddouble.c */
