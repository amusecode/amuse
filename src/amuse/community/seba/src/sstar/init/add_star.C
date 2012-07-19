
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// addstar:  Add physical stellar data to existing node structure.
////
//// Options:    -c    comment to put in the starbase log structure [none]
////             -M    mass scaling (system mass unit, in solar units)
////             -Q    virial ratio (used if time scaling is omitted) [0.5]
////             -R    dynamical size scaling for stars
////                       (system length unit, in parsecs)
////             -s    initial types of stars [main_sequence]
////             -T    dynamical time scaling  (system time unit, in Myr)
////             -t    stellar age (in millions of years) [0]
////
//++ Notes:
//++  Command-line arguments overrule parameters from input snapshot.
//++
//++  If no time unit is explicitly specified it is derived
//++  from the other units.
//++
//++  The mass of the system can be set independently from the total
//++  mass of the input N-body system (e.g. when modeling a large
//++  cluster by a smaller N-body system).
//++
//++  Should run makemass first to establish stellar masses.
//++
//++ Example of usage:
//++  makenode -n 10 | makemass -u 100 -l 10 | add_star -T 1
//++
//++ See also:  mkmass
//++            mknode
//++
//// Version 1.0: Apr 1993   Piet Hut, Steve McMillan, Jun Makino
//// Version 2.0: May 1993   Simon Portegies Zwart
//                           spz@grape.c.u-tokyo.ac.jp

#include "single_star.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  addstar  -- for all particles, add a star part using "new star()".
//-----------------------------------------------------------------------------

void  addstar(node * b, real t_current, stellar_type type, real z,int id, bool verbose)
{
  //    if(!((star*)b->get_starbase())->get_seba_counters()) {
  //      cerr << "Initialize SeBa counters" << endl;
  //      ((star*)b->get_starbase())->set_seba_counters(new seba_counters);
  //    }
    if (b->get_oldest_daughter() != NULL) {
	real m_min = VERY_LARGE_NUMBER;
	real m_max = -VERY_LARGE_NUMBER;
	real m_av = 0;
	int  n_star = 0;

	for_all_daughters(node, b, bi) {
	    addstar(bi, t_current, type, z, id, verbose);

	    if (verbose) {
		real m = bi->get_starbase()
			   ->conv_m_dyn_to_star(bi->get_mass());
		m_min = Starlab::min(m_min, m);
		m_max = Starlab::max(m_max, m);
		m_av += m;
		n_star++;		
	    }
	}
	if (verbose) {
	    if (n_star > 0) m_av /= n_star;
	    //  PRC(n_star); PRC(m_min); PRC(m_max); PRL(m_av);
	}
    } else {

//	cerr << "addstar: ";
//	PRL(b->get_starbase()->get_element_type());

	// NAS is the default return value of get_element_type, in
	// cases where only a starbase (but no star) exists.

	if (b->get_starbase()->get_element_type() != NAS) return;
	
	//int id = b->get_index();
        real t_cur=0, t_rel=0, m_rel=1, m_env=0, m_core=0.0;//0.01;
	real T_eff, L_eff;
	real p_rot=0, b_fld=0;
	real m_tot;

	// Create a (single) star part, using the unformation obtained from
	// the star story read in by get_node(), or created from scratch.

	stellar_type local_type = NAS;
	starbase * old_starbase = b->get_starbase();
	story * s = old_starbase->get_star_story();

	real mco_core = 0;
	extract_story_chapter(local_type, z, t_cur, t_rel,
			      m_rel, m_env, m_core, mco_core, T_eff, L_eff,
			      p_rot, b_fld, *s);

	m_tot = m_env+m_core;
	old_starbase->set_star_story(NULL);
	delete s;

	if (local_type != NAS) {	// Is the star properly defined?

	    // Use the data extracted from the star story.  Note that the
	    // input values of the type and t_rel arguments are IGNORED.

	  //cerr << "addstar: " << type_string(local_type)
	  //   << " t_rel = " << t_rel << " t_cur = " << t_cur << endl;
	
	    type = local_type;
	    single_star* new_star = new_single_star(type, id, z, t_cur,
						    t_rel,
						    m_rel, m_tot, m_core,
						    mco_core,
						    p_rot, b_fld, b);

	    //new_star->set_current_time(t_current); // added (SPZ:2/1998)
	    //new_star->dump(cerr);
	
	} else {		// No star story present, or at least no
	     			// proper definition for a single star.

	    // Create a default star of the specified type.

	    real t_cur=0, t_rel = 0, m_rel = 1, m_core = 0.0;//0.01;
	    real p_rot=0, b_fld=0;
	    real m_tot;
	    t_cur = t_current;
	//	real z=cnsts.parameters(Zsun);
	    //cerr << "addstar2: " << type_string(type) << " t_rel = " << t_rel
	    // << " t_cur = " << t_cur << endl;

	    starbase * old_starbase = b->get_starbase();
	    //id = b->get_index();
	    m_rel = m_tot = b->get_starbase()
			     ->conv_m_dyn_to_star(b->get_mass());

	    stellar_type local_type = type;

	    // Treat by dynamics pre-requested black holes
	    if(getiq(b->get_log_story(), "black_hole")==1) {
	      local_type = Black_Hole;
	      m_core = m_tot;
	      mco_core = m_tot;
	    }
	    else if(m_tot<cnsts.parameters(minimum_main_sequence)) {
	      local_type = Brown_Dwarf;
	      m_core = 0.01*m_tot;
	      mco_core = 0;
	    }
	
	    single_star* new_star = new_single_star(local_type, id, z,
						    t_cur, t_rel,
						    m_rel, m_tot, m_core,
						    mco_core,
						    p_rot, b_fld, b);

	}

	// Should not be needed since new_single_star does the job.
	// (SPZ:25 Nov 1998)
        // delete old_starbase;
    }
}

#else

main(int argc, char ** argv)
{
    int  c;
    bool  t_flag = FALSE;
    bool  M_flag = FALSE;
    bool  R_flag = FALSE;
    bool  T_flag = FALSE;
    bool  Q_flag = FALSE;
    bool  s_flag = FALSE;
    bool  c_flag = FALSE;
    real  m_tot = -1;
    real  r_vir = -1;		// NB code length unit; may or may not
				// actually be the virial radius
    real  t_vir = -1;
    real  t_rel = 0;
    real  q_vir = 0.5;
    real  T_start = 0;
    stellar_type type = Main_Sequence;
    char * star_type_string;

    char  *comment;
    extern char *poptarg;
    const char * param_string = "M:R:Q:T:t:s:c:";

    check_help();

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'M': M_flag = true;
		      m_tot = atof(poptarg);
		      break;
            case 'R': R_flag = TRUE;
                      r_vir = atof(poptarg);
                      break;
	    case 'Q': Q_flag = TRUE;
	              q_vir = atof(poptarg);
	              break;
            case 'T': T_flag = TRUE;
                      t_vir = atof(poptarg);
                      break;
	    case 't': T_flag = TRUE;
		      t_rel = atof(poptarg);
		      break;
            case 's': s_flag = true;
		      star_type_string = poptarg;
                      type = extract_stellar_type_string(star_type_string);
//            case 's': s_flag = TRUE;
//                      type = (stellar_type)atoi(poptarg);
//                      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	    }

    node *b;
    b = get_node();

    if (c_flag == TRUE)
      b->log_comment(comment);

    b->log_history(argc, argv);

    // See if any scaling parameters are specified in the input snapshot.

    real new_r_vir=-1, new_t_vir=-1, new_m_tot=-1;

    // System mass in solar units:
    real old_m_tot = b->get_starbase()->conv_m_dyn_to_star(1);

#define Rsun_pc 2.255e-8	// R_sun/1 parsec = 6.960e+10/3.086e+18;

    // Code length unit in parsecs:
    real old_r_vir = b->get_starbase()->conv_r_dyn_to_star(1) * Rsun_pc;

    // Code time unit in Myr:
    real old_t_vir = b->get_starbase()->conv_t_dyn_to_star(1);

    // Set new parameters from old ones or command line.

    if (old_m_tot > 0 && !M_flag)
	new_m_tot = old_m_tot;
    else if (M_flag)
	new_m_tot = m_tot;
    else
	err_exit("addstar: No mass scaling available.");

    if (old_r_vir > 0 && !R_flag)
	new_r_vir = old_r_vir;
    else
	new_r_vir = r_vir;

    if (old_t_vir > 0 && !(T_flag || Q_flag))
	new_t_vir = old_t_vir;
    else
	new_t_vir = t_vir;

    //    if(T_flag) {
    //      T_start = t_rel;
    //    }

    // Try to derive quantities not yet known.

    bool check_consistent = true;

    if (new_t_vir <= 0) {
	if (new_m_tot > 0 && new_r_vir > 0) {

	    // If no time unit is explicitly specified, derive it from
	    // the other units.

	    // Standard (N-body / Heggie & Mathieu) time unit, in Myr:

	    new_t_vir = 21.0 * sqrt(q_vir/new_m_tot) * pow(new_r_vir, 1.5);
	    check_consistent = false;
	}
	else
	    err_exit("addstar: Unable to set time scaling");
    }

    if (new_r_vir <= 0) {
	if (new_m_tot > 0 && new_t_vir > 0) {

	    // If no length unit is explicitly specified, derive it from
	    // the other units.

	    new_r_vir = pow(new_t_vir / (21.0 * sqrt(q_vir/new_m_tot)), 2.0/3);
	    check_consistent = false;
	}
	else
	    err_exit("addstar: Unable to set radius scaling");
    }

    if (check_consistent) {

	// Check that the various units are consistent.

	real tv = 21.0 * sqrt(q_vir*pow(new_r_vir, 3)/new_m_tot);
	if (abs(1-new_t_vir/tv) > 1.e-4) {
	    cerr << "Warning: inconsistent units: ";
	    PRC(new_t_vir);
	    PRL(tv);
	}
    }

    // Finally, set the scalings and create star parts.

    if (old_m_tot > 0 && !M_flag) {

	// Existing dynamical masses are to be interpreted as
	// solar units.  Set scaling accordingly before setting
	// up star parts.

	b->get_starbase()->set_stellar_evolution_scaling(1,1,1);
	addstar(b, T_start, type, true);
//	addstar(b, T_start, Main_Sequence, true);

	// New scaling assumes that total dynamical mass will be
	// rescaled to 1.

	b->get_starbase()->set_stellar_evolution_scaling(new_m_tot,
							 new_r_vir,
							 new_t_vir);
    }
    else {

	// New scaling assumes that total dynamical mass will be
	// rescaled to 1.

	b->get_starbase()->set_stellar_evolution_scaling(new_m_tot,
							 new_r_vir,
							 new_t_vir);
//	addstar(b, T_start, Main_Sequence, false);
	addstar(b, T_start, type, false);
    }

#if 0
    real m_sum = 0;
    for_all_daughters(node, b, bi) {
      m_sum += bi->get_mass();
    }
    b->set_mass(m_sum);

    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
    if(old_mtot!=m_sum) {
      real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
      real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
      b->get_starbase()->set_stellar_evolution_scaling(m_sum,
						       old_r_vir,
						       old_t_vir);
    }
#endif

    put_node(b);	
    delete b;
}

#endif

/* endof: addstar.c */
