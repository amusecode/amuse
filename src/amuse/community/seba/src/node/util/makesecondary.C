
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Create binary secondary components for randomly selected stars
//// in the input snapshot, placing the results in a binary tree.
//// Note that only top-level nodes are affected, and the binary
//// fraction is determined relative to the number of top-level nodes.
//// This function may be called several times to create several
//// classes of binaries in a single system.
////
//// Only secondary masses are set here; orbital parameters are set
//// by dyn::makebinary.
////
//// No attempt is made to maintain any existing mass or energy scaling.
//// Use scale after this function is called, but before makebinary, if
//// it is desired to specify energies in kT units.
////
//// Usage: makesecondary [OPTIONS] < input > output
////
//// Options:
////          -C    force output in 'col' format (dyn version only) [no]
////                (note: loses binary data)
////          -f    specify binary fraction for top-level stars having
////                masses greater than or equal to the value specified
////                with the "-m" option [0.1]
////          -i    use (a,b) as component indices [false]
////          -I    don't limit masses to primary mass range [false]
////          -l    specify lower limit on mass ratio or secondary mass [0]
////          -M    specify upper limit for primaries to be binaries [inf]
////          -m    specify lower limit for primaries to be binaries [0]
////          -q    select choice of minimum mass ratio [false]:
////                if true, secondary mass ratio is chosen uniformly
////                on [lower_limit, upper_limit];
////                if false, secondary mass is chosen uniformly on
////                [mmin, primary_mass], where mmin and mmax are
////                specified on the command line
////          -S    split primary star [false]
////          -s    specify random seed [random from system clock]
////          -u    specify upper limit on mass ratio or 
////                secondary mass [1 or m_primary]
////
//// Written by Steve McMillan and Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//   Steve McMillan, July 1996
//   Simon Portegies Zwart, December 1997

#ifndef DYN
#   include "node.h"
#   define NODE node
#else
#   include "dyn.h"
#   define NODE dyn
#endif

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

local void name_from_components(NODE *od, NODE *yd)
{
    char name1[256], name[256];
    strcpy(name1, od->format_label());
    sprintf(name, "(%s,%s)", name1, yd->format_label());
    od->get_parent()->set_name(name);
    od->get_parent()->set_index(-1);
}

// add_secondary: Add a secondary component to a primary to make
//		  a binary.  Only the mass, name and tree structure
//		  are modified at this stage.
//
//		  The old default naming was for primary and secondary
//		  components to have names "na" and "nb", respectively,
//		  where "n" was the name of the original star.  However,
//		  this leads to unique_id problems with worldlines, so
//		  numeric names are now the default, with (a,b) the
//		  alternative.  (Steve, 5/01)

local real add_secondary(NODE* original, real mass_ratio,
			 bool force_index, int &nmax,
			 bool split)
{
    NODE* primary = new NODE;
    NODE* secondary = new NODE;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    if(!split) {
	primary->set_mass(original->get_mass());
	secondary->set_mass(mass_ratio*original->get_mass());
	original->inc_mass(secondary->get_mass());
    }
    else {
	real m_total = original->get_mass();
	real m_prim = m_total / (1+mass_ratio);
	real m_sec = m_total - m_prim;
	if(m_prim<m_sec) {
	    real tmp = m_sec;
	    m_sec = m_prim;
	    m_prim = tmp;
	}

	primary->set_mass(m_prim);
	secondary->set_mass(m_sec);
    }

    // There are two possible naming conventions.  If force_index is
    // false, then we take the original name and add "a" and "b" for
    // the component names.  If force_index is true, the first component
    // inherits the original index and the second component gets an
    // index offset by some "resonable" amount.  In either case, the new
    // CM name is (cpt1,cpt1), as usual.

    if (force_index) {

	primary->set_index(original->get_index());
	secondary->set_index(original->get_index()+nmax);

    } else {

	if (original->get_name() == NULL) {

	    // Make a name for use by the components.

	    char tmp[128];
	    if (original->get_index() >= 0)
		sprintf(tmp, "%d", original->get_index());
	    else
		sprintf(tmp, "X");
	    original->set_name(tmp);
	}

	if (original->get_name() != NULL) {

	    // Label components "a" and "b".

	    primary->set_name(original->get_name());
	    secondary->set_name(original->get_name());
	    strcat(primary->get_name(), "a");
	    strcat(secondary->get_name(), "b");
	}
    }

    // Make standard CM name.

    name_from_components(primary, secondary);

    return secondary->get_mass();
}

local void warning(real mmin, real lower_limit) {

  cerr << "WARNING:" << endl;
  cerr << "       makesecondary:" << endl;
  cerr << "       actual minimum mass " << mmin
       << " is larger than lower limit " << lower_limit << ";" << endl;
  cerr << "       continuing execution with true minimum mass" << endl;

}

local real random_mass_ratio(const real qmin, 
			     const real qmax) {

  // For now, equal probability in q between qmin and qmax.

  return randinter(qmin, qmax);
}

local void name_recursive(NODE *b)
{
    // Make a standard CM name, making sure that we build names
    // from the bottom up.

    NODE *od = b->get_oldest_daughter();
    if (od) {
	NODE *yd = od->get_younger_sister();
	name_recursive(od);
	name_recursive(yd);
	name_from_components(od, yd);
    }
}

local int check_indices(NODE *b)
{
    int nmax = 0;
    bool all_ok = true;

    for_all_leaves(NODE, b, bb) {

	bb->set_name(NULL);

	if (bb->get_index() >= 0)
	    nmax = Starlab::max(nmax, bb->get_index());
        else
	    all_ok = false;
    }

    if (!all_ok) {

	cerr << "    makesecondary: assigning numeric indices to nodes" << endl;

	for_all_leaves(NODE, b, bb) {
	    if (bb->get_index() <= 0)
		bb->set_index(++nmax);
	}
    }

    // Redetermine all center-of-mass names.

    for_all_nodes(NODE, b, bb)
	if (bb!= b && bb->is_parent())
	    name_recursive(bb);

    return nmax;
}

local void makesecondary(NODE* b, real binary_fraction,
			 real min_mprim, real max_mprim,
			 real lower_limit, real upper_limit, 
			 bool force_index, bool q_flag, bool ignore_limits,
			 bool split)
{
    cerr << "makesecondary:" << endl;
    PRI(4); PRL(binary_fraction);
    PRI(4); PRL(min_mprim);
    PRI(4); PRL(max_mprim);
    PRI(4); PRL(q_flag);
    PRI(4); PRL(lower_limit);
    PRI(4); PRL(upper_limit);
    PRI(4); PRL(ignore_limits);
    PRI(4); PRL(split);

    int nmax = 0;
    if (force_index) {

	// All stars will are to have numeric labels.  First check
	// that the original labels are numeric.  If not, force all
	// labels into numeric form.  Function check_indices() returns
	// the maximum (corrected) index found.

	nmax = check_indices(b);

	// Round nmax up to something "reasonable".

	nmax = (int)(pow(10., ((int)log10((real)nmax)) + 1.0) + 0.1);
    }

    // For now, use a flat distribution in secondary mass ratio.
    // Assume that the probability of a star being the primary of
    // a binary is independent of mass.

    real m_prim, q;
    real m_primaries=0, m_secondaries=0, m_total=0;
    real sum = 0;

    int n_binaries = 0, n_target = 0;
    for_all_daughters(NODE, b, bi) {
	m_prim = bi->get_mass();
	if (m_prim >= min_mprim && m_prim <= max_mprim) n_target++;
    }
    real x_target = n_target*binary_fraction;
    n_target = (int)x_target;

    if (q_flag) {

        // Choose the secondary mass ratio uniformly on
	// [lower_limit, upper_limit].

	for_all_daughters(NODE, b, bi) {	// top level nodes only
	    if(!bi->is_parent()) {

		m_prim = bi->get_mass();
		m_primaries += m_prim;

		if (m_prim >= min_mprim && m_prim <= max_mprim) {
		    sum += binary_fraction;
		    if (sum >= 0.9999999999) {	// avoid roundoff problems!
			sum -= 1;		// thanks to Karim-Pierre MAALEJ
						// 9/04
			q = random_mass_ratio(lower_limit, upper_limit);
			
			m_secondaries += add_secondary(bi, q, force_index,
						       nmax, split);
			n_binaries++;
		    }
		}
	    }
	}

    } else {

	// Secondary mass ratio chosen uniformly on [mmin, mmax],
	// where mmin and mmax are the lower and upper limits
	// specified  on the command line.

	real mmin = VERY_LARGE_NUMBER, mmax = 0, qmin = 0, qmax;

	if (ignore_limits) {

	    mmin = lower_limit;
	    mmax = upper_limit;

	} else {

	    for_all_daughters(NODE, b, bi) {
		if (bi->is_leaf()) {
		    mmin = Starlab::min(mmin, bi->get_mass());
		    mmax = Starlab::max(mmax, bi->get_mass());
		}
	    }

	    if (mmin < lower_limit) {
		warning(mmin, lower_limit);
		sprintf(tmp_string,
			"         -l = %5.4f adopted", mmin); 
		b->log_comment(tmp_string);
	    }

	    if (mmax > upper_limit) {
		warning(mmax, upper_limit);
		sprintf(tmp_string,
			"         -u = %5.4f adopted", mmax); 
		b->log_comment(tmp_string);
	    }
	}
	    
	// Note that the following loop is essentially the same code as above.

	for_all_daughters(NODE, b, bi) {	// top level nodes only
	    if(!bi->is_parent()) {

		m_prim = bi->get_mass();
		m_primaries += m_prim;

		if (m_prim >= min_mprim && m_prim <= max_mprim) {
		    sum += binary_fraction;
		    if (sum >= 0.9999999999) {	// avoid roundoff problems!
			sum -= 1;
			qmin = mmin/bi->get_mass();
			qmax = Starlab::min(mmax,
					    bi->get_mass())/bi->get_mass();
			q = random_mass_ratio(qmin, qmax);
			m_secondaries += add_secondary(bi, q, force_index,
						       nmax, split);
			n_binaries++;
		    }
		}
	    }
	}
    }

    if(split) m_primaries -= m_secondaries;
    PRI(4); PRC(n_target); PRL(n_binaries);

    // Recompute the total system mass.

    real m_tot = total_mass(b);
    b->set_mass(m_tot);	
    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);

    if(old_mtot > 0 && old_mtot < m_tot) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_tot,
							 old_r_vir,
							 old_t_vir);
    }

    if(split) {
	sprintf(tmp_string,
		"       total mass in primaries = %8.2f Solar", m_primaries); 
	b->log_comment(tmp_string);
    }

    sprintf(tmp_string,
	    "       total mass in secondaries = %8.2f Solar", m_secondaries); 
    b->log_comment(tmp_string);

    // Update mass info, if any.

    if (find_qmatch(b->get_log_story(), "initial_mass"))
	putrq(b->get_log_story(), "initial_mass", m_tot);

    // Erase any info on virial radius or total energy, as we
    // do not recompute these quantities here.

    if (find_qmatch(b->get_log_story(), "initial_total_energy"))
        rmq(b->get_log_story(), "initial_total_energy");
    if (find_qmatch(b->get_log_story(), "initial_rvirial"))
        rmq(b->get_log_story(), "initial_rvirial");
}

int main(int argc, char ** argv)
{
#ifdef DYN
    bool C_flag = false;
#endif

    bool split = false;
    bool ignore_limits = false;
    real binary_fraction = 0.1;
    real lower_limit = 0.0;
    bool u_flag = false;
    real upper_limit = VERY_LARGE_NUMBER;
    real max_mprim = VERY_LARGE_NUMBER;
    real min_mprim = 0;
    bool force_index = true;
    real q_flag = false;	     // flag for mass-ratio selection
    int random_seed = 0;
    char seedlog[64];

    check_help();

    extern char *poptarg;
    int c;

#ifndef DYN
    const char *param_string = "qf:M:m:il:nu:Ss:I";
#else
    const char *param_string = "Cqf:M:m:il:nu:Ss:I";
#endif

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.19 $", _SRC_)) != -1)
	switch(c) {

#ifdef DYN
	    case 'C': C_flag = true;
		      break;
#endif
	    case 'f': binary_fraction = atof(poptarg);
		      break;
	    case 'i': force_index = true;
		      break;
	    case 'I': ignore_limits = true;
		      break;
	    case 'l': lower_limit = atof(poptarg);
		      break;
	    case 'M': max_mprim = atof(poptarg);
		      break;
	    case 'm': min_mprim = atof(poptarg);
		      break;
            case 'q': q_flag = true;
		      break;
	    case 'S': split = true;
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
	    case 'u': upper_limit = atof(poptarg);
	    	      u_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (!u_flag && q_flag) upper_limit = 1;

    if (binary_fraction < 0 || binary_fraction > 1)
        err_exit("makesecondary: Illegal binary fraction");
    if (q_flag && lower_limit > 1)
        err_exit("makesecondary: Illegal minimum mass for secondary");
    if (q_flag && lower_limit < 0) {
        cerr << "Negative mass ratio not allowed; use -l 0" << endl;
	lower_limit = 0;
    }
    if (q_flag && upper_limit > 1) {
        cerr << "Maximum mass ratio > 1 not allowed; use -u 1" << endl;
	upper_limit = 1;
    }

    if (lower_limit < 0)
        err_exit("Invalid lower limit");

    if (upper_limit < 0)
        err_exit("Invalid upper limit");
    else if (lower_limit > upper_limit)
        err_exit("Invalid selection of upper and lower limits to mass-ratio");

    if (max_mprim < min_mprim)
        err_exit("Invalid upper and lower primary mass selections");

    NODE *b;
#ifndef DYN
    b = get_node();
#else
    b = get_dyn();
#endif

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d", actual_seed);
    b->log_comment(seedlog);

    makesecondary(b, binary_fraction, min_mprim, max_mprim, 
		lower_limit, upper_limit, force_index, q_flag,
		ignore_limits, split);

#ifdef DYN
    b->set_col_output(false);
    if (C_flag) b->set_col_output(true);
    put_dyn(b);
#else
    put_node(b);
#endif

    rmtree(b);
    return 0;
}
#endif
