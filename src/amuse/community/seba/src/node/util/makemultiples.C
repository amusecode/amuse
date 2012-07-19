
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Create binary secondary, triple tertiary, quadruple quadruple,
//// quintuple quintuple, etc, etc, etc. components for randomly
//// selected stars in the input snapshot, placing the results in a
//// hierarchical tree.  Only the masses of the hierarchical components 
//// are set here; orbital parameters are set by dyn::makeorbits.  No
//// attempt is made to maintain any existing mass or energy scaling.
//// Use scale after this function is called, but before makeorbits, if
//// it is desired to specify energies in kT units. If it is desired to
//// scale orbital parameters in stellar units call add_stars before
//// makeorbits.
////
//// Usage: makemultiples [OPTIONS]
////
//// Options:
////               -d    specify depth: [2] 2: allow binary hierarchies;
////                     3: allow triples, etc.
////               -f    specify binary fraction [0.1]
////               -h    speficy fraction of hierarchies [0.1]
////                     fraction of triples, quadruples etc, are all 
////                     relative to the fraction of the lower order system.
////                     i.e., if -f 0.5 -h 0.2, then 20% of the binaries 
////                     will be triples, 20% of the triples will be 
////                     quadruples, etc.
////               -i    use (a,b) as component indices [false]
////               -l    specify lower limit on mass ratio secondary mass [0]
////               -M    specify upper limit for primaries to be binaries [inf]
////               -m    specify lower limit for primaries to be binaries [0]
////               -p    split fraction of primaries [0.5]
////               -q    select choice of minimum mass ratio [false]:
////                     if true, secondary mass ratio is chosen uniformly
////                     on [lower_limit, upper_limit];
////                     if false, secondary mass is chosen uniformly on
////                     on [mmin, primary_mass], where mmin and mmax are
////                     specified on the command line
////               -S    split primary star [false]
////               -s    specify random seed [random from system clock]
////               -u    specify upper limit on mass ratio or 
////                     secondary mass [1 or m_primary]
////
//// Written by Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//   Simon Portegies Zwart, Amsterdam, 6 June 2003

#include "node.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

local void name_from_components(node *od, node *yd)
{
    char name1[256], name[256];
    strcpy(name1, od->format_label());
    sprintf(name, "(%s,%s)", name1, yd->format_label());
    od->get_parent()->set_name(name);
    od->get_parent()->set_index(-1);
}

// split_particle: Split the specified node into a binary with the specified
//                 parameters.  All unspecified orbit elements are chosen
//                 randomly.  Newly created leaves have names "n1" and "n2",
//                 where "n" is the name of the leaf being split.

local real split_particle(node* bi, real mass_ratio, real fmultiple, 
			  int depth, int &nmax, bool force_index, bool split,
			  real primary_fraction)
{
  if(depth<=1)
    return true;

    if (bi->get_oldest_daughter() != NULL) {
	cerr << "Can't split a binary node!\n";
	return false;
    } else if (mass_ratio <= 0 || mass_ratio > 1) {
	cerr << "Must specify mass ratio in (0,1]!\n";
	return false;
    }

    // Update the binary tree structure:

    node* d1 = new node;
    node* d2 = new node;

    bi->set_oldest_daughter(d1);

    d1->set_parent(bi);
    d2->set_parent(bi);

    d1->set_younger_sister(d2);
    d2->set_elder_sister(d1);

    // Set new masses and radii:

    real m_total = bi->get_mass();
    real m1;
    if(split) 
      m1 = m_total / (1 + mass_ratio);
    else 
      m1 = m_total;
    real m2 = m1 * mass_ratio;

    // By convention, the first component has larger mass.

    if (m1 < m2) {
	real temp = m1;
	m1 = m2;
	m2 = temp;
    }

    d1->set_mass(m1);
    d2->set_mass(m2);

    // Naming convention:

    if (force_index) {

      	d1->set_index(bi->get_index());
      	d2->set_index(bi->get_index()+nmax);	
      //      d1->set_name(bi->get_name());
      //      d2->set_name(bi->get_name());
      //      strcat(d1->get_name(), "1");
      //      strcat(d2->get_name(), "2");

    } else {

	if (bi->get_name() == NULL) {

	    // Make a name for use by the components.

	    char tmp[128];
	    if (bi->get_index() >= 0)
		sprintf(tmp, "%d", bi->get_index());
	    else
		sprintf(tmp, "X");
	    bi->set_name(tmp);
	}

	if (bi->get_name() != NULL) {

	    // Label components "a" and "b".

	    d1->set_name(bi->get_name());
	    d2->set_name(bi->get_name());
	    strcat(d1->get_name(), "a");
	    strcat(d2->get_name(), "b");
	}
    }
    //#if 0
    if (bi->get_name() == NULL)
	if (bi->get_index() >= 0) {
	    char tmp[64];
	    sprintf(tmp, "%d", bi->get_index());
	    bi->set_name(tmp);
	}

    d1->set_name(bi->get_name());
    d2->set_name(bi->get_name());
    if (force_index) {
      strcat(d1->get_name(), "1");
      strcat(d2->get_name(), "2");
    } else {
      strcat(d1->get_name(), "a");
      strcat(d2->get_name(), "b");
    }
    //#endif
    name_from_components(d1, d2);
    
    if(randinter(0, 1)<fmultiple) {
      real q = randinter(0, 1);
      if(randinter(0, 1)<primary_fraction) {
      split_particle(d1, q, fmultiple, depth-1, nmax, force_index, split,
		     primary_fraction);
      }
      else {
	split_particle(d2, q, fmultiple, depth-1, nmax, force_index, split,
		     primary_fraction);
      }
    }

    return d2->get_mass();
}

#if 0
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

local real add_secondary(node* original, real mass_ratio,
			 real fmultiple, int max_depth,
			 bool force_index, int &nmax,
			 bool split)
{

  if (!split_particle(original, mass_ratio, fmultiple, max_depth, nmax, force_index, split, primary_fraction)
)
    err_exit("split_particles: fatal error");

  return 1;
}
#endif

local void warning(real mmin, real lower_limit) {

  cerr << "WARNING:" << endl;
  cerr << "       local void makesecondary(node*, real, real, bool)" << endl;
  cerr << "       actual minimum mass " << mmin
       << " is larger than lower limit " << lower_limit << ";" << endl;
  cerr << "       continuing execution with true minimum mass" << endl;

}

local real random_mass_ratio(const real qmin, 
			     const real qmax) {

  // For now, equal probability in q between qmin and qmax.

  return randinter(qmin, qmax);
}


local void name_recursive(node *b)
{
    // Make a standard CM name, making sure that we build names
    // from the bottom up.

    node *od = b->get_oldest_daughter();
    if (od) {
	node *yd = od->get_younger_sister();
	name_recursive(od);
	name_recursive(yd);
	name_from_components(od, yd);
    }
}

local int check_indices(node *b)
{
    int nmax = 0;
    bool all_ok = true;

    for_all_leaves(node, b, bb) {

	bb->set_name(NULL);

	if (bb->get_index() >= 0)
	    nmax = Starlab::max(nmax, bb->get_index());
        else
	    all_ok = false;
    }

    if (!all_ok) {

	cerr << "    makesecondary: assigning numeric indices to nodes" << endl;

	for_all_leaves(node, b, bb) {
	    if (bb->get_index() <= 0)
		bb->set_index(++nmax);
	}
    }

    // Redetermine all center-of-mass names.

    for_all_nodes(node, b, bb)
	if (bb!= b && bb->is_parent())
	    name_recursive(bb);

    return nmax;
}

local void makemultiples(node* b, real binary_fraction,
			 real min_mprim, real max_mprim,
			 real lower_limit, real upper_limit, 
			 real fmultiple, int max_depth,
			 bool force_index, bool q_flag, bool ignore_limits,
			 bool split, real primary_fraction)
{
    PRI(4); PRL(binary_fraction);
    PRI(4); PRL(min_mprim);
    PRI(4); PRL(max_mprim);
    PRI(4); PRL(q_flag);
    PRI(4); PRL(lower_limit);
    PRI(4); PRL(upper_limit);
    PRI(4); PRL(ignore_limits);
    PRI(4); PRL(split);
    PRI(4); PRL(primary_fraction);

    renumber(b, 1, true, true, true);


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

    // Removed sum-method as it is unreliable for most fractions 
    // 180504 SPZ&JvdB
    //  real sum = 0;
    b->set_mass(0);

    real m_prim, q;
    real m_primaries=0, m_secondaries=0, m_total=0;

    if (q_flag) {

        // Choose the secondary mass ratio uniformly on
	// [lower_limit, upper_limit].

	for_all_daughters(node, b, bi) {
	  if(!bi->is_parent()) {

		m_prim = bi->get_mass();
		m_primaries += m_prim;

		if (m_prim >= min_mprim && m_prim <= max_mprim) {
		    if(randinter(0, 1)<binary_fraction) {
		      // sum += binary_fraction;
		      // if (sum >= 0.9999999999) {// avoid roundoff problems!
		      //    sum = 0;
			q = random_mass_ratio(lower_limit, upper_limit);
			
			m_secondaries += split_particle(bi, q, fmultiple, 
							max_depth, nmax, 
							force_index, split,
							primary_fraction);
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

	    for_all_daughters(node, b, bi) {
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

	for_all_daughters(node, b, bi) {
	  if(!bi->is_parent()) {

		m_prim = bi->get_mass();
		m_primaries += m_prim;

		if (m_prim >= min_mprim && m_prim <= max_mprim) {
		    if(randinter(0, 1)<binary_fraction) {
		      // sum += binary_fraction;
		      //  if (sum >= 0.9999999999) {	// avoid roundoff problems!
		      //    sum = 0;
			qmin = mmin/bi->get_mass();
			qmax = Starlab::min(mmax, bi->get_mass())/bi->get_mass();
			q = random_mass_ratio(qmin, qmax);
			m_secondaries += split_particle(bi, q, fmultiple, 
							max_depth, nmax, 
							force_index, split,
							primary_fraction);
		    }
		}
	  }
	}
    }

    m_total = 0;
    for_all_leaves(node, b, bb) {
      m_total += bb->get_mass();
    }

    // this is not smart, but it works..
    real m_sub;
    for_all_nodes(node, b, bb) {
      m_sub = 0;
      if(bb->is_parent()) {
	for_all_leaves(node, bb, bi) {
	  m_sub += bi->get_mass();
	}
	bb->set_mass(m_sub);
      }
    }


    m_secondaries = m_total - m_primaries;
    PRL(m_secondaries);PRC(m_total);PRL(m_primaries);

    if(split) {
	m_total = m_primaries;
	m_primaries -= m_secondaries;
    }
    b->set_mass(m_total);	
    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);

    if(old_mtot>0 && old_mtot<m_total) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_total,
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
	putrq(b->get_log_story(), "initial_mass", m_total);

    for_all_daughters(node, b, bi) {
      construct_node_name(bi);
    }
}

int main(int argc, char ** argv)
{

    bool split = false;
    bool ignore_limits = false;
    real binary_fraction = 0.1;
    real primary_fraction = 0.5; // split half the primaries
    real fmultiple = 0.1;
    real lower_limit = 0.0;
    bool u_flag = false;
    real upper_limit = VERY_LARGE_NUMBER;
    real max_mprim = VERY_LARGE_NUMBER;
    real min_mprim = 0;
    bool force_index = true;
    real q_flag = false;	     // flag for mass-ratio selection
    int random_seed = 0;
    char seedlog[64];

    int max_depth = 1;  // only allow binaries;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "d:qf:h:M:m:il:nu:p:Ss:I";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c) {

	    case 'f': binary_fraction = atof(poptarg);
		      break;
	    case 'h': fmultiple = atof(poptarg);
		      break;
	    case 'd': max_depth = atoi(poptarg);
		      break;
	    case 'i': force_index = !force_index;
		      break;
	    case 'I': ignore_limits = true;
		      break;
	    case 'l': lower_limit = atof(poptarg);
		      break;
	    case 'M': max_mprim = atof(poptarg);
		      break;
	    case 'm': min_mprim = atof(poptarg);
		      break;
            case 'p': primary_fraction = atoi(poptarg);
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

    if(fmultiple>0) 
      max_depth = Starlab::max(2, max_depth);

    if (!u_flag && q_flag) upper_limit = 1;

    if (binary_fraction < 0 || binary_fraction > 1)
        err_exit("makemultiples: Illegal binary fraction");
    if (q_flag && lower_limit > 1)
        err_exit("makemultuples: Illegal minimum mass for secondary");
    if (q_flag && lower_limit < 0) {
        cerr << "Negative mass ratio not allowed; use -l 0" << endl;
	lower_limit = 0;
    }
    if (q_flag && upper_limit > 1) {
        cerr << "Maximum mass ratio > 1 not allowed; use -u 1" << endl;
	upper_limit = 1;
    }
    if(primary_fraction<0 || primary_fraction>1) 
      err_exit("makemultuples: Illegal primary_fraction");

    if (lower_limit < 0)
        err_exit("Invalid lower limit");

    if (upper_limit < 0)
        err_exit("Invalid upper limit");
    else if (lower_limit > upper_limit)
        err_exit("Invalid selection of upper and lower limits to mass-ratio");

    if (max_mprim < min_mprim)
        err_exit("Invalid upper and lower primary mass selections");

    node *b;
    b = get_node();

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d", actual_seed);
    b->log_comment(seedlog);

    makemultiples(b, binary_fraction, min_mprim, max_mprim, 
		  lower_limit, upper_limit, 
		  fmultiple, max_depth, 
		  force_index, q_flag,
		  ignore_limits, split, primary_fraction);

    put_node(b);
    rmtree(b);
    return 0;
}
#endif
