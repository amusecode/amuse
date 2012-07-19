
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// dyn_stats.C:  Helper functions used mainly by sys_stats.
//
// Global functions:
//
//	real print_binary_params()
//	real get_total_energy()
//	real get_semi_major_axis()
//	real get_period()
//	void get_total_energy_and_period()
//	void initialize_kepler_from_dyn_pair()
//	void print_binary_from_dyn_pair()
//	real print_structure_recursive()
//	real print_structure_recursive()
//	void compute_core_parameters()
//
// Note that all "verbose" arguments to these functions are bool.

#include "dyn.h"

real print_binary_params(kepler* k, real m1, real kT,
			 real dist_from_center,
			 bool verbose,			// default = true
			 bool long_binary_output,	// default = true
			 int init_indent,		// default = 0
			 int indent)			// default = BIN_INDENT

// Print out the orbital parameters of a binary, regardless of energy.
// Return the total energy of the binary system.

{
    real m = k->get_total_mass();
    real m2 = m - m1;
    real mu = m1 * m2 / m;
    real mu_kT;

    if (kT > 0) mu_kT = mu / kT;

    if (verbose) {

	for (int i = 0; i < init_indent-2; i++) cerr << " ";
	cerr << "  ";

	if (!long_binary_output) {

	    int p = cerr.precision(LOW_PRECISION);

	    cerr << "a= " << k->get_semi_major_axis()
		 << " e= " << k->get_eccentricity()
		 << " r= " << k->get_separation();

	    if (kT > 0) cerr << " E/kT= " << mu_kT*k->get_energy();

	    cerr << " D= " << dist_from_center;
	    cerr << " n= " << k->get_normal_unit_vector();

	    cerr << " n= " << k->get_normal_unit_vector();

	    // cerr << endl;	// newline should be added by calling
				// function (e.g. print_pert)

	    cerr.precision(p);

	} else {

	    cerr << "a = " << k->get_semi_major_axis()
		 << "  e = " << k->get_eccentricity()
		 << endl;

	    PRI(indent); cerr << "r = " << k->get_separation()
			      << "  D = " << dist_from_center
			      << endl;
	    
	    PRI(indent); cerr << "n = " << k->get_normal_unit_vector()
			      << endl;

	    PRI(indent); cerr << "P = " << k->get_period()
			      << "  peri = " <<  k->get_periastron();
	    if (k->get_energy() < 0)
		cerr << "  apo = " <<  k->get_semi_major_axis()
		    			* (1 + k->get_eccentricity());
	    cerr << endl;

	    PRI(indent); cerr << "E = " << mu*k->get_energy()
			      << "  E/mu = " << k->get_energy();
	    if (kT > 0) cerr << "  E/kT = " << mu_kT*k->get_energy();
	    cerr << endl;

	    PRI(indent); cerr << "masses  " << m1 << "  " << m2
			      << "  (total = " << m << ")" << endl;
	}

    } else {

	cerr << "  "  << m << "  " << m1 << "  " << m2
	     << "  "  << k->get_semi_major_axis()
	     << "  " << k->get_eccentricity()
	     << "  " << mu*k->get_energy()
	     << "  " << k->get_energy();
	if (kT > 0) cerr << "  " << -mu_kT*k->get_energy();
	cerr << "  " << dist_from_center
	     << endl;
    }

    return mu*k->get_energy();
}

real get_total_energy(dyn* bi, dyn* bj)
{
    real M = bi->get_mass() + bj->get_mass();
    vec dx = bj->get_pos() - bi->get_pos();
    vec dv = bj->get_vel() - bi->get_vel();
    real mu = (bi->get_mass() * bj->get_mass()) / M;
    return mu*(0.5*dv*dv - M / abs(dx));
}

real get_semi_major_axis(dyn* bi, dyn* bj)
{
    real M = bi->get_mass() + bj->get_mass();
    vec dx = bj->get_pos() - bi->get_pos();
    vec dv = bj->get_vel() - bi->get_vel();
    real E = 0.5*dv*dv - M / abs(dx);
    if (E == 0)
	return VERY_LARGE_NUMBER;
    else
	return -0.5*M/E;
}

real get_period(dyn* bi, dyn* bj)
{
    real M = bi->get_mass() + bj->get_mass();
    real E = -M / abs(bi->get_pos() - bj->get_pos())
			+ 0.5 * square(bi->get_vel() - bj->get_vel());
    if (E < 0) {

	real a = -0.5 * M / E;
	real P = TWO_PI * sqrt(a*a*a/M);

    	return P;

    } else

	return VERY_LARGE_NUMBER;
}

void get_total_energy_and_period(dyn* bi, dyn* bj, real& E, real& P)
{
    real M = bi->get_mass() + bj->get_mass();

    E = -M / abs(bi->get_pos() - bj->get_pos())
			+ 0.5 * square(bi->get_vel() - bj->get_vel());
    if (E < 0) {

	real a = -0.5 * M / E;
	P = TWO_PI * sqrt(a*a*a/M);

    } else

	P = VERY_LARGE_NUMBER;
}

void initialize_kepler_from_dyn_pair(kepler& k, dyn* bi, dyn* bj,
				     bool minimal)
{
    k.set_time(bi->get_system_time());
    k.set_total_mass(bi->get_mass()+bj->get_mass());
    k.set_rel_pos(bj->get_pos()-bi->get_pos());
    k.set_rel_vel(bj->get_vel()-bi->get_vel());
    k.initialize_from_pos_and_vel(minimal, false);
}

void print_binary_from_dyn_pair(dyn* bi, dyn* bj,
				real kT,		// default = 0
				vec center,		// default = (0,0,0)
				bool verbose,		// default = true
				bool long_binary_output) // default = true
{
    // Print out the relative orbital parameters of a dyn pair,
    // regardless of their relative energy.

    kepler k;
    initialize_kepler_from_dyn_pair(k, bi, bj);

    dyn *primary = bi, *secondary = bj;
//    if (bj->get_mass() > bi->get_mass()) {
//	primary = bj;
//	secondary = bi;
//    }

    if (verbose) PRI(4);
    cerr << "(";
    primary->pretty_print_node(cerr);
    cerr << ",";
    secondary->pretty_print_node(cerr);
    cerr << "):  ";

    real dist_from_center = abs(primary->get_pos() - center);

    int init_indent = 11 - strlen(bi->format_label());
    init_indent -= strlen(bj->format_label()) ;

    print_binary_params(&k, primary->get_mass(), kT,
			dist_from_center, verbose, long_binary_output,
			init_indent);

    // No newline -- should be added by calling function.
}

real print_structure_recursive(dyn* bi,
			       void (*dstar_params)(dyn*),
			       int& n_unp, real& e_unp,
			       real kT,			// default = 0
			       vec center,		// default = (0,0,0)
			       bool verbose,		// default = true
			       bool long_binary_output,	// default = true
			       int indent)		// default = 0
{

    // Recursively print out the orbital parameters of a hierarchical system,
    // regardless of energy.  Update some unperturbed counters and return the
    // total energy of the system.  On entry, bi is the top-level CM.

    real eb = 0;
    if (bi->get_oldest_daughter()) {

	dyn* od = bi->get_oldest_daughter();
	dyn* yd = od->get_younger_sister();

	kepler k;
	k.set_time(0);
	k.set_total_mass(bi->get_mass());
	k.set_rel_pos(od->get_pos()-yd->get_pos());
	k.set_rel_vel(od->get_vel()-yd->get_vel());
	k.initialize_from_pos_and_vel(true, false);  // minimal, non-verbose

	dyn* primary = od;
	if (yd->get_mass() > od->get_mass()) primary = yd;

	int init_indent = BIN_INDENT;
	if (indent > 0) {
	    for (int i = 0; i < indent; i++) cerr << " ";
	    if (!od->get_kepler())
		cerr << "    ";
	    else
		cerr << "  U ";
	    cerr << bi->format_label();
	    init_indent -= 4 + strlen(bi->format_label()) + indent;
	}

	real dist_from_center = abs(primary->get_pos() - center);

	eb += print_binary_params(&k, primary->get_mass(), kT,
				  dist_from_center, verbose,
				  long_binary_output, init_indent);
	// cerr << endl;

	if (dstar_params != NULL && od->is_leaf() && yd->is_leaf())
	    dstar_params(bi);
		    
	real e = od->print_pert(long_binary_output);
	if (e != 0) {
	    n_unp++;
	    e_unp += e;
	}

	for_all_daughters(dyn, bi, bb) {
 	    eb += print_structure_recursive(bb, dstar_params,
					    n_unp, e_unp,
					    kT, center,
					    verbose, long_binary_output,
					    indent+2);
	}
    }

    return eb;
}

real print_structure_recursive(dyn* bi,
			       real kT,			// default = 0
			       vec center,		// default = (0,0,0)
			       bool verbose,		// default = true
			       bool long_binary_output,	// default = true
			       int indent)		// default = 2
{
    // Simpler interface onto print_structure_recursive()

    int idum = 0;
    real edum = 0;

    return print_structure_recursive(bi, NULL, idum, edum, kT, center,
				     verbose, long_binary_output, indent);
}

#define TTOL 1.e-12				// arbitrary tolerance

void compute_core_parameters(dyn* b, int k,
			     bool allow_n_sq_ops,
			     vec& center,
			     real& rcore, int& ncore, real& mcore)
{
    vec vel;

    // Write densities to particle dyn stories, if necessary.

    if (!twiddles(getrq(b->get_dyn_story(), "density_time"),
		  b->get_system_time(), TTOL)) {

	if (b->n_leaves() < 12) {	// don't even try...
	    ncore = 0;
	    rcore = mcore = 0;
	    return;
	}

	if (!allow_n_sq_ops) {
	    cerr << "compute_core_parameters:  no densities available "
		 << "and allow_n_sq_ops set false" << endl;
	    ncore = 0;
	    rcore = mcore = 0;
	    return;
	}

	cerr << "\n  compute_core_parameters:  computing densities\n";
	compute_density(b, k);		// (N^2 ops on the front end...)
    }

    // Compute max (not mean) density center, if necessary.

    if (!twiddles(getrq(b->get_dyn_story(), "density_center_time"),
		  b->get_system_time(), TTOL)
	|| !strcmp(getsq(b->get_dyn_story(), "density_center_type"), "mean")) {

	// Find point with maximum density (changed by SM&SPZ Oct 11, 2001)

	compute_max_cod(b, center, vel);
	// compute_mean_cod(b, center, vel);
    }

    center -= b->get_pos();		// quantities include
    vel -= b->get_pos();		// the root node

    real total_weight = 0;
    real sum = 0;
    for_all_daughters(dyn, b, bi) {
	real density = getrq(bi->get_dyn_story(), "density");

	if (density > 0) {
	    real dens2 = density * density;		// weight factor

	    total_weight += dens2;
	    sum += dens2 * square(bi->get_pos() - center);
	}
    }

    real rcore2 = 0;
    if (total_weight > 0 && sum > 0)
	rcore2 = sum/total_weight;

    rcore = sqrt(rcore2);
    ncore = 0;
    mcore = 0;

    if (rcore2 > 0) {
	for_all_daughters(dyn, b, bj)
	    if (square(bj->get_pos() - center) <= rcore2) {
		ncore++;
		mcore += bj->get_mass();
	    }
    }
}

