
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out various diagnostic statistics on the input system.  These
//// include:
////
////         system time, number, mass, mass distribution;
////         relaxation time;
////         system energy;
////         core parameters;
////         lagrangian radii for quartiles [default], for ten-percentiles,
////         and for "special" choice of Lagrangian masses, currently 0.005,
////         0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9;
////         mass distribution by lagrangian zone;
////         anisotropy by lagrangian zone;
////         binary parameters.
////
//// In addition, Lagrangian radii are written to the root dyn story.
//// If N^2_ops	is selected, then core parameters and particle densities
//// are also written to the root and particle dyn stories.
////
//// If N^2_ops is selected, Lagrangian radii are computed relative to
//// the density center, as are binary radial coordinates.  Otherwise,
//// the modified center of mass (center of mass with outliers excluded)
//// is used.
////        	
//// If GRAPE is available, the hdyn version of this function will
//// compute the energies even if the "-n" flag is set false.
////
//// Usage: sys_stats [OPTIONS] < input > output
////
//// Options:
////              -b    specify level of binary statistics [2]
////                    0: none
////                    1 (or no arg): short binary output
////                    2: full binary output
////              -e    recalculate the total energy (even if -n is no) [yes]
////              -l    specify percentile choice [2]
////                    0: quartiles (1-3)
////                    1: 10-percentiles (10-90)
////                    2 (or no arg): nonlinear Lagrangian masses
////                    (0.5, 1, 2, 5, 10, 25, 50, 75, 90%)
////              -n    perform/don't perform actions requiring O(N^2)
////                    operations (e.g. computation of energy and core
////                    radius; see -e) [no]
////              -o    pipe system to cout [no]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.


// ***  MUST make sure that the definitions of virial radius and virial  ***
// ***  equilibrium are consistent between scale and sys_stats.          ***

//   Note:  Because of the reference to libstar.a in sys_stats, the
//          only global function in this file is sys_stats itself.
//          Other "stats-related" functions are found in dyn_stats.C.
//	    Thus, any "dyn-only" tools referencing just the helper
//	    functions will not have to load the "star" libraries.

#include "dyn.h"
#include "star/sstar_to_dyn.h"

#ifndef TOOLBOX

local bool check_for_doubling(dyn* b)
{
    if (getrq(b->get_log_story(), "fraction_of_masses_doubled") > 0)
	return true;
    else
	return false;
}

local void print_relaxation_time(dyn* b,
				 real& potential_energy,
				 real& kinetic_energy,
				 real& r_virial,
				 void (*compute_energies)(dyn*, real,
							  real&, real&, real&,
							  bool))
{
    // Print the relaxation time.  Energies are returned as a side effect!

    real t_relax;
    real total_mass = 0;
    int N = b->n_leaves();

    for_all_leaves(dyn, b, bj)
	total_mass += bj->get_mass();

    real e_total;
    compute_energies(b, 0.0, potential_energy, kinetic_energy, e_total, true);

    // Potential energy here is internal potential energy.

    r_virial = -0.5 * total_mass * total_mass / potential_energy;

    // Note suspect choice of Lambda...

    t_relax = 9.62e-2 * sqrt(r_virial * r_virial * r_virial / total_mass)
              * N / log10(0.4 * N);

    PRI(4); PRC(t_relax); PRL(r_virial);
}

local void print_numbers_and_masses(dyn* b, bool& mass_spectrum)
{
    // Numbers of stars and nodes:

    real total_cpu_time = getrq(b->get_log_story(), "total_cpu_time");
    if (total_cpu_time > -VERY_LARGE_NUMBER) {
	PRI(4); PRL(total_cpu_time);
    } else {
	PRI(4); cerr << "total_cpu_time = (undefined)" << endl;
    }

    int N = b->n_leaves();
    int N_top_level = b->n_daughters();
    cerr << "    N = " << N << "  N_top_level = " << N_top_level << endl;

    // Masses and averages:

    real total_mass_nodes = 0;
    for_all_daughters(dyn, b, bi)
	total_mass_nodes += bi->get_mass();

    real total_mass_leaves = 0;
    real m_min = VERY_LARGE_NUMBER, m_max = -m_min;
    for_all_leaves(dyn, b, bj) {
	total_mass_leaves += bj->get_mass();
	m_min = Starlab::min(m_min, bj->get_mass());
	m_max = Starlab::max(m_max, bj->get_mass());
    }
    real m_av = total_mass_leaves / Starlab::max(1, N);

    cerr << "    total_mass = " << b->get_mass();
    cerr << "  nodes: " << total_mass_nodes;
    cerr << "  leaves: " << total_mass_leaves << endl;

    cerr << "    m_min = " << m_min << "  m_max = " << m_max
	 << "  m_av = " << m_av << endl;

    mass_spectrum = (m_min < m_max);
}


local void print_parameters_for_massive_black_holes(dyn *b,
						    real kT,
						    vec center,
						    bool verbose)
{
    int n_bh = 0;
    for_all_leaves(dyn, b, bi) {
	if (find_qmatch(bi->get_log_story(), "black_hole")) {

	    if (n_bh++ == 0 && verbose)
		cerr << endl
		     << "  Orbital parameters for massive black holes:"
		     << endl;

	    // Original printed out binary parameters for bi relative to b,
	    // which is usually root.  Now print basic info for top-level
	    // hole, and binary info if really a binary...  (Steve, 3/06)

	    if (bi->is_top_level_node()) {
		cerr << "    (single) " << bi->format_label()
		     << ":  m = " << bi->get_mass();
		cerr << ",  |r-CM| = " << abs(bi->get_pos() - b->get_pos());
		cerr << "  |v-CM| = " << abs(bi->get_vel() - b->get_vel())
		     << endl;
	    } else {
		bool long_binary_output = true;	
		print_binary_from_dyn_pair(bi,bi->get_binary_sister(),
					   kT, center, verbose,
					   long_binary_output);
	    }
	}
    }

//    if (!n_bh) {
//	cerr << "           "
//	     << "   ---   "
//	     << "No massive black holes"
//	     << "   ---   "
//	     << endl;
//    }
}


local int which_zone(dyn* bi, vec& center_pos, int n_lagr, real* r_lagr)
{
    vec dr = bi->get_pos() - center_pos;
    real dr2 = dr*dr;

    for (int j = 0; j < n_lagr; j++)
      if (r_lagr[j] > dr2) return j;

    return n_lagr;
}

local void print_numbers_and_masses_by_radial_zone(dyn* b, int which)
{
    if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

        // Make a list of previously computed lagrangian radii.

        int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	real *r_lagr = new real[n_lagr];

	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	// Numbers of top-level nodes (leaves in multiple case):

	int *N = new int[n_lagr+1];

	// Note that r_lagr omits the 0% and 100% radii.

	// Masses and averages:

	real *total_mass_nodes = new real[n_lagr+1];
	real *m_min = new real[n_lagr+1];
	real *m_max = new real[n_lagr+1];

	// Initialization:

	for (int i = 0; i <= n_lagr; i++) {
	    N[i] = 0;
	    total_mass_nodes[i] = 0;
	    m_min[i] = VERY_LARGE_NUMBER;
	    m_max[i] = -VERY_LARGE_NUMBER;
	    if (i < n_lagr) r_lagr[i] *= r_lagr[i];
	}

	// Use the same center as lagrad, or the geometric center
	// in case of error (shouldn't happen, as lagr_pos is written
	// whenever r_lagr is).

	vec center_pos = 0;

	if (find_qmatch(b->get_dyn_story(), "lagr_pos"))
	    center_pos = getvq(b->get_dyn_story(), "lagr_pos");
	else
	    warning(
	    "print_numbers_and_masses_by_radial_zone: lagr_pos not found");

	int N_total = 0;

	// PRL(which);

	for_all_daughters(dyn, b, bi)

	    // Count single stars and binaries separately.
	    // Count single and doubled-mass stars separately also.

	    if ( (which == 0 && bi->get_oldest_daughter() == NULL)
		|| (which == 1 && bi->get_oldest_daughter() != NULL)
                || (which == 2 && getiq(bi->get_log_story(),
					"mass_doubled") == 0 )
		|| (which == 3 && getiq(bi->get_log_story(),
					"mass_doubled") == 1) ) {
		 //		 || (which == 2)
		 //		 || (which == 3) ) {

	        // Find which zone we are in.

	        int i = which_zone(bi, center_pos, n_lagr, r_lagr);

		// Update statistics.

		if (which != 1) {
		    N[i]++;
		    N_total++;
		} else {
		    N[i] += bi->n_leaves();
		    N_total += bi->n_leaves();
		}

		total_mass_nodes[i] += bi->get_mass();
		m_min[i] = Starlab::min(m_min[i], bi->get_mass());
		m_max[i] = Starlab::max(m_max[i], bi->get_mass());
	    }

	if (N_total > 0) {

	    int i;
	    for (i = 0; i <= n_lagr; i++) {
	        cerr << "    zone " << i << "  N = " << N[i];
	        if (N[i] > 0)
		    cerr << "  m_av = " << total_mass_nodes[i]
						/ Starlab::max(N[i], 1)
		         << "  m_min = " << m_min[i]
		         << "  m_max = " << m_max[i];
		cerr << endl;
	    }

	    cerr << "\n  Cumulative:\n";

	    int  Nc = 0;
	    real m_totc = 0;
	    real m_minc = VERY_LARGE_NUMBER;
	    real m_maxc = -VERY_LARGE_NUMBER;

	    for (i = 0; i <= n_lagr; i++) {
	        Nc += N[i];
		m_totc += total_mass_nodes[i];
		m_minc = Starlab::min(m_minc, m_min[i]);
		m_maxc = Starlab::max(m_maxc, m_max[i]);
		cerr << "    zone " << i << "  N = " << Nc;
		if (Nc > 0)
		    cerr << "  m_av = "  << m_totc / Starlab::max(Nc, 1)
		         << "  m_min = " << m_minc
		         << "  m_max = " << m_maxc;
		cerr << endl;
	    }

	} else

	    cerr << "    (none)\n";

	delete [] r_lagr;
	delete [] N;
	delete [] total_mass_nodes;
	delete [] m_min;
	delete [] m_max;
    }
}

local void print_anisotropy_by_radial_zone(dyn* b, int which)
{
    if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

        // Make a list of previously computed lagrangian radii.

        int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	real *r_lagr = new real[n_lagr];
	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	real *total_weight = new real[n_lagr+1];
	real *vr2_sum = new real[n_lagr+1];
	real *vt2_sum = new real[n_lagr+1];

	// Note that r_lagr omits the 0% and 100% radii.

	// Initialization:

	for (int i = 0; i <= n_lagr; i++) {
	    total_weight[i] = 0;
	    vr2_sum[i] = 0;
	    vt2_sum[i] = 0;
	    if (i < n_lagr) r_lagr[i] *= r_lagr[i];
	}

	// Use the same center as lagrad, or the geometric center
	// in case of error (shouldn't happen, as lagr_pos is written
	// whenever r_lagr is).

	vec center_pos = 0, center_vel = 0;
	if (find_qmatch(b->get_dyn_story(), "lagr_pos"))
	    center_pos = getvq(b->get_dyn_story(), "lagr_pos");
	else
	    warning("print_anisotropy_by_radial_zone: lagr_pos not found");

	if (find_qmatch(b->get_dyn_story(), "lagr_vel"))
	    center_vel = getvq(b->get_dyn_story(), "lagr_vel");

	int N_total = 0;

	// PRL(which);

	for_all_daughters(dyn, b, bi)

	    // Count single stars and binaries separately.

	    if ( (which == 0 && bi->get_oldest_daughter() == NULL)
		|| (which == 1 && bi->get_oldest_daughter() != NULL)
		|| (which == 2 && getiq(bi->get_log_story(),
					"mass_doubled") == 0)
		|| (which == 3 && getiq(bi->get_log_story(),
					"mass_doubled") == 1) ) {

	        // Find which zone we are in.

	        int i = which_zone(bi, center_pos, n_lagr, r_lagr);

		// Update statistics.

		N_total++;

		real weight = bi->get_mass();
		weight = 1;			// Heggie/Binney & Tremaine

		vec dr = bi->get_pos() - center_pos;
		vec dv = bi->get_vel() - center_vel;

		real r2 = square(dr);
		real v2 = square(dv);
		real vr = dr*dv;

		real vr2 = 0;
		if (r2 > 0) vr2 = vr * vr / r2;

		total_weight[i] += weight;
		vr2_sum[i] += weight * vr2;
		vt2_sum[i] += weight * (v2 - vr2);		
	    }

	if (N_total > 0) {

	    int k;
	    for (k = 0; k <= n_lagr; k += 5) {
	        if (k > 0) cerr << endl;
		cerr << "           ";
		for (int i = k; i < k+5 && i <= n_lagr; i++) {
		    if (vr2_sum[i] > 0)
		        cerr << " " << 1 - 0.5 * vt2_sum[i] / vr2_sum[i];
		    else
		        cerr << "   ---   ";
		}
	    }
	    cerr << endl;

	    cerr << "\n  Cumulative anisotropy:\n";

	    real vr2c = 0;
	    real vt2c = 0;

	    for (k = 0; k <= n_lagr; k += 5) {
	        if (k > 0) cerr << endl;
		cerr << "           ";
		for (int i = k; i < k+5 && i <= n_lagr; i++) {
		    vr2c += vr2_sum[i];
		    vt2c += vt2_sum[i];
		    if (vr2c > 0)
		        cerr << " " << 1 - 0.5 * vt2c / vr2c;
		    else
		        cerr << "   ---   ";
		}
	    }
	    cerr << endl;

	} else

	    cerr << "    (none)\n";

	delete [] r_lagr;
        delete [] total_weight;
	delete [] vr2_sum;
	delete [] vt2_sum;
    }
}

// local real system_energy(dyn* b)
// {
//     real kin = 0;
//     for_all_leaves(dyn, b, bj)
// 	kin += bj->get_mass()
// 	        * square(something_relative_to_root(bj, &dyn::get_vel));
//     kin *= 0.5;
//
//     real pot = 0.0;
//     for_all_leaves(dyn, b, bi) {
// 	real dpot = 0.0;
// 	for_all_leaves(dyn, b, bj) {
// 	    if (bj == bi) break;
// 	    vec dx = something_relative_to_root(bi, &dyn::get_pos)
// 			  - something_relative_to_root(bj, &dyn::get_pos);
// 	    dpot += bj->get_mass() / abs(dx);
// 	}
// 	pot -= bi->get_mass() * dpot;
//     }
//
//     return kin + pot;
// }

local real top_level_kinetic_energy(dyn* b)
{
    // Compute energy relative to the system center of mass.

    vec cmv = vec(0);
    real mass = 0;

    for_all_daughters(dyn, b, bb) {
	mass += bb->get_mass();
	cmv += bb->get_mass()*bb->get_vel();
    }
    cmv /= mass;

    real kin = 0;
    for_all_daughters(dyn, b, bb)
        kin += bb->get_mass()*square(bb->get_vel()-cmv);

    return 0.5*kin;
}

local void print_energies(dyn* b,
			  real& potential_energy,	// top-level
			  real& kinetic_energy,		// top-level
			  real& kT,
			  void (*compute_energies)(dyn*, real,
						   real&, real&, real&,
						   bool))
{
    // Energies (top-level nodes):

    real e_total;

    if (kinetic_energy == 0)
	compute_energies(b, 0.0,
			 potential_energy, kinetic_energy,
			 e_total,
			 true);				// "true" ==> top-level

    real total_int_energy = potential_energy + kinetic_energy;	// NB internal
								// pot and kin
    real external_pot = get_external_pot(b);
    e_total = total_int_energy + external_pot;

    int ppp = cerr.precision(STD_PRECISION);

    vec com_pos, com_vel;
    compute_com(b, com_pos, com_vel);

    // Center of mass should be the root node.  Any difference
    // represents a numerical error.

    real kin_com = 0.5*b->get_mass()*square(com_vel);

    // CM quantities include the root node.  Want relative quantities here.

    com_pos -= b->get_pos();
    com_vel -= b->get_vel();

    real pot_int = potential_energy, kin_int = kinetic_energy;
    kT = kin_int / (1.5*b->n_daughters());
    real virial_ratio = -kin_int / (pot_int + get_external_virial(b));

    real pot_tot = potential_energy, kin_tot = kinetic_energy;
    e_total = total_int_energy;					// = int PE + KE

    // Recompute energies, resolving any top-level nodes.

    for_all_daughters(dyn, b, bb) {
	if (bb->is_parent()) {

	    // Tricky!!  Do the calculation only if a CM is detected.

	    compute_energies(b, 0.0, pot_tot, kin_tot, e_total, false);
						    // "false" ==> all leaves
	    break;
	}
    }

    pot_tot += external_pot;
    kin_tot += kin_com;
    e_total += external_pot + kin_com;

    int p = cerr.precision(INT_PRECISION);
    cerr << endl
	 << "    Energies: " << pot_tot << " " << kin_tot << " " << e_total
	 << endl;
    cerr.precision(p);

    if (external_pot) {
	cerr << "        external potential = " << external_pot << " (";
	print_external(b->get_external_field());
	cerr << ")" << endl;
    }

    cerr << "    CM kinetic energy = " << kin_com << endl;

    cerr << "    top-level internal potential = " << pot_int
	 << ",  kinetic = " << kin_int
	 << endl
	 << "    top-level total energy = " << pot_int + kin_int
	 << ",  virial_ratio = " << virial_ratio
	 << endl;

    // PRL(get_external_virial(b));

    cerr.precision(ppp);
}

void search_for_binaries(dyn* b,
			 real energy_cutoff,		// default = 0
			 real kT,			// default = 0
			 vec center,			// default = (0,0,0)
			 bool verbose,			// default = true
			 bool long_binary_output)	// default = true
{
    // Search for bound binary pairs among top-level nodes.

    bool found = false;

    for_all_daughters(dyn, b, bi)
	for_all_daughters(dyn, b, bj) {
	    if (bj == bi) break;

	    real E = get_total_energy(bi, bj);

	    // Convention: energy_cutoff is in terms of kT if set,
	    //		   in terms of E/mu if kT = 0.

	    if ((kT > 0 && E < -energy_cutoff*kT)
		|| (kT <= 0 && E * (bi->get_mass() + bj->get_mass())
		    	< -energy_cutoff * bi->get_mass() * bj->get_mass())) {

		print_binary_from_dyn_pair(bi, bj,
					   kT, center, verbose,
					   long_binary_output);
		cerr << endl;

		found = true;
	    }
	}

    if (!found) cerr << "    (none)" << endl;
}


// Histogram routines -- accumulate and print out key binary statistics.

// Convenient to keep these global:

#define P_MAX	0.1
#define E_MIN	0.5
#define NP	12
#define NR	10
#define NE	12

static int p_histogram[NP][NR];	// first index:  log period
				// second index: radius

static int e_histogram[NP][NR];	// first index:  log E/kT
				// second index: radius

bool have_kT = false;

local void initialize_histograms()
{
    for (int ip = 0; ip < NP; ip++)
	for (int jr = 0; jr < NR; jr++)
	    p_histogram[ip][jr] = e_histogram[ip][jr] = 0;
}

local void bin_recursive(dyn* bi,
			 vec center, real rcore, real rhalf,
			 real kT)
{
    // Note: inner and outer multiple components are counted
    //	      as separate binaries.

    have_kT = false;
    if (kT > 0) have_kT = true;

    for_all_nodes(dyn, bi, bb) {
	dyn* od = bb->get_oldest_daughter();

	if (od) {

	    // Get period of binary with CM bb.

	    bool del_kep = false;

	    if (!od->get_kepler()) {
		new_child_kepler(bb);
		del_kep = true;
	    }

	    real P = od->get_kepler()->get_period();
	    real R = abs(something_relative_to_root(bb, &dyn::get_pos)
			   - center);
	    real mu = od->get_mass() * od->get_younger_sister()->get_mass()
			/ bb->get_mass();
	    real E = -mu*od->get_kepler()->get_energy();
	    if (kT > 0) E /= kT;

	    if (del_kep) {
		delete od->get_kepler();
		od->set_kepler(NULL);
	    }

	    // Bin by P and R and by E and R.

	    if (P > 0) {

		// P bins are half dex, decreasing.
		// P > P_MAX in bin 0, P_MAX >= P > 0.316 P_MAX in bin 1, etc.

		int ip = (int) (1 - 2*log10(P/P_MAX));
		if (ip < 0) ip = 0;
		if (ip >= NP) ip = NP - 1;

		int ir = 0;
		if (R > rcore) {

		    // R bin 0 is the core.
		    // Linear binning in rcore < r < rhalf, logarithmic
		    // outside rhalf (x2); bin NR/2 has rhalf <= R < 2*rhalf.

		    if (R < rhalf) {
			real x = R/rhalf;
			ir = (int) (NR*x/2);
			if (ir < 1) ir = 1;
		    } else {
			real rlim = 2*rhalf;
			ir = NR/2;
			while (R > rlim) {
			    if (ir >= NR-1) break;
			    rlim *= 2;
			    ir++;
			}
		    }
		}

		int ie = 0;
		real elim = E_MIN;

		// E bin 0 is 0 < E < 0.5 (kT), increase by factors of 2.

		while (elim < E) {
		    ie++;
		    elim *= 2;
		    if (ie >= NE-1) break;
		}

		// Simply count binaries, for now.

		p_histogram[ip][ir]++;
		e_histogram[NE-1-ie][ir]++;	// histogram will be
						// printed backwards
	    }
	}
    }	
}

local void print_histogram(int histogram[][NR], int n, const char *label)
{
    cerr << endl
	 << "  Binary distribution by " << label << " and radius:"
	 << endl;

    PRI(13+5*n);
    cerr << "total" << endl;

    for (int jr = NR-1; jr >= 0; jr--) {

	if (jr == NR/2)
	    cerr << "   Rhalf->";
	else if (jr == 0)
	    cerr << "   core-> ";
	else
	    PRI(10);

	int sum = 0;
	for (int i = n-1; i >= 0; i--) {
	    fprintf(stderr, "%5d", histogram[i][jr]);
	    sum += histogram[i][jr];
	}
	fprintf(stderr, "  %5d\n", sum);		// total by radius
    }

    cerr << endl << "   total  ";

    // Totals by column.

    int total = 0;
    for (int i = n-1; i >= 0; i--) {
	int sum = 0;
	for (int jr = NR-1; jr >= 0; jr--) sum += histogram[i][jr];
	fprintf(stderr, "%5d", sum);
	total += sum;
    }

    fprintf(stderr, "  %5d\n", total);			// grand total
}

// Print out whatever histograms we have just accumulated.

local void print_binary_histograms()
{
    // Print period histogram...

    print_histogram(p_histogram, NP, "period");

    // ...and add caption.

    cerr << endl;
    PRI(14); cerr << "<";
    PRI((NP/2-1)*5-6); cerr << "period (0.5 dex)";
    PRI((NP/2)*5-13); cerr << "> " << P_MAX << endl;

    // Print energy histogram...

    print_histogram(e_histogram, NE, "energy");

    // ...and add caption.

    cerr << endl;
    PRI(11); cerr << "< "; fprintf(stderr, "%3.1f", E_MIN);
    PRI((NE/2-1)*5-4);
    if (have_kT)
	cerr << "E/kT (x 2)";
    else
	cerr << "E/mu (x 2)";
    PRI((NE/2)*5-8); cerr << ">" << endl;
}


local void print_binaries(dyn* b, real kT,
			  vec center, real rcore, real rhalf,
			  int verbose,
			  bool long_binary_output = true,
			  void (*dstar_params)(dyn*) = NULL)
{
    
    // Print out the properties of "known" binaries (i.e. binary subtrees),
    // and get some statistics on binaries in the core.  Note that *all*
    // binary subtrees are printed, regardless of energy.
    // Node b is the root node.

    // Printout occurs in four passes:	(0) perturbed multiples
    //					(1) unperturbed multiples
    //					(2) perturbed binaries
    //					(3) unperturbed binaries

    real eb = 0;
    real nb = 0;
    int n_unp = 0;
    real e_unp = 0;

    real mbcore = 0;
    int nbcore = 0;

    initialize_histograms();

    bool found = false;

    for (int pass = 0; pass < 4; pass++)
    for_all_daughters(dyn, b, bi) {

	dyn* od = bi->get_oldest_daughter();

	if (od &&
	    (pass == 0 && bi->n_leaves()  > 2 && od->get_kepler() == NULL) ||
	    (pass == 1 && bi->n_leaves()  > 2 && od->get_kepler() != NULL) ||
	    (pass == 2 && bi->n_leaves() == 2 && od->get_kepler() == NULL) ||
	    (pass == 3 && bi->n_leaves() == 2 && od->get_kepler() != NULL)) {

	    if (verbose > 0) {
		if (od->get_kepler() == NULL)
		    cerr << "    ";
		else
		    cerr << "  U ";
		cerr << bi->format_label();
	    }

	    int init_indent = BIN_INDENT - 4 - strlen(bi->format_label());

	    if (bi->n_leaves() > 2) {

		if (verbose > 0) {
		    for (int i = 0; i < init_indent; i++) cerr << " ";
		    cerr << "is a multiple system" << endl;
		}

		// Recursively print out the structure (indent = 0, note):

		eb += print_structure_recursive(bi,
						dstar_params,
						n_unp, e_unp,
						kT, center, verbose>0,
						long_binary_output);

		bin_recursive(bi, center, rcore, rhalf, kT);

	    } else {

		// if (verbose > 0) {
		//    if (od->get_kepler()) {
		//	PRI(init_indent);
		//	cerr << "is unperturbed (has a kepler pointer)"
		//	     << endl;
		//
		//	init_indent = BIN_INDENT;
		//    }
		// }

		bool reset = false;
		if (od->get_kepler() == NULL) {
		    new_child_kepler(bi);	// attach temporary kepler
		    reset = true;		// to oldest daughter
		}

		real dist_from_center = abs(bi->get_pos() - center);

		eb += print_binary_params(od->get_kepler(), od->get_mass(),
					  kT, dist_from_center, verbose>0,
					  long_binary_output, init_indent);
		nb++;

		// Indirect reference to dstar output:

		if (dstar_params != NULL)
		    dstar_params(bi);

		if (rcore > 0) {
		    if (dist_from_center <= rcore) {
			nbcore++;
			mbcore += bi->get_mass();
		    }
		}

		bin_recursive(bi, center, rcore, rhalf, kT);

		if (reset) {
		    delete od->get_kepler();
		    od->set_kepler(NULL);
		}

		// Virtual reference to additional binary output.
		// Must do this *after* temporary keplers have been removed.

		real e = od->print_pert(long_binary_output);
		if (e != 0) {
		    n_unp++;
		    e_unp += e;
		}
	    }

	    found = true;
	}
    }

    if (!found)

	cerr << "    (none)\n";

    else {

	if (verbose > 0) cerr << endl;
	cerr << "  Total binary energy ("
	     << nb << " binaries) = " << eb << endl;

	if (n_unp != 0) {
	    cerr << "  Total unperturbed binary energy ("
		 << n_unp << " binaries) = " << e_unp << endl;
	}

	if (rcore > 0)
	    cerr << "  nbcore = " << nbcore
		 << "  mbcore = " << mbcore << endl;

	if (verbose > 1)
	  print_binary_histograms();	// histograms are presently only
					// computed if binary output is
					// enabled
    }
}

local void print_core_parameters(dyn* b, bool allow_n_sq_ops,
				 vec& density_center, real& rcore)
{
    real mcore;
    int ncore;

    compute_core_parameters(b, 12, allow_n_sq_ops,
			    density_center, rcore, ncore, mcore);

    // Note: density_center is relative to the root node.

    cerr << "    density_center = " << density_center+b->get_pos() << endl;
    cerr << "    rcore = " << rcore << "  ncore = " << ncore
	 << "  mcore = " << mcore << endl;

    // Place the data in the root dyn story.  Use HIGH_PRECISION for
    // the time stamp.

    putrq(b->get_dyn_story(), "core_radius_time", b->get_system_time(),
	  HIGH_PRECISION);
    putrq(b->get_dyn_story(), "core_radius", rcore);
    putiq(b->get_dyn_story(), "n_core", ncore);
    putrq(b->get_dyn_story(), "m_core", mcore);
}

// Linear interpolation routine.

local real linear_interpolation(const real x,
				const real x1, const real x2,
				const real y1, const real y2) {

        real a = (y2-y1)/(x2-x1);
	real b = y1 - a*x1;

	real y = a*x + b;
	return y;
}

//-----------------------------------------------------------------------------
//
// Local functions for dominated ("disk") motion:

local dyn* dominant_mass(dyn* b, real f)
{
    // Return a pointer to the most massive top-level node contributing
    // more than fraction f of the total mass.

    real total_mass = 0, max_mass = 0;
    dyn* b_dom;

    for_all_daughters(dyn, b, bi) {
	total_mass += bi->get_mass();
	if (bi->get_mass() > max_mass) {
	    max_mass = bi->get_mass();
	    b_dom = bi;
	}
    }

    if (max_mass < f*total_mass)
	return NULL;
    else
	return b_dom;
}

typedef  struct
{
    real  radius;
    real  mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return +1;
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return -1;
    else
        return 0;
}

local void print_dominated_lagrangian_radii(dyn* b, dyn* b_dom)
{
    // Compute and print 10, 20, 30, ..., 90% radii of current system,
    // taking particle b_dom as the center.

    if (!b_dom) return;

    int n = 0;
    { for_all_daughters(dyn, b, bi)
	if (bi != b_dom) n++;
    }

    // Set up an array of (radius, mass) pairs.  Also find the total
    // mass of all nodes under consideration.

    rm_pair_ptr rm_table = new rm_pair[n];

    if (rm_table == NULL) {
	cerr << "print_dominated_lagrangian_radii: "
	     << "not enough memory left for rm_table\n";
	return;
    }

    real total_mass = 0;
    int i = 0;

    for_all_daughters(dyn, b, bi)
	if (bi != b_dom) {
	    vec dr = bi->get_pos() - b_dom->get_pos();

	    // "Radius" is 2-D (x-y) only.

	    dr[2] = 0;

	    total_mass += bi->get_mass();
	    rm_table[i].radius = abs(dr);
	    rm_table[i].mass = bi->get_mass();
	    i++;
	}

    // Sort the array by radius.

    qsort((void *)rm_table, (size_t)i, sizeof(rm_pair), compare_radii);

    // Print out the percentile radii:

    cerr << endl << "  Lagrangian radii: ";

    real cumulative_mass = 0.0;
    i = 0;

    for (int k = 1; k <= 9; k++) {

        while (cumulative_mass < 0.1*k*total_mass)
	    cumulative_mass += rm_table[i++].mass;

	cerr << " " << rm_table[i-1].radius;
	if (k == 5) cerr << endl << "                    ";
    }
    cerr << endl;
}

local void print_dominated_velocity_dispersions(dyn* b, dyn* b_dom)
{

    // Compute and print system velocity dispersions in the radial,
    // transverse, and vertical (z) directions.

    // Assume that dominated motion is primarily in the x-y plane,
    // and compute the transverse dispersion relative to keplerian
    // motion.

    if (!b_dom) return;

    real total_mass = 0;
    real v2r = 0, v2t = 0, v2z = 0;

    for_all_daughters(dyn, b, bi)
	if (bi != b_dom) {
	    vec dr = bi->get_pos() - b_dom->get_pos();
	    vec dv = bi->get_vel() - b_dom->get_vel();

	    total_mass += bi->get_mass();
	    v2z += bi->get_mass() * dv[2]*dv[2];

	    // Restrict motion to x-y plane and subtract off
	    // keplerian motion in the positive sense.

	    dr[2] = 0;
	    dv[2] = 0;
	    real r = abs(dr);

	    if (r > 0) {
		real vr = dr*dv/r;

		// Subtract off keplerian motion.

		real vkep = sqrt(b_dom->get_mass()/r);
		dv[0] -= -dr[1]*vkep/r;
		dv[1] -=  dr[0]*vkep/r;

		v2r += bi->get_mass() * vr*vr;
		v2t += bi->get_mass() * (dv*dv - vr*vr);
	    }
	}

    if (total_mass > 0)
	cerr << endl << "  Velocity dispersions:  "
	     << sqrt(v2r/total_mass) << "  "
	     << sqrt(v2t/total_mass) << "  "
	     << sqrt(v2z/total_mass)
	     << endl;
}

//-----------------------------------------------------------------------------

bool parse_sys_stats_main(int argc, char *argv[],
			  int &which_lagr,
			  bool &binaries,
			  bool &long_binary_output,
			  bool &B_flag,
			  bool &calc_e,
			  bool &n_sq,
			  bool &out,
			  int  &verbose,
			  const char *cvs_id,
			  const char *source)
{
    // Parse the sys_stats command-line (used by both dyn and hdyn versions).

    // Set standard defaults for standalone tools:

    which_lagr = 2;			// print nonlinear Lagrangian zones

    binaries = true;			// print binary output
    B_flag = false;			// force binary evolution (hdyn version)
    calc_e = true;			// compute energy
    n_sq = false;			// don't allow n^2 operations
    out = false;			// write cin to cout
    long_binary_output = true;		// long binary output
    verbose = 2;			// extended output

    extern char *poptarg;
    int c;
    const char *param_string = "0b.Bnel.ov";
					// only the "verbose > 0" options work!

    while ((c = pgetopt(argc, argv, param_string, cvs_id, source)) != -1)
	switch(c) {

	    case '0': break;			// for hdyn compatibility

	    case 'b': {int b = 1;
		      if (poptarg) b = atoi(poptarg);
		      if (b == 0) {
			  binaries = false;
			  long_binary_output = false;
		      } else if (b == 1) {
			  binaries = true;
			  long_binary_output = false;
		      } else {
			  binaries = true;
			  long_binary_output = true;
		      }}
	    case 'B': B_flag = true;
		      break;
	    case 'e': calc_e = false;
		      break;
	    case 'l': if (poptarg)
			  which_lagr = atoi(poptarg);
		      else
			  which_lagr = 2;
		      break;
	    case 'n': n_sq = !n_sq;
		      break;
	    case 'o': out = true;
		      break;
	    case 'v': verbose = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      return false;
	}

    return true;
}


#define TTOL 1.e-12				// arbitrary tolerance

void sys_stats(dyn* b,
	       real energy_cutoff,			// default = 1
	       int  verbose,				// default = 2
	       bool binaries,				// default = true
	       bool long_binary_output,			// default = false
	       int  which_lagr,				// default = 2
	       bool print_time,				// default = false
	       bool compute_energy,			// default = false
	       bool allow_n_sq_ops,			// default = false
	       void (*compute_energies)(dyn*, real,	// default = calculate_
					real&, real&,	//	       energies
					real&, bool),
	       void (*dstar_params)(dyn*),		// default = NULL
	       bool (*print_dstar_stats)(dyn*, bool,	// default = NULL
					 vec, bool))

// The last three arguments are a C-ish way of allowing this dyn function
// to output additional dyn/dstar information when called from kira and
// other hdyn programs.

{
    int p = cerr.precision(STD_PRECISION);

    if (print_time) {

	real time = getrq(b->get_dyn_story(), "t");
	if (time <= -VERY_LARGE_NUMBER)
	    time = b->get_system_time();

	if (time > -VERY_LARGE_NUMBER)
	    cerr << "Time = " << time << endl;
	else
	    cerr << "Time = (undefined)" << endl;

	if (b->get_use_sstar())
	    print_sstar_time_scales(b);
    }

    cerr << "\n  Overall parameters (sys_stats):\n";
    bool mass_spectrum = false;
    print_numbers_and_masses(b, mass_spectrum);

    vec com_pos, com_vel;
    compute_com(b, com_pos, com_vel);

    cerr << "    center of mass position = " << com_pos << endl
	 << "                   velocity = " << com_vel << endl;

    // CM quantities include the root node, but from here on we will work
    // with quantities relative to the root.

    com_pos -= b->get_pos();
    com_vel -= b->get_vel();

    bool heavy_stars = check_for_doubling(b);
    // PRI(4); PRL(heavy_stars);

    real kinetic_energy = 0, potential_energy = 0;

    real r_virial = -1;
    if (compute_energy)				// Need virial radius for tR
	print_relaxation_time(b, potential_energy, kinetic_energy,
			      r_virial, compute_energies);

    // NB:  If set, potential energy here is internal potential energy.

    if (r_virial == -1) {
        if (find_qmatch(b->get_dyn_story(), "energy_time")
	    && twiddles(getrq(b->get_dyn_story(), "energy_time"),
			b->get_system_time(), TTOL)) {

	    // Take energies from the dyn story.

	    if (find_qmatch(b->get_dyn_story(), "kinetic_energy")
		&& find_qmatch(b->get_dyn_story(), "kinetic_energy")) {

		kinetic_energy = getrq(b->get_dyn_story(), "kinetic_energy");
		potential_energy = getrq(b->get_dyn_story(),
					 "potential_energy");
		if (potential_energy < 0)
		    r_virial = -0.5*b->get_mass()*b->get_mass()
					/ potential_energy;
	    }
	}
    }

    // Note: even if allow_n_sq_ops is false, we can still compute kT
    // from the kinetic energy.

    real kT = 0;
    int nd = b->n_daughters();

    const int BIG_ND = 10;

    if (compute_energy)					// recompute energies

        print_energies(b, potential_energy, kinetic_energy, kT,
		       compute_energies);

    else if (b->get_oldest_daughter())

	kT = top_level_kinetic_energy(b) / (1.5*nd);

    real vrms = sqrt(3*nd*kT/b->get_mass());

//  Hmmm.  Better print out (top-level) r_max and v_max, too.

    real r_max = 0, v_max = 0;
    for_all_daughters(dyn, b, bb) {
	real r2 = square(bb->get_pos() - com_pos);
	if (r2 > r_max) r_max = r2;
	real v2 = square(bb->get_vel() - com_vel);
	if (v2 > v_max) v_max = v2;
    }

    r_max = sqrt(r_max);
    v_max = sqrt(v_max);
    PRI(4); PRC(vrms); PRC(kT); PRC(r_max); PRL(v_max);

    vec center = com_pos;
    real rcore = 0, rhalf = 0;

    // "Dominant" mass means > 50% of the total mass of the system.

    dyn* b_dom = NULL;
    if (nd > BIG_ND) b_dom = dominant_mass(b, 0.5);

    if (!b_dom) {

	// No dominant mass.

	bool has_densities = (twiddles(getrq(b->get_dyn_story(), "density_time"),
				       b->get_system_time(), TTOL));

	if (verbose > 0 && (has_densities || allow_n_sq_ops)) {

	    // cerr << "Compute densities" << endl;	// ???

	    // Function compute_core_parameters will look in the dyn story
	    // for densities and use them if they are current.  It will
	    // only do O(N^2) operations if no current densities are found
	    // and allow_n_sq_ops is true.

	    cerr << "\n  Core parameters";
	    if (has_densities)
	      cerr << " (densities from input snapshot)";
	    cerr << ":\n";
	    
	    // Core parameters are always expressed relative to the
	    // (mean) density center.

	    print_core_parameters(b, allow_n_sq_ops, center, rcore);

	    if (rcore > 0 && r_virial > 0) {
		cerr << "\n  Dynamical King parameters:\n";
		print_fitted_king_model(rcore/r_virial, rcore_rvirial);
	    }
	}

	// If core parameters were computed, then center now is the mean
	// density center.  If not, then center is the center of mass.
	// However, this value of center is presently not used, as the
	// vector will be set equal to lagr_pos below.
	//
	// The vector lagr_pos is set when the Lagrangian radii are
	// computed; lagr_pos and the associated Lagrangian radii are
	// used internally by *all* functions which print quantities
	// by radial zone.

	// final variable which = 0 : single stars/CMs
	//			= 1 : binaries
	//			= 2 : light stars (heavy_stars = true)
	//			= 3 : heavy stars (heavy_stars = true)
	//			= 4 : most massive stars (mass_spectrum = true)

	if (verbose > 1 || (verbose > 0 && nd > BIG_ND)) {

	    cerr << endl << "  All single stars/CMs:";
	    rhalf = print_lagrangian_radii(b, which_lagr, verbose>0, 0);
	    PRI(4); PRL(rhalf);

	    if (rhalf > 0) {
	      real density = 1.5*b->get_root()->get_mass()
				/ (4*M_PI*pow(rhalf, 3));
	      putrq(b->get_root()->get_dyn_story(), "kira_rhalf", rhalf);
	      putrq(b->get_root()->get_dyn_story(),
		    "kira_half_density", density);
	      set_new_rhalf();
	    }
	}

	// PRL(heavy_stars);

	if (heavy_stars && verbose > 1) {
	    cerr << endl << "  \"Light\" stars only:";
	    print_lagrangian_radii(b, which_lagr, verbose>0, 2);
	    cerr << endl << "  \"Heavy\" stars only:";
	    print_lagrangian_radii(b, which_lagr, verbose>0, 3);
	}

	if (mass_spectrum) {

	    set_lagr_cutoff_mass(b, 0.90);	// 90th percentile

	    // Print out Lagrangian radii of the most massive stars, as
	    // defined by set_lagr_cutoff_mass() above.

	    cerr << endl << "  Most massive stars (cutoff_mass = "
		 << get_lagr_cutoff_mass() << "):";
	    real rhalf_massive = print_lagrangian_radii(b, which_lagr,
							verbose>0, 4);
	    PRI(4); PRL(rhalf_massive);

	    if (verbose > 1) {
		cerr << endl
		     << "  Mass distribution by Lagrangian zone (singles):"
		     << endl;
		print_numbers_and_masses_by_radial_zone(b, 0);

		cerr << endl
		     << "  Mass distribution by Lagrangian zone (multiples):"
		     << endl;
		print_numbers_and_masses_by_radial_zone(b, 1);
	    }
	}

	if (verbose > 1)
	  print_parameters_for_massive_black_holes(b, kT, center, verbose>1);

	if (verbose > 1) {
	    cerr << "\n  Anisotropy by Lagrangian zone (singles):\n";
	    print_anisotropy_by_radial_zone(b, 0);

	    cerr << "\n  Anisotropy by Lagrangian zone (multiple CMs):\n";
	    print_anisotropy_by_radial_zone(b, 1);

	    if (heavy_stars) {
	        // if (verbose)
		cerr << endl
		     << "  Mass distribution by Lagrangian zone (light stars):"
		     << endl;
		print_numbers_and_masses_by_radial_zone(b, 2);

		// if (verbose)
		cerr << endl
		     << "  Mass distribution by Lagrangian zone (heavy stars):"
		     << endl;
		print_numbers_and_masses_by_radial_zone(b, 3);

		// if (verbose)
		cerr << "\n  Anisotropy by Lagrangian zone (light stars):\n";
		print_anisotropy_by_radial_zone(b, 2);

		// if (verbose)
		cerr << "\n  Anisotropy by Lagrangian zone (heavy stars):\n";
		print_anisotropy_by_radial_zone(b, 3);
	    }
	}

    } else {

	// Node b_dom dominates the motion ("central black hole").

	center = b_dom->get_pos();
	putvq(b->get_dyn_story(), "lagr_pos", center);

	// Assume we have a disk system of some sort, and that the
	// primary plane is x-y.  All functions relating to this
	// option are currently local, but may be made global later,
	// if appropriate.

	cerr << endl << "  Motion dominated by node "
	     << b_dom->format_label() << endl
	     << "      (mass = " << b_dom->get_mass()
	     << ",  pos = " << center << ")" << endl;
	cerr << "  Taking x-y plane as plane of symmetry." << endl;

	print_dominated_lagrangian_radii(b, b_dom);
	print_dominated_velocity_dispersions(b, b_dom);
    }

    bool sstar = b->get_use_sstar();

    if (print_dstar_stats != NULL)
	sstar = !print_dstar_stats(b, mass_spectrum, center, verbose>0);

    if (sstar)
	  sstar_stats(b, mass_spectrum, center, verbose>0);

    if (binaries) {

	// Use the same center as used by the Lagrangian radii functions.

	if (find_qmatch(b->get_dyn_story(), "lagr_pos"))
	    center = getvq(b->get_dyn_story(), "lagr_pos");
	else
	    if (verbose > 1) warning("sys_stats: lagr_pos not found");

	set_kepler_tolerance(2);	// avold termination on error

	if (verbose > 0) {
	    cerr << endl;
	    cerr << "  Binaries/multiples";
	    if (!long_binary_output) cerr << " (short output)";
	    cerr << ":" << endl;
	}

	print_binaries(b, kT, center, rcore, rhalf,
		       verbose, long_binary_output, dstar_params);

	if (allow_n_sq_ops) {

	    if (verbose > 0) {
		cerr << "\n  Other bound pairs with ";
		if (kT > 0)
		    cerr << "|E| > " << energy_cutoff << " kT:\n";
		else
		    cerr << "|E/mu| > " << energy_cutoff << ":\n";
	    }
	    search_for_binaries(b, energy_cutoff, kT,
				center, verbose>0, long_binary_output);

	} else {

	    // Do an NN search for binaries (hdyn only).

	    int which = 0;
	    bool found = false;
	    for_all_daughters(dyn, b, bb)
		found |= bb->nn_stats(energy_cutoff, kT,   // Dummy function for
				      center, verbose>0,   // dyn; get NN binary
				      long_binary_output,  // info for hdyn
				      which++);
	    if (!found)
		cerr << "    (none)\n";
	}
    }

    cerr.precision(p);

    // Finally, refine estimates of the cluster mass (moved here from hdyn...)

    refine_cluster_mass(b, 1);
    cerr << endl;
}

#else

main(int argc, char **argv)
{
    check_help();

    int  verbose;
    bool binaries, long_binary_output, B_flag, out, n_sq, calc_e;
    int which_lagr;

    if (!parse_sys_stats_main(argc, argv,
			      which_lagr,
			      binaries, long_binary_output, B_flag,
			      calc_e, n_sq, out, verbose,
			      "$Revision: 1.46 $", _SRC_)) {
	get_help();
	exit(1);
    }

    // For dyn version, calc_e requires n_sq.

    // if (!n_sq) calc_e = false;

    // Loop over input until no more data remain.

    dyn *b;
    int i = 0;

    while (b = get_dyn())	{ // NB:  get_xxx() reads NaN as legal garbage...

	check_addstar(b);
	check_set_external(b, true);	// true ==> verbose output
	cerr << endl;

	if (i++ > 0) cerr << endl;
	sys_stats(b, 0.5, verbose, binaries, long_binary_output,
		  which_lagr, true, calc_e, n_sq);

	if (out) put_dyn(b);
	rmtree(b);
    }
}

#endif
