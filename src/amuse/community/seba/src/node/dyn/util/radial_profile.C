
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute the specified 1-dimensional radial profile of an N-body system.
//// The function fills in a profile array corresponding to the provided
//// radius array, then prints out radius and profile data in a form
//// suitable for plotting.  If a current density center is found in the root
//// dyn story it is used as the system center for the density computation.
//// Otherwise, if a valid center of mass is found, it is used.  If neither
//// center is found, the geometric center is used.
////
//// Usage: radial_profile [OPTIONS] < input > output
////
//// Options:     
////		  -n    specify number of radial bins [100]
////              -o    choose the quantity computed [0]
////                    0 ==> density of all stars,
////                    1 ==> density of stars having mass >= average,
////                    2 ==> velocity dispersion,
////                    3 ==> radial velocity dispersion,
////                    4 ==> 2, 3, and radial velocity dispersion too
////              -r    specify maximum radius [4]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//-----------------------------------------------------------------------------
//   version 1:  Jun 2003   Steve McMillan
//-----------------------------------------------------------------------------

#include "dyn.h"
#include <vector>
#include <algorithm>

#ifndef TOOLBOX

// New sorting code using STL.

class dr_pair {
    public:
	dyn *b;
	real r_sq;
	dr_pair(dyn *bb = NULL, real rr = 0) {b = bb; r_sq = rr;}
};

bool operator < (const dr_pair& x, const dr_pair& y) {return x.r_sq < y.r_sq;}

// get_density_profile:  Return the density profile of a system, subject to
// an optional user-provided selection function defining the particles to be
// included in the calculation.

int get_density_profile(dyn *b, vec cpos,
			int n_zones, real r[], real rho[],
			bool (*select)(dyn*))		    // default = NULL
{
    if (n_zones <= 0) return 1;

    // Initialize the density array.

    int j;
    for (j = 0; j < n_zones; j++) rho[j] = 0;

    // Set up an array of (dyn, radius) pairs, including a
    // possible selection function.

    vector<dr_pair> v;

    for_all_daughters(dyn, b, bb)
	if (!select || select(bb)) {
//	    PRC(bb->get_mass()); PRL(bb->format_label());
	    v.push_back(dr_pair(bb, square(bb->get_pos() - cpos)));
	}

    // Sort the array by radius (may repeat work done elsewhere...).

    sort(v.begin(), v.end());

    // Bin the ordered data.

    j = 0;
    real rj2 = r[j]*r[j];

    vector<dr_pair>::iterator vi;
    for (vi = v.begin(); vi < v.end(); vi++) {
	while (vi->r_sq > rj2) {
	    if (++j >= n_zones) break;
	    rj2 = r[j]*r[j];
	}
	if (j >= n_zones) break;
	rho[j] += vi->b->get_mass();
    }

    // Convert from mass to density.

    real v0 = 0;	// assume that the first zone extends in to r = 0
    for (j = 0; j < n_zones; j++) {
	real v1 = pow(r[j], 3);
	rho[j] /= (4*M_PI/3) * (v1 - v0);	// dM --> dM/dV
	v0 = v1;
    }
}

// get_profile:  Return the mass-weighted radial profile of some
// user-specified quantity Q(dyn*).

int get_profile(dyn *b, vec cpos,
		int n_zones, real r[],
		real q[], real (*Q)(dyn*))
{
    static dyn *bprev = NULL;
    static vec cprev = 0;
    static vector<dr_pair> v;

    if (!b) {
	bprev = NULL;		// reset
	return 1;
    }
    if (n_zones <= 0 || !Q) return 1;

    real mass[n_zones];

    // Initialize the mass and output arrays.

    int j;
    for (j = 0; j < n_zones; j++) q[j] = mass[j] = 0;

    // Set up an array of (dyn, radius) pairs, and sort it,
    // if necessary.  Check to see if the work has already
    // been done for this system and center.

    if (b != bprev || cpos != cprev) {

	for_all_daughters(dyn, b, bb)
	    v.push_back(dr_pair(bb, square(bb->get_pos() - cpos)));

	// Sort the array by radius (may repeat work done elsewhere...).

	sort(v.begin(), v.end());

	bprev = b;
	cprev = cpos;
    }

    // Bin the ordered data.

    j = 0;
    real rj2 = r[j]*r[j];

    static vector<dr_pair>::iterator vi;
    for (vi = v.begin(); vi < v.end(); vi++) {
	while (vi->r_sq > rj2) {
	    if (++j >= n_zones) break;
	    rj2 = r[j]*r[j];
	}
	if (j >= n_zones) break;
	mass[j] += vi->b->get_mass();
	q[j] += vi->b->get_mass() * Q(vi->b);
    }

    // Complete the calculation.

    for (j = 0; j < n_zones; j++) if (mass[j] > 0) q[j] /= mass[j];
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  get_density_profile() or get_profile() as a tool.
//-----------------------------------------------------------------------------

#define N_DEFAULT 25
#define R_DEFAULT 5

// Global data (but only within the main program...)!!

real m_bar = 0;
vec cpos, cvel;

// Sample selector function.

bool select(dyn *b)
{
    return (b->get_mass() >= m_bar);			// note mass() wrapper
}

// Sample "quantity" functions.

real Q1(dyn *b)					// v^2
{
    return speed_sq(b, cvel);
}

real Q2(dyn *b)					// vt^2
{
    return v_trans_sq(b, cpos, cvel);
}

real Q3(dyn *b)					// vr^2
{
    return v_rad_sq(b, cpos, cvel);
}

typedef real (*qfn)(dyn*);

main(int argc, char ** argv)
{
    int n_zones = N_DEFAULT;
    int option = 0;
    real r_max = R_DEFAULT;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "n:o:r:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {

	    case 'n': n_zones = atoi(poptarg);
		      break;

	    case 'o': option = atoi(poptarg);
		      break;

	    case 'r': r_max = atof(poptarg);
		      break;

            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            

    if (n_zones <= 1) n_zones = N_DEFAULT;
    cerr << "radial_profile: "; PRC(option); PRC(r_max); PRL(n_zones);

    dyn *b;
    while (b = get_dyn()) {

	real r[n_zones], q[n_zones], q2[n_zones], q3[n_zones];

	// Use the density center if known and up to date (preferred).
	// Otherwise, use modified center of mass, if known and up to date.

	get_std_center(b, cpos, cvel);

	cpos -= b->get_pos();			// std_center quantities
	cvel -= b->get_vel();			// include the root node

	// Compute m_bar and r_max, if necessary.

	int n = 0;
	m_bar = 0;
	bool need_r_max = (r_max <= 0);
	if (need_r_max) r_max = 0;

	for_all_daughters(dyn, b, bb) {
	    n++;
	    m_bar += bb->get_mass();
	    if (need_r_max) {
		real r_sq = square(bb->get_pos()-cpos);
		if (r_sq > r_max) r_max = r_sq;
	    }
	}
	m_bar /= n;
	if (need_r_max) r_max = sqrt(r_max);

	if (option <= 1) PRL(m_bar); 

	// Set up the radial array.

	for (int i = 0; i < n_zones; i++)
	    r[i] = (i+1) * r_max / n_zones;	// r[i] = outer edge of zone i

	// Compute and print the density or profile array.

	if (option == 0)
	    get_density_profile(b, cpos, n_zones, r, q);
	else if (option == 1)
	    get_density_profile(b, cpos, n_zones, r, q, select);
	else if (option == 2)
	    get_profile(b, cpos, n_zones, r, q, Q1);
	else if (option == 3)
	    get_profile(b, cpos, n_zones, r, q, Q2);
	else if (option == 4) {
	    get_profile(b, cpos, n_zones, r,  q, Q1);
	    get_profile(b, cpos, n_zones, r, q2, Q2);
	    get_profile(b, cpos, n_zones, r, q3, Q3);
	}

	real r0 = 0;
	for (int i = 0; i < n_zones; i++) {
	    real r1 = r[i];
	    cout << i << " " << (r0+r1)/2 << " " << q[i];
	    if (option == 4)
		cout << " " << q2[i] << " " << q3[i];
	    cout << endl;
	    r0 = r1;
	}

	rmtree(b);
    }
}

#endif
