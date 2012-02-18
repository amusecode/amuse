
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute the 1-dimensional radial density profile of an N-body system.
//// The function fills in the density array corresponding to the provided
//// radius array.  The tool prints out radius and density in a form
//// suitable for plotting.  If a current density center is found in the
//// root dyn story it is used as the system center for the density
//// computation.  Otherwise, if a valid center of mass is found, it is used.
//// If neither center is found, the geometric center is used.
////
//// Usage: radial_numdensity [OPTIONS] < input > output
////
//// Options:     
////		  -n    specify number of radial bins [100]
////              -r    specify maximum radius [take from data]
////
//// Written by Steve McMillan and Ernest Mamikonyan.
////
//// Report bugs to starlab@sns.ias.edu.

//-----------------------------------------------------------------------------
//   version 1:  Apr 2003   Steve McMillan
//   version 2:  Apr 2004   Ernest Mamikonyan
//-----------------------------------------------------------------------------

#include <vector>
#include <algorithm>
#include "dyn.h"

#ifndef TOOLBOX

// Prototype function.  Could also be implemented using get_density_profile()
// with a mass selection function that accepts all stars.


// Sorting code is taken almost verbatim from lagrad.C.

typedef struct {
    real r_sq;
    real mass;
} rm_pair;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

// essentially, bool rm_pair::operator<(const rm_pair&)
local bool compare_radii(const rm_pair& rm1, const rm_pair& rm2)
{
    return rm1.r_sq < rm2.r_sq;
}

vector<real>& get_radial_numdensities(dyn *b, vec cpos, vector<real>& r,
				   bool (*filter)(dyn*))
{
    //if (r.size() < 2) return 0;

    // Set up an array of (r_sq, mass) pairs.
    vector<rm_pair>* table = new vector<rm_pair>;

    for_all_daughters(dyn, b, bi)
        if (!filter || filter(bi))
	    table->push_back((rm_pair)
			     {square(bi->get_pos()-cpos), bi->get_mass()});

    // Sort the array by radius (may repeat work done elsewhere...).
    sort(table->begin(), table->end(), compare_radii);

    // Bin the (ordered) data.
    vector<real>& rho = *(new vector<real>(r.size()));
    real rj2 = r[0]*r[0];

    for (int i = 0, j = 0; i < table->size(); i++) {
	real ri_sq = (*table)[i].r_sq;
	while (ri_sq > rj2) {
	    if (++j >= r.size()) break;
	    rj2 = r[j]*r[j];
	}
	if (j >= r.size()) break;
	++rho[j];
    }
    delete table;

    real v0 = 0;	// assume that the first zone extends in to r = 0
    for (int j = 0; j < r.size(); j++) {
	real v1 = pow(r[j], 3);
	rho[j] /= (4*M_PI/3) * (v1 - v0);	// dM --> dM/dV
	v0 = v1;
    }
    return rho;
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  get_radial_densities() as a tool.
//-----------------------------------------------------------------------------

#define N_DEFAULT 100

main(int argc, char ** argv)
{
    int n_zones = N_DEFAULT;
    real r_max = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "n:r:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c) {

	    case 'n': n_zones = atoi(poptarg);
		      break;

	    case 'r': r_max = atof(poptarg);
		      break;

            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            

    if (n_zones <= 1) n_zones = N_DEFAULT;
    PRC(r_max); PRL(n_zones);

    vector<real> r(n_zones);


    while (dyn* b = get_dyn()) {

	// Use the density center if known and up to date (preferred).
	// Otherwise, use modified center of mass, if known and up to date.

	vec cpos, cvel;
	get_std_center(b, cpos, cvel);

	cpos -= b->get_pos();			// std_center quantities
	cvel -= b->get_vel();			// include the root node

	// See if we need to compute r_max.

	if (r_max <= 0) {
	    r_max = 0;
	    for_all_daughters(dyn, b, bb) {
		real r2 = square(bb->get_pos()-cpos);
		if (r2 > r_max) r_max = r2;
	    }
	    r_max = sqrt(r_max);
	}

	// Set up a linear radial array.

	for (int i = 0; i < n_zones; i++)
	    r[i] = (i+1) * r_max / n_zones;	// r[i] = outer edge of zone i

	// Compute and print the density array.

	vector<real>& rho = get_radial_numdensities(b, cpos, r);

	real r0 = 0;
	for (int i = 0; i < n_zones; i++) {
	    real r1 = r[i];
	    cout << i << " " << (r0+r1)/2 << " " << rho[i] << endl;
	    r0 = r1;
	}

	rmtree(b);
	delete &rho;
    }
}

#endif
