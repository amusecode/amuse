
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute the 1-dimensional radial velocity dispersion profile of
//// an N-body system.  The function fills in the vdisp array
//// corresponding to the provided radius array.  The tool prints out
//// radius and velocity dispersion in a form suitable for plotting.
//// If a current density center is found in the root dyn story it is
//// used as the system center for the density computation.  Otherwise,
//// if a valid center of mass is found, it is used.  If neither center
//// is found, the geometric center is used.
////
//// Usage: radial_vdisp [OPTIONS] < input > output
////
//// Options:     
////		  -n    specify number of radial bins [100]
////              -r    specify maximum radius [take from data]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//-----------------------------------------------------------------------------
//   version 1:  Apr 2003   Steve McMillan
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

// Prototype function.  Could also be implemented using get_radial_profile()
// with a target function of v*v.

// Sorting code is taken almost verbatim from lagrad.C and radial_density.C.

typedef struct {
    real r_sq;
    real mass;
    real v_sq;
} rmv, *rmv_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rmv_ptr)pi)->r_sq > ((rmv_ptr)pj)->r_sq)
        return +1;
    else if (((rmv_ptr)pi)->r_sq < ((rmv_ptr)pj)->r_sq)
        return -1;
    else
        return 0;
}

int get_radial_vdisp(dyn *b, vec cpos, vec cvel,
		     int n_zones, real r[], real v2[])
{
    if (n_zones < 2) return 1;

    // Set up an array of (r_sq, mass, vel) triples.

    int n = b->n_daughters();	// (NB implicit loop through the entire system)

    // Would be possible to determine n and set up the array simultaneously
    // using malloc and realloc.  Not so easy with new...
    
    rmv_ptr table = new rmv[n];

    if (!table) {
	cerr << "get_radial_vdisp: insufficient memory for table"
	     << endl;
	return 1;
    }

    int i = 0;
    for_all_daughters(dyn, b, bi) {
	table[i].r_sq = square(bi->get_pos() - cpos);
	table[i].mass = bi->get_mass();
	table[i].v_sq = square(bi->get_vel() - cvel);
	i++;
    }

    // Sort the array by radius (may repeat work done elsewhere...).

    qsort((void *)table, (size_t)i, sizeof(rmv), compare_radii);

    // Initialize the target arrays (mass and mv^2).

    real mass[n_zones];
    int j;
    for (int j = 0; j < n_zones; j++) mass[j] = v2[j] = 0;

    // Bin the (ordered) data.

    j = 0;
    real rj2 = r[j]*r[j];

    for (i = 0; i < n; i++) {
	real ri_sq = table[i].r_sq;
	while (ri_sq > rj2) {
	    j++;
	    if (j >= n_zones) break;
	    rj2 = r[j]*r[j];
	}
	if (j >= n_zones) break;
	mass[j] += table[i].mass;
	v2[j] += table[i].mass * table[i].v_sq;
    }

    // Convert from mv^2 to <v^2>.

    for (j = 0; j < n_zones; j++)
	if (mass[j] > 0) v2[j] /= mass[j];

    delete[] table;
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  get_radial_vdisp() as a tool.
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
		    "$Revision: 1.11 $", _SRC_)) != -1)
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

    dyn *b;

    while (b = get_dyn()) {

	real r[n_zones], v2[n_zones];

	// Use the density center if known and up to date (preferred).
	// Otherwise, use modified center of mass, if known and up to date.

	vec cpos, cvel;
	get_std_center(b, cpos, cvel);

	cpos -= b->get_pos();		// std_center quantities include
	cvel -= b->get_vel();		// the root node

	// See if we need to compute r_max.

	if (r_max <= 0) {
	    r_max = 0;
	    for_all_daughters(dyn, b, bb) {
		real r2 = square(bb->get_pos()-cpos);
		if (r2 > r_max) r_max = r2;
	    }
	    r_max = sqrt(r_max);
	}

	// Set up the radial array.

	for (int i = 0; i < n_zones; i++)
	    r[i] = (i+1) * r_max / n_zones;	// r[i] = outer edge of zone i

	// Compute and print the velocity dispersion array.

	get_radial_vdisp(b, cpos, cvel, n_zones, r, v2);

	real r0 = 0;
	for (int i = 0; i < n_zones; i++) {
	    real r1 = r[i];
	    cout << i << " " << (r0+r1)/2 << " " << v2[i] << endl;
	    r0 = r1;
	}

	rmtree(b);
    }
}

#endif
