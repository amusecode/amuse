
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute local densities around particles in the input N-body
//// system, based on k-th nearest neighbors. Save results in the 
//// particle dyn stories.
////
//// Usage: compute_density [OPTIONS] < input > output
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -k    specify k [12]
////              -v    verbose mode [off]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  May 1989   Piet Hut
//    version 2:  Nov 1994   Piet Hut
//.............................................................................
//  non-local function: 
//    compute_density
//.............................................................................
//  see also: density_center.C
//.............................................................................


// This program can be very slow, as it runs on the front end.
// It should be made to use the GRAPE if available...

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//
// compute_density  --  Compute the density for one or all particles.
//               note:
//                    the neighbor_dist_sq[q] array will contain the distances
//                    to the q-th nearest neighbor (in the form of the squares 
//                    thereof); therefore  neighbor_dist_sq[0] = 0.0  since
//                    this indicates the distance of a particle to itself.
//                    Including this trivial case simplifies the algorithm
//                    below by making it more homogeneous.
//-----------------------------------------------------------------------------
//              Litt.:
//                    Stefano Casertano and Piet Hut: Astroph.J. 298,80 (1985).
//                    To use their recipe, k >= 2 is required.
//-----------------------------------------------------------------------------

local void write_density(dyn *d, int k, real density)
{
    putrq(d->get_dyn_story(), "density_time", d->get_system_time(),
	  HIGH_PRECISION);
    putiq(d->get_dyn_story(), "density_k_level", k);
    putrq(d->get_dyn_story(), "density", density);
}

static int nwarn = 0;

void  compute_density(dyn * b,	      // pointer to N-body system or node
		      int k,	      // use kth nearest neighbor [12]
		      bool use_mass,  // use mass density [no]
		      dyn ** list,    // list of neighbor nodes to work from [0]
		      int n_list)     // number of nodes on the list [0]

// If list = NULL [default] or n_list = 0, use the entire system.

{
    int  q, r;
    real *neighbor_dist_sq;
    real *neighbor_mass;
    real  delr_sq;

    // Special cases and error checks.


    if (k == 0) {

	// Just write zero density for b and quit.

	write_density(b, 0, 0);
	return;
    }
    
    if (k <= 1) {
        if (nwarn++ < 5)
	    cerr << "compute_density: k = " << k << " <= 1" << endl;
	return;
    }

    if (list == NULL)
        n_list = 0;
    else if (n_list <= 0) {
        if (nwarn++ < 5)
	    cerr << "compute_density: n_list = " << n_list << endl;
	return;
    }

    int nsys = b->n_leaves() - 1;
    if (n_list > 0)
        nsys = n_list;

    if (k > nsys) {			// no k-th nearest neighbor exists
        if (nwarn++ < 5)
	    cerr << "compute_density: k = " << k << " >= nsys = "
		 << nsys << endl;
	return;
    }

    neighbor_dist_sq = new real[k+1];
    neighbor_mass = new real[k+1];

    // Set first body d for which density is to be computed.

    // for_all_leaves(dyn, b, d) {	// The density is now defined
    // for_all_daughters(dyn, b, d) {	// ONLY for top-level nodes.

    dyn* d;

    if (list)
        d = b;
    else
        d = b->get_oldest_daughter();	// b is root in this case

    while (d) {

	neighbor_dist_sq[0] = 0.0;
	neighbor_mass[0] = d->get_mass();
	for (q = 1; q <= k; q++) {
	    neighbor_dist_sq[q] = VERY_LARGE_NUMBER;
	    neighbor_mass[q] = 0;
	}

	// Search bodies dd to find k-th nearest neighbor.
	// May be slow if the original list is long -- qsort should
	// perhaps be used.

	// for_all_leaves(dyn, b, dd) {
	// for_all_daughters(dyn, b, dd) {

	dyn* dd;

	if (list)
	    dd = list[n_list-1];
	else
	    dd = b->get_oldest_daughter();	// b is root

	while (dd) {

	    if (d != dd) {

		delr_sq = square(something_relative_to_root(d, &dyn::get_pos)
			      - something_relative_to_root(dd, &dyn::get_pos));

		if (delr_sq < neighbor_dist_sq[k]) {

		    // Place dd on d's neighbor list.

		    for (q = k-1; q >= 0; q--) {
			if (delr_sq > neighbor_dist_sq[q]) {
			    for (r = k; r > q+1; r--) {
				neighbor_dist_sq[r] = neighbor_dist_sq[r-1];
				neighbor_mass[r] = neighbor_mass[r-1];
			    }
			    neighbor_dist_sq[q+1] = delr_sq;
			    neighbor_mass[q+1] = dd->get_mass();

			    break;
			}
		    }
		}
	    }

	    // Get next dd.

	    if (list) {
		if (--n_list < 0)
		    dd = NULL;
		else
		    dd = list[n_list];
	    } else
		dd = dd->get_younger_sister();
	}
	    
        real density =  0;

	if (neighbor_dist_sq[k] > 0) {
	    while (k > 0 && neighbor_dist_sq[k] >= VERY_LARGE_NUMBER)
		k--;
	    if (k > 1 && neighbor_dist_sq[k] > 0
		      && neighbor_dist_sq[k] < VERY_LARGE_NUMBER) {
		real volume = (4.0/3.0) * PI * pow(neighbor_dist_sq[k], 1.5);
		if (volume > 0) {

//		    density = (k - 1) / volume;		// ApJ, 298, 80 (1985)

		    real mass = 0;   
		    for(int m = 1; m < k; m++)		// exclude self and
			mass += neighbor_mass[m]; 	// outermost star
		    density = mass/volume;

		}
	    }
	}

	// Store the density in d's dyn story.

	write_density(d, k, density);

	// Get next d.

	if (list)
	    d = NULL;
	else
	    d = d->get_younger_sister();
    }

    // Timestamp the density computation if the entire system was involved.

    if (!list) {
	dyn* root = b->get_root();
	putrq(root->get_dyn_story(), "density_time", root->get_system_time(),
	      HIGH_PRECISION);
    }

    // Clean up.

    delete neighbor_dist_sq;
    delete neighbor_mass;
}

#else

main(int argc, char ** argv)
{
    int  k = 12;                // default

    char *comment;
    bool c_flag = false;        // if true: a comment given on command line
    bool verbose = false;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:k:v";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.15 $", _SRC_)) != -1)
	switch(c) {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'k': k = atoi(poptarg);
		      break;
	    case 'v': verbose = !verbose;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    // Loop over input until no more data remain.

    dyn *b;
    int i = 0;
    while (b = get_dyn()) {

	if (verbose) cerr << "snap #" << ++i
	                  << ", time = " << b->get_system_time() << endl;

        if (c_flag) b->log_comment(comment);
        b->log_history(argc, argv);

        compute_density(b, k);

        put_dyn(b);
	rmtree(b);
    }
}

#endif

// endof: compute_density.C
