
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Group-finding tool based on the HOP algorithm for N-body systems
//// (Eisenstein, D.J. & Hut, P., 1998, ApJ 498, 137).  This version
//// partitions particles according to HOP, but doesn't merge groups
//// or apply a density threshhold.
////
//// The algorithm requires that particles have densities associated
//// with them.  If none are found, they will be computed internally
//// (*without* using GRAPE!): use get_densities on a GRAPE-enabled
//// machine for a faster computation.
////            
//// Usage:    hop < input_snapshot > list_of_clump_memberships
////
//// Options:     -i    number the particles sequentially [don't number]
////              -n    number of neighbors searched for density extremum
////              -v    write snapshot at end [no]
////
//// Written by Simon Portegies Zwart.

//   Version 1.0:       Simon Portegies Zwart,          Hamilton, Aug 2004
//   Adapted to starlab:
//   Version 1.1:       Simon Portegies Zwart,          Haarlem, Mar 2005
//   Modified by	Steve McMillan,			Philadelphia, Sep 2006

#ifndef TOOLBOX

#else

#include "hop.h"

dyn* hop::densest_nth_nn(dyn *b)
{
    // Return a pointer to the node with the greatest density
    // among the nn_search nearest neighbors of b.

    nearest_neighbor *nn = new nearest_neighbor[nn_search];

    for_all_daughters(dyn, b->get_root(), bb) {
	real sep2 = square(bb->get_pos() - b->get_pos());
	for (int i = 0; i < nn_search; i++) {
	    if (sep2 < nn[i].d_nn_sq) {
		for (int j = nn_search-1; j > i; j--) {
		    nn[j] = nn[j-1];
		}
		nn[i].nn = bb;
		nn[i].d_nn_sq = sep2;
		break;
	    }
	}
    }
    real d, dmax = -VERY_LARGE_NUMBER;
    int imax;
    for (int i = 0; i < nn_search; i++) {
	d = getrq(nn[i].nn->get_dyn_story(), "density");
	if (d > dmax) {
	    dmax = d;
	    imax = i;
	}
    }

    return nn[imax].nn;
}

void hop::add_cluster_center(dyn* bc, dyn *bi)
{
    // Create/update the cluster labeled by the center bc, found
    // following the search starting from node bi.

    vector<cluster>::iterator ic; 
    bool already_known = false;

    for (ic = cl.begin(); ic < cl.end(); ic++)
	if (ic->get_h() == bc) {
	    already_known = true;
	    break;
	}

    if (already_known)
	ic->increment(bi);
    else {
	cluster ccl(bc, bi);
	cl.push_back(ccl);
    }
}

void hop::find_clump_center(dyn *bi)
{
    // Follow the density gradient to the densest neighbor of bi;
    // add it to the list of cluster centers if not already on it.

    // cerr << "Finding clump center for i = " <<  bi->get_index() << endl;

    dyn *onnd, *nnd = bi;
    do {
	onnd = nnd;
	nnd = densest_nth_nn(onnd);
    }
    while (onnd != nnd);

    // cerr << "clump center id = " << nnd->get_index() << endl;

    // Put the clump center on the list of clumps, update clump
    // properties, and associate bi with the clump.

    add_cluster_center(nnd, bi);
    putiq(bi->get_dyn_story(), "hop_clump_center_id", nnd->get_index());
}

void hop::find_primary_cluster(dyn *b)
{
    if(!nn_search)
	nn_search = int(3*sqrt((real)b->n_daughters()));

    for_all_daughters(dyn, b, bi) 
	find_clump_center(bi);

    // Need an extra routine to merge clumps, as in Eisenstein & Hut.
}

int main(int argc, char ** argv)
{
    int nth_neighbor = -1;
    bool use_nsqrt = true;

    extern char *poptarg;
    int c; 
    const char *param_string = "in:v";

    check_help();

    bool i_flag = false;
    bool v_flag = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'i': i_flag = true;
		      break;
	    case 'n': nth_neighbor = atoi(poptarg);
	              use_nsqrt = false;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (use_nsqrt)
	cerr << "hop: no nn-th neighbor specified, using sqrt(3n)" << endl;

    bool first = true;
    int i = 0;

    dyn* b;
    hop *h;

    while((b = get_dyn())) {

	// cerr << "hop: read snapshot #" << i++ << endl;

	if(use_nsqrt)
	    nth_neighbor = int(3*sqrt((real)b->n_daughters()));

	// See if we need densities (use 12th nearest neighbor of we do).

	bool compute_dens = true;
	story *s = b->get_oldest_daughter()->get_dyn_story();
	if (find_qmatch(s, "density_time"))
	    if (getrq(s, "density_time") == b->get_system_time())
		compute_dens = false;

	if (compute_dens) {
	    cerr << "hop: computing densities..." << endl;
	    compute_density(b, 12);
	}

	h = new hop();
	h->set_nn_search(nth_neighbor);
      
	h->find_primary_cluster(b);	// decompose system into clums
	h->put();			// print out info on the clumps


	//-------------------------------------------------------------

	// *** Detailed info on the members of each clump is available
	// *** as follows:

	vector<cluster> cl = h->get_cl();	// the list of clusters
	vector<cluster>::iterator i;

	for (i = cl.begin(); i != cl.end(); i++) {

	    cluster c = *i;			// the current cluster; has
						// mass, nstar, pos, vel, and
						// a list of nodes in "clump"

	    vector<dynptr> clump = c.get_clump();
	    vector<dynptr>::iterator ci;

	    for (ci = clump.begin(); ci != clump.end(); ci++) {

		dyn *bi = *ci;			// current clump member

//		PRC(ci-clump.begin());		// counter
//		PRL(bi->get_index());		// (etc.)
	    }
	}

	//-------------------------------------------------------------

	delete h;

	if (v_flag)
	    put_dyn(b);
    }
}

#endif
/* end of: hop.C */
