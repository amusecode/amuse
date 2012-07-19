
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// plot_stars:  plot part of an N-body system.

//.............................................................................
//    version 1:  Steve, Jun 2000, adapted from "starplot"
//.............................................................................
//  Plot positions of some bodies in an N-body system, in a primitive way:
//  the screen is divided in character positions, and the representation is:
//  every character stands for one body; every * for more than one body.
//  Very primitive, but independent of terminal or plotting software.
//.............................................................................

#include "dyn.h"

#ifndef TOOLBOX

local real find_neighbors(dyn *bi, dyn **list, int n)	// note: this n is n+1!
{
    dyn *root = bi->get_root();

    real *dr2 = new real[n];
    for (int i = 0; i < n; i++) {
	list[i] = NULL;
	dr2[i] = VERY_LARGE_NUMBER;
    }

    for_all_daughters(dyn, root, bj) {

	real sep2 = square(bj->get_pos() - bi->get_pos());
	int j;

	for (j = n-1; j >= 0; j--)
	    if (sep2 > dr2[j]) break;

	// bj should be in location j+1

	if (j < n-1) {
	    for (int k = n-1; k > j+1; k--) {
		list[k] = list[k-1];
		dr2[k] = dr2[k-1];
	    }
	    list[j+1] = bj;
	    dr2[j+1] = sep2;
	}
    }

    for (int i = n-1; i >= 0; i--)
	if (list[i]) return sqrt(dr2[i]);

    return 0;
}

#define  HBINS  41			// number of horizontal bins
#define  VBINS  21			// number of vertical bins

void plot_stars(dyn * bi,
		int n,			// default = 5
		int k)			// default = 3: plot (x, y)

// Draw the n stars closest to bi.

{
    if (n <= 0) return;
    bi = bi->get_top_level_node();	// deal only with top-level nodes

    // Make a list of the n (actual, not projected) nearest neighbors of bi.
    // Convenient to keep bi at location 0.

    dyn **list = new dynptr[n+1];
    real scale = find_neighbors(bi, list, n+1);
    if (scale <= 0) return;

    // Clear the display and make a box.

    char disp[HBINS][VBINS];		// pixels on the display

    for (int j = 0; j < VBINS; j++)
	for (int i = 0; i < HBINS; i++)
	    disp[i][j] = ' ';

    for (int j = 0; j < VBINS; j++)
	disp[0][j] = disp[HBINS-1][j] = '|';

    for (int i = 0; i < HBINS; i++)
	disp[i][0] = disp[i][VBINS-1] = '-';

    disp[0][0] = disp[0][VBINS-1]
	= disp[HBINS-1][0] = disp[HBINS-1][VBINS-1] = '+';

    // Add the stars.

    int  kx, ky;			// projected coordinates

    switch (k) {
	case 1: kx = 1; ky = 2; break;	// project along x
	case 2: kx = 2; ky = 0; break;	// project along y
	default:
	case 3: kx = 0; ky = 1; break;	// project along z
    }

    // Refine scale before proceeding.

    scale = 0;
    for (int m = 0; m <= n; m++) {
	if (list[m]) {
	    real x = list[m]->get_pos()[kx] - list[0]->get_pos()[kx];
	    real y = list[m]->get_pos()[ky] - list[0]->get_pos()[ky];
	    if (abs(x) > scale) scale = abs(x);
	    if (abs(y) > scale) scale = abs(y);
	}
    }

    real scale2 = 256;
    while (scale2 > scale) scale2 /= sqrt(2.0);
    scale2 *= sqrt(2.0);
    scale = scale2;

    // Create the plot.

    for (int m = 0; m <= n; m++) {
	if (list[m]) {

	    // Center of mass.

	    real x = list[m]->get_pos()[kx] - list[0]->get_pos()[kx];
	    real y = list[m]->get_pos()[ky] - list[0]->get_pos()[ky];
	    int i = (int)(0.5*(0.999999 + x/scale) * HBINS);
	    int j = (int)(0.5*(0.999999 + y/scale) * VBINS);
	    if (disp[i][j] == ' '
		|| i == 0 || i == HBINS-1
		|| j == 0 || j == VBINS-1) {
		if (m < 10)
		    disp[i][j] = '0' + m;
		else
		    disp[i][j] = 'a' + m - 10;
	    } else
		disp[i][j] = '*';

	    // Add next-level components, if any.

	    dyn *od = list[m]->get_oldest_daughter();
	    if (od) {
		x += od->get_pos()[kx];
		y += od->get_pos()[ky];
		int icm = i, jcm = j;
		i = (int)(0.5*(0.999999 + x/scale) * HBINS);
		j = (int)(0.5*(0.999999 + y/scale) * VBINS);
		if (disp[i][j] == ' '
		     || (disp[i][j] != '*' && i == icm && j == jcm))
		    disp[i][j] = 'X';
		else
		    disp[i][j] = '*';

		dyn *yd = od->get_younger_sister();
		x += yd->get_pos()[kx] - od->get_pos()[kx];
		y += yd->get_pos()[ky] - od->get_pos()[ky];
		i = (int)(0.5*(0.999999 + x/scale) * HBINS);
		j = (int)(0.5*(0.999999 + y/scale) * VBINS);
		if (disp[i][j] == ' ')
		    disp[i][j] = 'x';
		else
		    disp[i][j] = '*';
	    }
	}
    }

    // Print the results.

    int off = 10;
    cerr << endl;
    PRI(off+1); cerr << "system_time = " << bi->get_system_time() << endl;
    PRI(off+1); cerr << "plot scale  = " << 2*scale << endl << endl;

    PRL(n);

    int m = 0;
    for (int j = VBINS-1; j >= 0; j--) {
	PRI(off);
	for (int i = 0; i < HBINS; i++) cerr << disp[i][j];
	if (j < VBINS-1 && m <= n) {
	    dyn *l = list[m];
	    if (l && l->is_valid())
		fprintf(stderr, "   %3d: %7s %10.2e",
			m, l->format_label(), l->get_mass());
	    else
		cerr << "   invalid l = " << l << " for m = " << m << endl;
	    m++;
	}
	cerr << endl << flush;
    }
    cerr << endl << flush;

#if 1

    // Add an array of all interparticle separations if not too large.

    if (n <= 5) {

	// Make a list of all nodes and leaves involved.

	int nnodes = 1000;

	while (nnodes > 10) {
	  nnodes = 0;
	  for (m = 0; m <= n; m++) {
	    if (list[m]) {
	      for_all_nodes(dyn, list[m], bb) nnodes++;
	    }
	  }
	  if (nnodes > 10) n--;
	  if (n <= 0) break;
	}

	if (n > 0) {

	  PRL(nnodes);

	  dyn ** nodes = new dynptr[nnodes];

	  int i = 0;
	  for (m = 0; m <= n; m++) {
	    if (list[m]) {
	      for_all_nodes(dyn, list[m], bb) nodes[i++] = bb;
	    }
	  }

	  PRI(14);
	  for (i = 0; i < nnodes; i++) {

	    // Label placement is Steve's aesthetic judgement (6/00)...

	    int len = strlen(nodes[i]->format_label());
	    if (len > 9) len = 9;
	    int ind = 7-len;
	    if (len > 3) ind = 6 - (1+len)/2;

	    PRI(ind);
	    fprintf(stderr, "%*s", len, nodes[i]->format_label());
	    PRI(10-ind-len);
	  }
	  cerr << endl;

	  for (int j = 0; j < nnodes; j++) {
	    fprintf(stderr, "%13s ", nodes[j]->format_label());
	    for (i = 0; i < nnodes; i++)
		if (i == j)
		    fprintf(stderr, "          ");
		else
		    fprintf(stderr, "%10.2e",
		      abs(something_relative_to_root(nodes[i], &dyn::get_pos) -
			  something_relative_to_root(nodes[j], &dyn::get_pos)));
	    cerr << endl;
	  }
	  cerr << endl;
	}
    }
#endif
}

#else

main(int argc, char ** argv)
{
    int  k = 3;
    int n = 5;
    int seed = 0;
    dyn *b;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "a:n:s:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {
	    case 'a': k = atoi(poptarg);
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
    }            

    srandinter(seed);

    while (b = get_dyn()) {
	dyn *bi = b->get_oldest_daughter();
	int nd = (int)(randinter(0, 1) * b->n_daughters());
	for (int i = 0; i < nd-1; i++) bi = bi->get_younger_sister();
        plot_stars(bi, n, k);
	rmtree(b);
    }
}

#endif
