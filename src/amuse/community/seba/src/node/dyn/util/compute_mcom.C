
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Determine the "modified center of mass" position and velocity 
//// of the input N-body system, recursively defined as the center 
//// of mass of the f [default 90%] particles closest to the center
//// of mass.  Center of mass position and velocity are written to the
//// dyn story of the top-level node; they are also optionally returned
//// as function arguments in the library version.
//// Note: The computed center of mass is defined in absolute terms, 
//// and so includes the pos and vel of the parent node.
////
//// Usage: compute_mcom [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -f    specify the fraction of particles to retain [0.9]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifndef TOOLBOX

typedef  struct
{
    real r2;
    dyn  *p;
} rp_pair, *rp_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rp_pair_ptr) pi)->r2 > ((rp_pair_ptr) pj)->r2)
        return +1;
    else if (((rp_pair_ptr)pi)->r2 < ((rp_pair_ptr)pj)->r2)
        return -1;
    else
        return 0;
}

#define TTOL 1.e-12				// arbitrary tolerance

void compute_mcom(dyn *b,
		  vec& pos, vec& vel,
		  real f,			// default = 0.9
		  int n_iter)			// default = 2
{
    // See if mcom_pos is already defined and up to date, with the
    // same value of f.

    if (find_qmatch(b->get_dyn_story(), "mcom_pos")
	&& getrq(b->get_dyn_story(), "mcom_f") == f
	&& twiddles(getrq(b->get_dyn_story(), "mcom_time"),
		    b->get_system_time(), TTOL)) {
	pos = getvq(b->get_dyn_story(), "mcom_pos");
	vel = getvq(b->get_dyn_story(), "mcom_vel");
	return;
    }

    // See if com_pos is already defined and up to date, or compute it.
    // Com quantities now include the root node (Steve, 6/03).

    if (find_qmatch(b->get_dyn_story(), "com_pos")
	&& twiddles(getrq(b->get_dyn_story(), "com_time"),
		    b->get_system_time(), TTOL)) {
	pos = getvq(b->get_dyn_story(), "com_pos");
	vel = getvq(b->get_dyn_story(), "com_vel");
    } else
	compute_com(b, pos, vel);

    pos -= b->get_pos();	// move to relative quantitles
    vel -= b->get_vel();

    // Use pos and vel as the starting point for computing the modified com.

    int n = 0;
    for_all_daughters(dyn, b, bb) n++;

    rp_pair_ptr rp = new rp_pair[n];

    if (rp == NULL) {
	cerr << "compute_mcom: "
	     << "not enough memory to store rp_pair array\n";
	return;
    }

    bool loop;
    int count = n_iter;

    if (f > 1) f = 1;
    if (f == 1) count = 0;

    do {
	loop = false;

	// Set up an array of radii.

	int i = 0;
	for_all_daughters(dyn, b, bi) {
	    rp[i].r2 = square(bi->get_pos() - pos);
	    rp[i].p = bi;
	    i++;
	}

	// Sort the array by radius.

	qsort((void *)rp, (size_t)i, sizeof(rp_pair), compare_radii);

	// Compute mpos and mvel by removing outliers.

	// Currently, apply mass independent weighting to all stars, but
	// let the weighting function go smoothly to zero at the cutoff.

	real r_max2i = 1/rp[(int)(f*n)].r2;
	real weighted_mass = 0;
	vec new_pos = 0, new_vel = 0;

	for (i = 0; i < f*n; i++) {
	    dyn *bi = rp[i].p;
	    real weight = bi->get_mass()
				* Starlab::max(0.0, 1 - rp[i].r2*r_max2i);
	    weighted_mass += weight;
	    new_pos += weight * bi->get_pos();
	    new_vel += weight * bi->get_vel();
	}

	if (weighted_mass > 0) {
	    pos = new_pos / weighted_mass;
	    vel = new_vel / weighted_mass;
	    if (count-- > 0) loop = true;
	}

    } while (loop);

    // Include the parent quantities.

    pos += b->get_pos();
    vel += b->get_vel();

    // Use INT/HIGH_PRECISION here because these quentities may be used
    // in detailed calculations elsewhere.

    putrq(b->get_dyn_story(), "mcom_time", b->get_system_time(),
	  HIGH_PRECISION);
    putvq(b->get_dyn_story(), "mcom_pos", pos,
	  INT_PRECISION);
    putvq(b->get_dyn_story(), "mcom_vel", vel,
	  INT_PRECISION);
    putrq(b->get_dyn_story(), "mcom_f", f);
    putiq(b->get_dyn_story(), "mcom_n_iter", n_iter);

    delete [] rp;
}

void compute_mcom(dyn *b,
		  real f,			// default = 0.9
		  int n_iter)			// default = 2
{
    vec pos, vel;
    compute_mcom(b, pos, vel, f, n_iter);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use compute_mcom() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    dyn * b;
    bool  c_flag = FALSE;       // if TRUE: a comment given on command line

    check_help();

    extern char *poptarg;
    int c;
    real f = 0.9;
    const char *param_string = "c:f:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.14 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'f': f = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn()) == NULL)
       err_exit("compute_mcom: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        compute_mcom(b, f);

	// Write system to stdout and get next system (if any).

        put_dyn(b);
	rmtree(b);
	b = get_dyn();
    }
}

#endif

// endof: compute_mcom.C
