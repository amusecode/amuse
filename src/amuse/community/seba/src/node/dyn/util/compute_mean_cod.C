
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Determine the density center position and velocity for the input N-body
//// system.  The density center is defined analogously to the center of mass,
//// but instead of a mass weighting, the weighting used is proportional to
//// the square of the density (the original suggestion in the paper below(*),
//// to use a weighting factor linear in the density, may not converge well).
//// Densities are not computed here -- run compute_density before invoking
//// compute_mean_cod.  Density center position and velocity are written to
//// the dyn story of the top-level node; they are also optionally returned as
//// function arguments in the library version.
////
//// (*) Stefano Casertano and Piet Hut: Astroph.J. 298,80 (1985).
//// To use their recipe, k >= 2 is required.
////
//// Note: The computed density center is defined in absolute terms, and so
//// includes the pos and vel of the parent node.
////
//// Usage: compute_mean_cod [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////
//// Written by Piet Hut, Steve McMilland, and Jun Makino.
////
//// Report bugs to starlab@sns.ias.edu.
////
//-----------------------------------------------------------------------------
//   version 1:  Nov 1994   Piet Hut
//   version 2:  Jul 1996   Steve McMillan & Jun Makino
//.............................................................................
//.............................................................................
//   non-local function: 
//      compute_mean_cod
//.............................................................................
//   see also: density.C
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  compute_mean_cod -- Returns the position and velocity of the density
//		        center of an N-body system.
//-----------------------------------------------------------------------------

#define MAX_MSG_COUNT 5

static bool print_msg = true;
static int msg_count = 0;

#define TTOL 1.e-12				// arbitrary tolerance

void compute_mean_cod(dyn *b, vec& pos, vec& vel)
{
    // First see if the data are already known.

    if (find_qmatch(b->get_dyn_story(), "density_center_type"))
	if (streq(getsq(b->get_dyn_story(), "density_center_type"), "mean")) {
	    if (twiddles(getrq(b->get_dyn_story(), "density_center_time"),
			 b->get_system_time(), TTOL)) {
		pos = getvq(b->get_dyn_story(), "density_center_pos");
		vel = getvq(b->get_dyn_story(), "density_center_vel");
		return;
	    }
	}

    real total_weight = 0;

    pos = 0;
    vel = 0;
    
//  for_all_leaves(dyn, b, d)
    for_all_daughters(dyn, b, d) {

	real dens_time = getrq(d->get_dyn_story(), "density_time");

	if (print_msg
	    && !twiddles(dens_time, (real) b->get_system_time(), TTOL)) {
	    warning("compute_mean_cod: using out-of-date densities.");
	    PRL(d->format_label());
	    int p = cerr.precision(HIGH_PRECISION);
	    PRL(b->get_system_time());
	    PRL(dens_time);
	    cerr.precision(p);
	    if (++msg_count > MAX_MSG_COUNT) print_msg = false;
	}

	real this_density = getrq(d->get_dyn_story(), "density");

	if (this_density > 0) {
	    real dens2 = this_density * this_density;	    // weight factor
	    total_weight += dens2;
	    pos += dens2 * d->get_pos();
	    vel += dens2 * d->get_vel();
	} else if (this_density <= -VERY_LARGE_NUMBER) {
	    if (print_msg) {
		warning("compute_mean_cod: density not set.");
		PRL(d->format_label());
	    }
	    if (++msg_count > MAX_MSG_COUNT) print_msg = false;
	}
    }	

    if (total_weight > 0) {
	pos /= total_weight;
	vel /= total_weight;
    }

    // Include the parent quantities.

    pos += b->get_pos();
    vel += b->get_vel();

    putsq(b->get_dyn_story(), "density_center_type", "mean");
    putrq(b->get_dyn_story(), "density_center_time", b->get_system_time(),
	  HIGH_PRECISION);
    putvq(b->get_dyn_story(), "density_center_pos", pos);
    putvq(b->get_dyn_story(), "density_center_vel", vel);
}

void compute_mean_cod(dyn *b)
{
    vec pos, vel;
    compute_mean_cod(b, pos, vel);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use compute_mean_cod() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    dyn * b;
    bool  c_flag = FALSE;       // if TRUE: a comment given on command line

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.13 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn()) == NULL)
       err_exit("compute_mean_cod: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        compute_mean_cod(b);

        put_dyn(b);
	rmtree(b);
	b = get_dyn();
    }
}

#endif


// endof: compute_mean_cod.C
