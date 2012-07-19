
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Determine the max density center position and velocity for the 
//// input N-body system.  The max density center is defined as the 
//// position of the particle with the highest local density.
//// Densities are not computed here -- run compute_density before 
//// invoking compute_max_cod.  Density center position and velocity
//// are written to the dyn story of the top-level node; they are also
//// optionally returned as function arguments in the library version.
//// Note: The computed density center is defined in absolute terms,
//// .and so includes the pos and vel of the parent node.
////
//// Usage: compute_max_cod [OPTIONS] < input > output
////
//// Options:   
////		   -c    add a comment to the output snapshot [false]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

//-----------------------------------------------------------------------------
//   version 1:  May 1989   Piet Hut
//   version 2:  Nov 1994   Piet Hut
//   version 3:  Jul 1996   Steve McMillan & Jun Makino
//.............................................................................
//   non-local function: 
//      compute_max_cod
//.............................................................................
//   see also: density.C
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  compute_max_cod -- Returns the position and velocity of the density
//		       center of an N-body system.
//-----------------------------------------------------------------------------

#define MAX_MSG_COUNT 5
#define TTOL 1.e-12				// arbitrary tolerance

static bool print_msg = true;
static int msg_count = 0;

void compute_max_cod(dyn *b, vec& pos, vec& vel)
{
    // First see if the data are already known.

    if (find_qmatch(b->get_dyn_story(), "density_center_type"))
	if (streq(getsq(b->get_dyn_story(), "density_center_type"), "max")) {
	    if (twiddles(getrq(b->get_dyn_story(), "density_center_time"),
			 b->get_system_time(), TTOL)) {
		pos = getvq(b->get_dyn_story(), "density_center_pos");
		vel = getvq(b->get_dyn_story(), "density_center_vel");
		return;
	    }
	}

    real max_density = 0;

    pos = 0;
    vel = 0;

//  for_all_leaves(dyn, b, d)
    for_all_daughters(dyn, b, d) {

	real dens_time = getrq(d->get_dyn_story(), "density_time");

	if (!twiddles(dens_time, b->get_system_time(), TTOL) && print_msg) {
	    warning("compute_max_cod: using out-of-date densities.");
	    PRC(d->format_label()), PRL(dens_time);
	    if (++msg_count > MAX_MSG_COUNT) print_msg = false;
	}

	real this_density = getrq(d->get_dyn_story(), "density");

	if (this_density > 0) {
	    if (max_density < this_density) {
		max_density = this_density;
		pos = d->get_pos();
		vel = d->get_vel();
	    }
	} else if (this_density <= -VERY_LARGE_NUMBER) {
	    if (print_msg) {
		warning("compute_max_cod: density not set.");
		PRL(d->format_label());
	    }
	    if (++msg_count > MAX_MSG_COUNT) print_msg = false;
	}
    }

    // Include the parent quantities.

    pos += b->get_pos();
    vel += b->get_vel();

    putsq(b->get_dyn_story(), "density_center_type", "max");
    putrq(b->get_dyn_story(), "density_center_time", b->get_system_time(),
	  HIGH_PRECISION);
    putvq(b->get_dyn_story(), "density_center_pos", pos);
    putvq(b->get_dyn_story(), "density_center_vel", vel);
}

void compute_max_cod(dyn *b)
{
    vec pos, vel;
    compute_max_cod(b, pos, vel);
}

#else

/*-----------------------------------------------------------------------------
 *  main  --  driver to use compute_max_cod() as a tool
 *-----------------------------------------------------------------------------
 */
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
		    "$Revision: 1.12 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn()) == NULL)
       err_exit("compute_max_cod: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        compute_max_cod(b);

        put_dyn(b);
	rmtree(b);
	b = get_dyn();
    }
}

#endif

// endof: compute_max_cod.C
