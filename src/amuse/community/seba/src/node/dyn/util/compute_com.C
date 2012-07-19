
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Determine the center of mass position and velocity of the 
//// input N-body system.  Center of mass position and velocity are
//// written to the dyn story of the top-level node; they are also
//// optionally returned as function arguments in the library version.
//// Note: The computed center of mass is defined in absolute terms,
//// and so includes the pos and vel of the parent node.
////
//// Usage: compute_com [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifndef TOOLBOX

#define DEBUG false

#define TTOL 1.e-12				// arbitrary tolerance

void compute_com(dyn *b, vec& com_pos, vec& com_vel)
{
    // First see if the data are already known.

    if (twiddles(getrq(b->get_dyn_story(), "com_time"),
		 b->get_system_time(), TTOL)) {
	com_pos = getvq(b->get_dyn_story(), "com_pos");
	com_vel = getvq(b->get_dyn_story(), "com_vel");
	// cerr << "    compute_com: using saved CM quantities" << endl;
	return;
    }

    real total_mass = 0;
    com_pos = 0;
    com_vel = 0;
    
    for_all_daughters(dyn, b, d) {
	total_mass += d->get_mass();
	com_pos += d->get_mass() * d->get_pos();
	com_vel += d->get_mass() * d->get_vel();
    }	

    com_pos /= total_mass;
    com_vel /= total_mass;

    // Include the parent quantities.

    com_pos += b->get_pos();
    com_vel += b->get_vel();

    // Use HIGH_PRECISION here because these quantities may be used
    // in detailed calculations elsewhere.

    if (DEBUG) {
	cerr << "compute_com: writing stories..." << endl << flush;
	PRL(0);
    }
    putrq(b->get_dyn_story(), "com_time", b->get_system_time(),
	  HIGH_PRECISION);
    if (DEBUG) {
	PRL(1);
	PRL(b);
	PRL(b->format_label());
	PRL(b->get_dyn_story());
	PRL(com_pos);
	b->print_dyn_story(cerr);
    }
    putvq(b->get_dyn_story(), "com_pos", com_pos,
	  HIGH_PRECISION);
    if (DEBUG) PRL(2);
    putvq(b->get_dyn_story(), "com_vel", com_vel,
	  HIGH_PRECISION);
}

void compute_com(dyn *b)
{
    vec tmp_pos, tmp_vel;
    compute_com(b, tmp_pos, tmp_vel);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use compute_com() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char *comment;
    dyn  *b;
    bool c_flag = FALSE;       // if TRUE: a comment given on command line

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.16 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn()) == NULL)
       err_exit("compute_com: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        compute_com(b);

	// Write system to stdout and get next system (if any).

        put_dyn(b);
	rmtree(b);
	b = get_dyn();
    }
}

#endif

// endof: compute_com.C
