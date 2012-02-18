
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Bring all positions and velocities to center-of-mass frame.
//// Forces all daughter nodes to have a CM at rest at the origin,
//// and zeroes the root pos and vel.  Uses compute_com and writes
//// to the root dyn story.  Does not correct the virial radius in
//// the case of a tidal field.
////
//// Usage: to_com [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 1992   Piet Hut
//   version 2:  May 2003   Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void dyn::to_com()	// function is a special case of set_com
			// -- retain for compatibility
{
    set_com();		// default: pos = 0, vel = 0
}

#else

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        b->to_com();

	put_dyn(b);
	rmtree(b);
    }
}

#endif

// endof: to_com

