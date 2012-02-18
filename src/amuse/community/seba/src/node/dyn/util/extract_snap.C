
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Find and print the last snapshot, or the first snapshot
//// following a specified time.  If -n is set to 0 and -v is set,
//// then all we do is count snapshots.
////
//// Usage: extract_snap [OPTIONS] < input > output
////
//// Options:    
////		  -c    add a comment to the output snapshot [false]
////              -t    specify time [last snapshot]
////              -n    specify number of snapshots to extract [1]
////              -v    set verbose mode [off]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.


#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = FALSE;	// if TRUE, a comment given on command line

    real  t_extract;
    bool  t_flag = FALSE;	// if TRUE, a time was specified

    bool  v_flag = FALSE;	// if TRUE, print snap times as read

    int   n = 1;
    bool  n_flag = false;	// if TRUE, a number was specified

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:n:t:v";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {
	    case 'c': c_flag = TRUE;
	    	      comment = poptarg;
		      break;
	    case 't': t_flag = TRUE;
		      t_extract = atof(poptarg);
		      break;
	    case 'n': n_flag = TRUE;
		      n = atoi(poptarg);
		      break;
	    case 'v': v_flag = TRUE;
		      break;
	    case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
		      exit(1);
	}            

    if (n < 0) err_exit("n < 0 specified.");
    if (n > 1 && !t_flag) {
	warning("n > 1 but no time specified -- 0 assumed.");
	t_extract = 0;
	t_flag = true;
    }

    node *b = NULL, *bp = NULL;
    int i = 0;

    // Need to ensure that the data are not changed by this function.
    // Currently, get_dyn() ends by calling check_and_correct_node(),
    // which may change pos and vel to force the system to the center
    // of mass frame.  Can fix (1) by suppressing this call in this
    // case, (2) setting a tolerance to avoid forcing if the system is
    // already "close" to the CM frame, or (3) use nodes here, which
    // will never interpret pos or vel and hence will not change them.
    // The only number we really need below is sytsem_time.  Options
    // (1) and (2) are awkward, because this function is a member
    // function and changing its arguments would propagate through the
    // other classes.  However, we could add a static local option to
    // the dyn instance of the function.  Option (3) seems cleanest.
    // To revert, change node to dyn and use b->get_system_time() to
    // determine time.
    //						(Steve, 7/04)

    while (b = get_node()) {

	real time;
	if (find_qmatch(b->get_dyn_story(), "real_system_time"))
	    time = getrq(b->get_dyn_story(), "real_system_time");
	else
	    time = getrq(b->get_dyn_story(), "system_time");

	// real time = b->get_system_time();

	i++;

	if (v_flag) cerr << "Snap time #" << i << " = " << time << endl;

	if (t_flag && time >= t_extract) {

	    if (n > 0) {

		if (c_flag == TRUE)
		    b->log_comment(comment);

		b->log_history(argc, argv);
		put_node(b);
	    }

	    if (--n <= 0) exit(0);
	}

	if (bp) rmtree(bp);
	bp = b;
    }

    if (n > 0 && !t_flag) {
	bp->log_history(argc, argv);
	put_node(bp);
    }

}

#endif
