
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Select a specific snapshot from the input stream.
////
//// Usage: getsnap [OPTIONS] < input > output
////
//// Options:    
////		  -c    add comment to snapshots [false]
////              -t    specify time of snapshot [0]
////
//// Written by Simon Portieges Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  July 2002   Simon Portieges Zwart

#include "dyn.h"

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    real  t = 0;   // default: first snapshot at time = 0
    bool  c_flag = FALSE;
    char  *comment;
    
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:t:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
    switch(c)
    {
	case 'c': c_flag = TRUE;
	          comment = poptarg;
	          break;
	case 't': t = atof(poptarg);
	          break;
	case '?': params_to_usage(cerr, argv[0], param_string);
		  get_help();
		  exit(1);
	      }            

    dyn *b;
    
    int node_read = 0;
    while (b = get_dyn())
    {
        node_read++;
        if (c_flag == TRUE)
	b->log_comment(comment);
        b->log_history(argc, argv);
	if(b->get_system_time()==t) {
	  put_dyn(b);
	  exit(1);
	}
	rmtree(b);
    }
}

/* endof: snapprune.c */
