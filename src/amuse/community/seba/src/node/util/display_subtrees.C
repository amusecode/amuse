
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Output all subtrees in input snapshot(s), neglecting the root
//// and top-level nodes.
////
//// Usage: mass_dist [OPTIONS]
////
//// Options:
////         -c    add a comment to snapshot [false]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Nov 1994   Piet Hut

#include "node.h"

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver, to directly display subtrees.
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    check_help();

    extern char *poptarg;
    int  c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            
    
    node *root;

    while (root = get_node()) {
	for_all_daughters(node, root, daughter)
	    if (daughter->is_parent()) {
		if (c_flag == TRUE)
		    daughter->log_comment(comment);
		daughter->log_history(argc, argv);

		put_node(daughter);
	    }
	rmtree(root);
    }
}

#endif

/*===========================================================================*/

/* endof: display_subtrees.C */

