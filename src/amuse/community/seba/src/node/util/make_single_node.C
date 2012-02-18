
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Create a single node, to which other nodes can be added later,
//// using `add_daughter_node'.
////
//// Usage: make_single_node [OPTIONS]
////
//// Options:
////             -c  comment            optional comment
////             -m  mass               optional mass [1]
////             -i                     label nodes
////
//// Written by Piet Hut and Peter Teuben.
////
//// Report bugs to starlab@sns.ias.edu.


//   version 1:  Dec 1994   Piet Hut
//               1998-11-25 Peter Teuben        fixed help


#include "node.h"

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver to create directly a single node
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    bool  i_flag = FALSE;
    real m = 1;               // default mass: unity
    char  *comment;

    check_help();

    extern char *poptarg;
    int  c;
    const char *param_string = "c:im:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c) {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = TRUE;
		      break;
	    case 'm': m = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            
    
    node * root;

    root = new node();
    root->set_root(root);
    root->set_mass(m);

    if (c_flag == TRUE)
        root->log_comment(comment);
    root->log_history(argc, argv);

    if (i_flag)
        root->set_label(1);

    put_node(root);
    rmtree(root);
}

#endif

/*===========================================================================*/

/* endof: mk_single_node.C */

