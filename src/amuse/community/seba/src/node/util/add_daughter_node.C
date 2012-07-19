
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add one extra node to the tree structure in the input snapshot,
//// under a specified node.
////
//// Usage: mass_dist [OPTIONS]
////
//// Options:
////          -c     add a comment to snapshot [false]
////          -e     echo tree structure [false]
////          -i     specify index of node to add to [root]
////          -j     specify index for new node [none]
////          -m     specify node mass [1]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 1994   Piet Hut

#include "node.h"

//===========================================================================

#ifdef TOOLBOX

//---------------------------------------------------------------------------
//  main  --  driver to directly add one extra daughter node.
//-------------------------------------------------------------------------

main(int argc, char ** argv)
{
    int  i;
    int  j;
    bool  c_flag = FALSE;
    bool  e_flag = FALSE;     // echo flag: if true, echo tree structure
    bool  i_flag = FALSE;
    bool  j_flag = FALSE;
    real m = 1;               // default mass: unity
    char  *comment;

    check_help();

    extern char *poptarg;
    int  c;
    const char *param_string = "c:ei:j:m:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'e': e_flag = TRUE;
		      break;
	    case 'i': i_flag = TRUE;
		      i = atoi(poptarg);
		      break;
	    case 'j': j_flag = TRUE;
		      j = atoi(poptarg);
		      break;
	    case 'm': m = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}
    
    node * root;    // root node
    node * p;       // parent node
    node * d;       // older daughter node
    node * y;       // younger daughter node
    node * n;       // new daughter node

    root = get_node();

    if (root == NULL)
	err_exit("add_daughter_node: no input nodes provided");

    if (c_flag == TRUE)
        root->log_comment(comment);
    root->log_history(argc, argv);

    if (i_flag == FALSE)   // default parent: root
	p = root;
    else
	p = node_with_index(i, root);

    if (p == NULL)
	err_exit("add_daughter_node: no such parent");

    n = new node();
    n->set_mass(m);
    if (j_flag)
        n->set_label(j);

    d = p->get_oldest_daughter();
    if (d == NULL)
	p->set_oldest_daughter(n);
    else {
	y = d->get_younger_sister();
	while (y) {
	    d = y;
	    y = d->get_younger_sister();
	}
	d->set_younger_sister(n);
    }

    if (e_flag)
	root->pretty_print_tree(cerr);

    put_node(root);
    rmtree(root);
	}

#endif
