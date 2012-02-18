
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Create a linked list of equal-mass nodes.
////
//// Usage: makenode [OPTIONS]
////
//// Options:
////         -m       specify total mass [1]
////         -n       specify number of nodes [1]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//	      Steve McMillan, July 1996

#include "node.h"

#ifndef TOOLBOX

node * mknode_mass(int n, real m)
{
    node *root = mknode(n);	// note that sequential labels are always set.
    root->set_root(root);

    root->set_mass(m);
    for_all_daughters(node, root, b) b->set_mass(m/n);

    return root;
}

#else

int main(int argc, char ** argv)
{
    int  n = 1;
    real m = 1.0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "m:n:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {

	    case 'm': m = atof(poptarg);
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    if (m <= 0) err_exit("mknodes: M > 0 required!");
    if (n <= 0) err_exit("mknodes: N > 0 required!");

    node * root = mknode_mass(n, m);

    root->log_history(argc, argv);

    put_node(root);
    rmtree(root);
    return 0;
}

#endif

/* end of: mknode.c */
