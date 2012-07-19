
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out the tree structure of the input snapshot(s).
////
//// Options:
//// None.
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 1994   Piet Hut

#include "node.h"

//===========================================================================

#ifdef TOOLBOX

//---------------------------------------------------------------------------
//  main  --  driver to directly print out a tree structure
//---------------------------------------------------------------------------

main(int argc, char ** argv)
{
    node *root;    // root node

    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.5 $", _SRC_);

    while (root = get_node()) {
	root->pretty_print_tree();
	rmtree(root);
    }
}

#endif
