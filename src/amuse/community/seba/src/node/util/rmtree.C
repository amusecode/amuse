
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  rmtree.C: recursively deletes a node together with all of its offspring
//.............................................................................
//    version 1:  Nov 1994   Piet Hut
//.............................................................................

#include "node.h"

void rmtree(node* b,
	    bool delete_b)	// default = true
{
    node* d = b->get_oldest_daughter();
    while (d) {
	node* tmp = d->get_younger_sister();
	rmtree(d);
	d = tmp;
    }

    if (delete_b) delete b;	// optionally leave node itself untouched
}

/* endof: rmtree.c */
