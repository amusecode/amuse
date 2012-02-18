
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Flatten a dyn tree to a single-level linked list under the root node.
////
//// Usage: flatten [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -C    force col output [take from input format]
////              -v    print diagnostic info [no info]
////
//// Written by Piet Hut
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  March 1994   Piet Hut

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  unbundle_node  --  ...
//                 note:
//                      the positions and velocities of the daughters of the
//                      node-to-be-unbundled are increased by the values of the
//                      corresponding quantities of the node-to-be-unbundled;
//                      the other physical quantities of the daughters remain
//                      unchanged, while those of the node-to-be-unbundled are
//                      to be discarded.
//                 note:
//                      the parent, daughters, and granddaughters really should
//                      be synchronized in time; to be implemented.
//-----------------------------------------------------------------------------

local int unbundle_node(dyn *ud)
{
    dyn *parent = ud->get_parent();

    if (parent == NULL)
	err_exit("unbundle(): no parent for this node");

    dyn *od = ud->get_oldest_daughter();

    if (od == NULL)
	return 0;                         // nothing left to unbundle

    vec pos = ud->get_pos();
    vec vel = ud->get_vel();

    for_all_daughters(dyn, ud, dd) {
	dd->inc_pos(pos);
	dd->inc_vel(vel);
    }

    // pointer adjustments:

    dyn *elder_sister = ud->get_elder_sister();

    if (elder_sister == NULL)
	parent->set_oldest_daughter(od);
    else {
	elder_sister->set_younger_sister(od);
	od->set_elder_sister(elder_sister);
    }

    dyn *d = od;
    dyn *pd;

    while (d) {
	d->set_parent(parent);
	pd = d;
	d = d->get_younger_sister();
    }
    
    dyn *younger_sister = ud->get_younger_sister();

    if (younger_sister) {
	younger_sister->set_elder_sister(pd);
	pd->set_younger_sister(younger_sister);
    }

    delete ud;
    return 1;
}

//-----------------------------------------------------------------------------
//  flatten_node  --  
//-----------------------------------------------------------------------------

int dyn::flatten_node()
{
    int n_flat = 0;

    for_all_daughters(dyn, this, d) {
	if (d->is_grandparent())
	    n_flat += d->flatten_node();

	n_flat += unbundle_node(d);
    }
    return n_flat;
}

//=============================================================================

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  flatten_node() as a tool.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool  c_flag = false;
    bool  C_flag = false;
    bool  verbose = false;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:Cv";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'v': verbose = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag)
            b->log_comment(comment);
        b->log_history(argc, argv);

	if (C_flag) b->set_col_output(true);

        int n_flat = b->flatten_node();

	if (verbose)
	    cerr << "time = " << b->get_system_time() << ", "
		 << n_flat << " nodes flattened" << endl;

	put_dyn(b);
	rmtree(b);
    }
}

#endif

// endof: flatten.C

