
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Sequentially merge all snapshots in the input stream.  Do *not*
//// rescale masses, positions, or velocities, but always set the
//// system center of mass to 0.  Root node information is retained
//// for all snapshots read.
////
//// Usage: merge_snaps [OPTIONS] < input > output
////
//// Options:
////            -n      don't renumber the particles [do renumber]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    bool renumber = true;
    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    const char *param_string = "n";

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1) {
	switch (c) {
	    case 'n':	renumber = false;
			break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    dyn *b, *root = new dyn;
    if (!root) err_exit("merge_snaps: can't create root node.");
    int count = 0;
    root->log_history(argc, argv);

    while (b = get_dyn()) {

	b->offset_com();		// add root pos/vel to daughters

	// Attach top-level nodes of b to root.

	dyn *last = root->get_oldest_daughter();
	while (last && last->get_younger_sister())
	    last = last->get_younger_sister();

	for_all_daughters(dyn, b, bb) {

	    // First node is a special case.

	    if (!last)
		root->set_oldest_daughter(bb);
	    else
		last->set_younger_sister(bb);

	    bb->set_elder_sister(last);
	    bb->set_parent(root);
	    last = bb;
	}

	// Merge the root log stories, with a prefix indicating which
	// initial snapshot this was.  This code is not completely general,
	// but should be sufficient for simple (typical) stories.

	story *sr = root->get_log_story();
	story *s = b->get_log_story();
	for (story * d = s->get_first_daughter_node(); d != NULL;
	     d = d->get_next_story_node())
	    if (!d->get_chapter_flag()) {
		char tmp[1024];
		sprintf(tmp, "  +%4.4d:  ", count);
		strncat(tmp, d->get_text(), 1013);
		tmp[1023] = '\0';
		add_story_line(sr, tmp);
	    }

	count++;
	delete b;
    }

    // Recompute the total mass and force the center of mass to 0.
    // Better renumber too, unless specifically suppressed..

    real mass = 0;
    int index = 0;
    for_all_daughters(dyn, root, bb) {
	mass += bb->get_mass();
	if (renumber) bb->set_index(++index);
    }
    root->set_mass(mass);

    // Set com using set_com() to update dyn story too.

    root->set_com();

    if (root->get_system_time() == 0.0) {

	// Rewrite initial_mass if it exists.  Delete initial_total_energy
	// and initial_rvirial.

	putrq(root->get_log_story(), "initial_mass", mass);
	rmq(root->get_log_story(), "initial_total_energy");
	rmq(root->get_log_story(), "initial_rvirial");
    }

    put_dyn(root);
    rmtree(root);
}

#endif
