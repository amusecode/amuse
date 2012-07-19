
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Set masses of specified particle(s) in input snapshot to specified values.
////
//// Usage:  setmass -l l1 -m mass1 -l l2 -m mass2 ...
////
//// Options:
////         -l    specify label of next particle to modify [no default]
////         -m    specify new mass for particle [no default]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// Steve McMillan, July 1998

#include "node.h"

#ifdef TOOLBOX

// leaf_with_label: Return a pointer to the leaf with the specified name.

local node* leaf_with_label(node* b, char* label)
{
    for_all_leaves(node, b, bi) {
	if (bi->get_name() != NULL) {
	    if (strcmp(bi->get_name(), label) == 0)
		return bi;
	} else if (bi->get_index() >= 0) {
	    char id[64];
	    sprintf(id, "%d", bi->get_index());
	    if (strcmp(id, label) == 0)
		return bi;
	}
    }
    return NULL;
}

local void set_mass_by_label(node* b, char* label, real mass)
{
    // Locate the particle to be split.

    node* bi = leaf_with_label(b, label);

    if (bi == NULL)
	cerr << "Warning: particle \"" << label << "\" not found.\n";
    else
	bi->set_mass(mass);

}

int main(int argc, char ** argv)
{
    char label[64];
    label[0] = '\0';

    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.6 $", _SRC_);

    node *b;
    b = get_node();

    b->log_history(argc, argv);

    // Parse the command line by hand, modifying the system as we go.

    int i = 0;
    while (++i < argc)
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {

		case 'l': strcpy(label, argv[++i]);
			  break;

		case 'm': set_mass_by_label(b, label, (real)atof(argv[++i]));
			  break;

		default:  get_help();
	    }

    // Erase any info on initial mass, virial radius or total energy,
    // as we do not recompute these quantities here.

    if (find_qmatch(b->get_log_story(), "initial_mass"))
	rmq(b->get_log_story(), "initial_mass");
    if (find_qmatch(b->get_log_story(), "initial_total_energy"))
	rmq(b->get_log_story(), "initial_total_energy");
    if (find_qmatch(b->get_log_story(), "initial_rvirial"))
	rmq(b->get_log_story(), "initial_rvirial");

    put_node(b);
    return 0;
}

#endif
