// Additional functions to implement the interface for smallN.  Do
// this here in a separate file in order to keep interface.cc looking
// "standard" and to avoid changing any of the core smallN source
// files.

#include "smallN_interface.h"

int add_particle(hdyn *b, real mass, real radius,
		 vec pos, vec vel, int index_to_set)
{
    hdyn *bnew = new hdyn;
    bnew->set_mass(mass);
    bnew->set_radius(radius);
    bnew->set_pos(pos);
    bnew->set_vel(vel);
    if (!b->get_oldest_daughter()) {
	bnew->set_index(index_to_set < 0 ? 1 : index_to_set);
	b->set_oldest_daughter(bnew);
	bnew->set_parent(b);
    } else {
	int max_index = -1;
	bool set_index_OK = (index_to_set >= 0);
	for_all_daughters(hdyn, b, bb) {
	    if (bb->get_index() > max_index) max_index = bb->get_index();
	    if (bb->get_index() == index_to_set) set_index_OK = false;
	    if (bb->get_younger_sister() == NULL) {
		bnew->set_index(set_index_OK ? index_to_set : max_index + 1);
		bnew->set_parent(b);
		bnew->set_older_sister(bb);
		bb->set_younger_sister(bnew);
	    }
	}
    }
    return bnew->get_index();
}

// For small systems, probably OK just to search the list to find a
// particle with a given index.

hdyn *particle_with_index(hdyn *b, int index)
{
    for_all_daughters(hdyn, b, bb)
	if (bb->get_index() == index) return bb;
    return NULL;
}

int remove_particle(hdyn *bb)
{
    if (!bb) return -1;

    hdyn *b = bb->get_parent();
    if (bb == b->get_oldest_daughter())
	b->set_oldest_daughter(bb->get_younger_sister());
    if (bb->get_older_sister())
	bb->get_older_sister()->set_younger_sister(bb->get_younger_sister());
    if (bb->get_younger_sister())
	bb->get_younger_sister()->set_older_sister(bb->get_older_sister());

    rmtree(bb);
    return 0;
}
