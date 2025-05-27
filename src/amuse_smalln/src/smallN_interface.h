#include "hdyn.h"

// Functions in smallN_interface.cc:

int add_particle(hdyn *b, real mass, real radius,
		 vec pos, vec vel, int index_to_set);

hdyn *particle_with_index(hdyn *b, int index);

int remove_particle(hdyn *bb);
