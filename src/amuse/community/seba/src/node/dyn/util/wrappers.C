
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// wrapper.C:  Wrapper functions to provide a uniform interference to the
//	       starcluster package.  To be implemented (Steve, 6/03).

// "Bound" quantities need tidal data.
// "Core" quantities need density data.
//  Energy and virial radius need potential (and binary) data.

#include "dyn.h"

int  bound_number(dyn *b){return 0;}

real bound_mass(dyn *b){return 0;}

real total_energy(dyn *b){return 0;}

vec  total_angular_momentum(dyn *b, vec x, vec v)	// vec defaults = 0
{
    vec J = 0;
    for_all_daughters(dyn, b, bb)
	J += angular_momentum(bb, x, v);
    return J;
}

// Not currently implemented:

real core_radius(dyn *b){return 0;}

real core_mass(dyn *b){return 0;}

int  core_number(dyn *b){return 0;}

real virial_radius(dyn *b){return 0;}

real tidal_radius(dyn *b){return 0;}

real core_density(dyn *b){return 0;}
