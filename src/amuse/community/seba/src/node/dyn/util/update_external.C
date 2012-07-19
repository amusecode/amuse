
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  update_external.C: manually update parameters of external
//		       (currently only Plummer) potential.
//
//.............................................................................
//    version 1:  Aug 2007, Steve McMillan
//.............................................................................
//
// Global functions:
//
//	void update_external
//
//.........................................................................

#include "dyn.h"

// NOTES:  0. Currently only modifies a Plummer external field.
//	      Terminology: Plummer potential is -M / sqrt(r^2 + a^2).
//
//	   1. This is intended as a simple function to allow easy user
//	      access to the external field.
//
//	   2. There is currently no provision for automatic correction
//	      of the energy when the potential changes.

static bool initialized = false;
static real initial_M = 0;
static real initial_a_sq = 0;

void update_external(dyn *b)			// root node
{
    if (!b->get_plummer()) return;		// Plummer only

    if (!initialized) {
	initial_M = b->get_p_mass();		// initial M
	initial_a_sq = b->get_p_scale_sq();	// initial a^2
	initialized = true;
    }

    real t = b->get_system_time();		// system time

    real M = initial_M;
    real a_sq = initial_a_sq;

    //--------------------------------------------------------------------
    //
    // Modify the global Plummer parameters at time t.  Specify the
    // *relative* scaling of M and a^2.  (Default is to do nothing.)

    // To use, replace #if 0 by #if 1 and apply the desired scaling.

#if 0

    M = initial_M * (1 - t/50.);
    a_sq = initial_a_sq * (1 + t/10.);

#endif

    //--------------------------------------------------------------------

    if (M < 0) M = 0;
    if (a_sq < 0) a_sq = 0;

    if (M != initial_M || a_sq != initial_a_sq) {
	cerr << endl << "update_plummer: new "; PRC(M); PRL(a_sq);
	b->set_p_mass(M);
	b->set_p_scale_sq(a_sq);
    }
}
