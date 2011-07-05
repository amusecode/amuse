
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //          _\|/_
//=======================================================//            /|\ ~


//// Self-contained few-body integrator, using Starlab's
//// smallN_evolve, a fourth-order Hermite scheme incorporating a
//// modified unperturbed treatment of close approaches and time
//// symmetrization.

// This version integrates the input system in isolation, stopping
// when:
//
//	xxx
//
// Tidal errors due to the external system are currently not computed
// or absorbed, but should be (somewhere).

#include  "smallN2.h"

bool integrate(hdyn *b,
	       real t_end, real dt_struc,		// units = dt_dyn
	       real dt_dyn,				// dynamical time scale
							// default = 1
	       real eps2,				// default = 0
	       bool verbose,				// default = true
	       real break_r2,				// default = infinity
	       real dt_log,				// default = infinity
	       real dt_energy,				// default = infinity
	       real dt_snap)				// default = infinity
{
  // SmallN_evolve expects a flat tree.  Do that now.

  b->flatten_node();

  // Note: every time we return from smallN_evolve, we change the data
  // in the simulation, even though the function as written should be
  // exactly restartable.  The reason appars to be tha difference in
  // precision between data that remain in the processor (80 bits?)
  // and those that get restored from memory (64 bits).  As a
  // workaround, we always call kira_smallN for a specified length of
  // time dt_dyn (specified by the calling function); so long as
  // dt_struc is a multiple of this, we should have reproducible
  // results...

  t_end *= dt_dyn;
  dt_struc *= dt_dyn;

  bool over = false;
  while (b->get_time() < t_end) {
    real t_struc = b->get_time() + dt_struc;
    while (b->get_time() < t_struc) {
      // int status = 
      smallN_evolve(b, b->get_time() + dt_dyn, break_r2,
		    false, dt_log, dt_energy, dt_snap);
    }
    if (over = check_structure(b, eps2, verbose)) break;    // eps2 = 0 always?
  }
  cerr << "  integrate: "; PRL(over);
  return over;
}
