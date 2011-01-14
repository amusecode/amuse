
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //          _\|/_
//=======================================================//            /|\ ~


//// SmallN2: Self-contained few-body integrator, using a fourth-order
//// Hermite scheme incorporating a modified unperturbed treatment of
//// close approaches and time symmetrization.  The program reads a
//// snapshot from standard input and optionally writes snapshot(s) to
//// standard output.  Optional periodic log output is sent to
//// standard error.  This version runs as a standalone program, or
//// can be incorporated into kira to handle multiple systems
//// containing perturbed hard binaries.
////
//// Usage: smallN2 [OPTIONS] < input > ouptut
////
//// Options:
////         -a    set accuracy parameter [0.03]
////         -d    set log output interval [0 --> no output]
////         -D    set snap output interval [0 --> no output]
////         -e    set external softening squared [0]
////         -E    set energy output interval [0 --> no output]
////         -g    set unperturbed limit [1.e-5]
////         -n    set number of symmetrization iterations [1]
////         -r    set termination radius [infinite]
////         -s    set structure check interval [at t_end]
////         -t    set termination time t_end [200]
////
//// Author: Steve McMillan.
////
//// Report bugs to steve@physics.drexel.edu

// This version of smallN integrates the input system in isolation,
// stopping when:
//
//	xxx
//
// Tidal errors due to the external system are currently not computed
// or absorbed, but should be (somewhere).

// NEED 1. better treatment of stability (Mardling)
//	2. "quarantining" of marginal cases that slip through
//	3. removal of escapers -- simplifies the quarantining logic??
//	4. time and size limits
//	5. reporting and storage of details

#include  "smallN2.h"

int main(int argc, char *argv[])
{
    extern char *poptarg;
    int c;
    char* param_string = "a:d:D:e:E:g:n:r:s:t:";

    real dt_log = 0;
    bool d_set = false;
    real dt_snap = 0;
    bool D_set = false;
    real dt_energy = 0;
    bool E_set = false;
    real break_r2 = VERY_LARGE_NUMBER;
    real t_end = 200;
    real dt_struc = 0;
    real eps2 = 0;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'a': set_smallN_eta(atof(poptarg));
		      break;
	    case 'd': dt_log = atof(poptarg);
	    	      d_set = true;
		      break;
	    case 'D': dt_snap = atof(poptarg);
	    	      D_set = true;
		      break;
	    case 'e': eps2 = atof(poptarg);
		      break;
	    case 'E': dt_energy = atof(poptarg);
	    	      E_set = true;
		      break;
	    case 'g': set_smallN_gamma(atof(poptarg));
		      break;
	    case 'n': set_smallN_niter(atoi(poptarg));
		      break;
	    case 'r': break_r2 = pow(atof(poptarg),2);
		      break;
 	    case 's': dt_struc = atof(poptarg);
		      break;
	    case 't': t_end = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      exit(1);
        }            

    PRL(get_smallN_eta());
    PRL(get_smallN_gamma());
    PRL(get_smallN_niter());
    PRL(break_r2);
    PRL(t_end);

    hdyn *b = get_hdyn();
    b->log_history(argc, argv);
    int seed = getiq(b->get_log_story(), "random number generator seed");
    cerr << "initial random seed = " << seed << endl;

    b->set_label("root");
    b->flatten_node();		// just in case...

    kira_options ko;
    ko.perturber_criterion = 2;
    b->set_kira_options(&ko);

    for_all_nodes(hdyn, b, bi) bi->set_time(0);

    // Use smallN_evolve() as the integrator.  Later, may want to add an
    // option to use the autoscaling features of integrate_multiple.


    if (dt_struc <= 0) dt_struc = t_end;

    // Arbitrarily set r_esc at 5 times the initial size of the system.
    // Currently not used...

    real m = 0;
    vec cmx = 0, cmv = 0;
    for_all_leaves(hdyn, b, bb) {
      m += bb->get_mass();
      cmx += bb->get_mass()*bb->get_pos();
      cmv += bb->get_mass()*bb->get_vel();
    }
    cmx /= m;
    cmv /= m;
    real rmax2 = 0;
    for_all_leaves(hdyn, b, bb) {
      real r2 = square(bb->get_pos() - cmx);
      if (r2 > rmax2) rmax2 = r2;
    }
//     real r_esc = 5*sqrt(rmax2);
//     dyn *bstore[2*b->n_leaves()];
//     int nstore = 0;

    int p = cerr.precision(INT_PRECISION);
    bool over = integrate(b, t_end, dt_struc, 1.0,
			  eps2, break_r2, dt_log, dt_energy, dt_snap);

    b->flatten_node();
    cerr << "initial random seed = " << seed << endl
	 << "final state (id, m, x, v):" << endl;
    cerr.precision(6);
    for_all_daughters(hdyn, b, bb)
	cerr << bb->format_label()
	     << " " << bb->get_mass()
	     << " " << bb->get_pos()
	     << " " << bb->get_vel() << endl;

    cerr.precision(p);

    if (!over) {
      real t_sys = b->get_system_time();
      b->set_system_time(b->get_time());
      my_sys_stats(b);
      b->set_system_time(t_sys);
      print_recalculated_energies((dyn*)b, 1);
    }
}
