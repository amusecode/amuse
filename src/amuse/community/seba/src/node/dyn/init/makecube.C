
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct a simple homogeneous cube, with (default)
////
////              M = 1, T/U = -1/2, E = -1/4.
////
//// If the "-u" flag is set, the particles are left unscaled, with
//// masses 1/n, positions uniformly distributed in [-L, L], and
//// velocities uniformly distributed in a range giving approximate
//// virial equilibrium.
////
//// Usage:  makecube [OPTIONS]
////
//// Options:
////              -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -i    number the particles sequentially [don't number]
////              -l    write cube size to dyn story [don't write]
////              -L    specify cube size (+/-L) [1]
////              -n    specify number of particles [no default]
////              -o    echo value of random seed [don't echo]
////              -s    specify random seed [random from system clock]
////              -u    leave unscaled [scale to E=-1/4, M = 1, R = 1]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

local void makecube(dyn *root, int n,
		    real L = 1,
		    int u_flag = false)
{
    if (L <= 0) return;

    real pmass = 1.0 / n;

    // Factor scaling the velocity places the system in approximate
    // virial equilibrium without scaling.

    // Not clear what the proper velocity distribution should be...

    real vfac = 0.7/sqrt(L);

    for_all_daughters(dyn, root, bi) {
	bi->set_mass(pmass);	
        bi->set_pos(L*vec(randinter(-1,1),randinter(-1,1),randinter(-1,1)));
        bi->set_vel(vfac*vec(randinter(-1,1),randinter(-1,1),randinter(-1,1)));
    }

    // Transform to center-of-mass coordinates, and optionally
    // scale to standard parameters.

    root->to_com();
    root->set_mass(1);
    putrq(root->get_log_story(), "initial_mass", 1.0);

    if (!u_flag && n > 1) {

        real potential, kinetic;

	// Note: scale_* operates on internal energies.

	get_top_level_energies(root, 0.0, potential, kinetic);
	scale_virial(root, -0.5, potential, kinetic);	// scales kinetic
	real energy = kinetic + potential;
	scale_energy(root, -0.25, energy);		// scales energy
	putrq(root->get_log_story(), "initial_total_energy", -0.25);
	putrq(root->get_log_story(), "initial_rvirial", 1.0);
	putrq(root->get_dyn_story(), "total_energy", -0.25);
    }
}

#define  SEED_STRING_LENGTH  60

main(int argc, char ** argv) {
    int  i;
    int  n;
    int  input_seed, actual_seed;
    bool c_flag = false;
    bool C_flag = false;
    int  n_flag = FALSE;
    int  s_flag = FALSE;
    int  i_flag = FALSE;
    int  l_flag = FALSE;
    int  o_flag = FALSE;
    int  u_flag = FALSE;   // if true, leave scaling untouched

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    real  cube_size = 2;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:CilL:n:os:u";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.15 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'l': l_flag = true;
		      break;
	    case 'L': cube_size = 2*atof(poptarg);
		      break;
	    case 'n': n_flag = true;
		      n = atoi(poptarg);
		      break;
	    case 'o': o_flag = true;
                      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'u': u_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            
    
    if (!n_flag) {
        cerr << "makecube: must specify the number # of";
	cerr << " particles with -n#\n";
	exit(1);
    }
    
    if (n < 1) {
        cerr << "makecube: n < 1 not allowed\n";
	exit(1);
    }

    dyn *b, *by, *bo;
    b = new dyn();
    b->set_root(b);

    bo = new dyn();
    if (i_flag) bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);
    if (l_flag)
	putrq(b->get_log_story(), "cube_size", cube_size);

    // Create the tree.

    for (i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag) by->set_label(i+1);
	by->set_parent(b);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    if (C_flag) b->set_col_output(true);
    if (c_flag) b->log_comment(comment);

    b->log_history(argc, argv);

    if (!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makecube: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    if (n > 0) makecube(b, n, cube_size/2, u_flag);

    put_dyn(b);
}

#endif

