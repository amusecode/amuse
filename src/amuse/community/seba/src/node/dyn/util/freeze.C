
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Reduce all top-level velocities, leaving positions unchanged.
////
//// Usage: freeze [OPTIONS] < input > output
////
//// Options:     
////              -c    add a comment to the output snapshot [false]
////              -f    specify freeze factor [0]
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 1992   Piet Hut
//               Oct 2001   Steve McMillan

#include "dyn.h"

#ifdef TOOLBOX

local void freeze(dyn * b, real fac)
{
    for_all_daughters(dyn, b, bi)
	bi->scale_vel(fac);
}

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    real fac = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:f:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.11 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'f': fac = atof(poptarg);
	    	      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    dyn *b;

    while (b = get_dyn()) {
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        freeze(b, fac);

	// Clean up total and kinetic energies.

	story *sk = find_qmatch(b->get_dyn_story(), "kinetic_energy");
	story *se = find_qmatch(b->get_dyn_story(), "total_energy");

	if (sk) {

	    if (se) {
		real k = getrq(b->get_dyn_story(), "kinetic_energy");
		real e = getrq(b->get_dyn_story(), "total_energy");
		putrq(b->get_dyn_story(), "total_energy", e-k);
	    }

	    putrq(b->get_dyn_story(), "kinetic_energy", 0);
	    
	} else if (se)

	    rm_daughter_story(b->get_dyn_story(), se);

	put_dyn(b);
	rmtree(b);
    }
}

#endif

// endof: freeze.C
