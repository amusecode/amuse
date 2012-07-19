
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute the rms distance in phase space between two (flat)
//// N-body systems.
////
//// Usage: dyndiff [OPTIONS] < input > output
////
//// Options:  
////		  -r    use spatial differences only.
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  Nov 1993   Piet Hut            email: piet@guinness.ias.edu
//                           Institute for Advanced Study, Princeton, NJ, USA
//.............................................................................
//  non-local function: 
//    dyndiff
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

//-----------------------------------------------------------------------------
//  dyndiff  --  determine the rms difference between the 6N coordinates of
//               two N-body systems
//               NOTE: simple version for flat tree only.
//-----------------------------------------------------------------------------

real  dyndiff(dyn * b1, dyn * b2, bool r_flag)
{
    //  quick fix to determine  n  for a flat tree:

    int  n1, n2;
    dyn * bi1;
    dyn * bi2;
    for (n1 = 0, bi1 = b1->get_oldest_daughter(); bi1 != NULL;
         bi1 = bi1->get_younger_sister())
        n1++;
    for (n2 = 0, bi2 = b2->get_oldest_daughter(); bi2 != NULL;
         bi2 = bi2->get_younger_sister())
        n2++;
    if (n1 != n2)
        {
        cerr << "dyndiff: N1 = " << n1 << " != N2 = " << n2 << endl;
	exit(1);
	}

    real sum_sq_diff = 0;     // sum of the squares of the 6N differences

    for (bi1 = b1->get_oldest_daughter(), bi2 = b2->get_oldest_daughter();
	 bi1 != NULL;
         bi1 = bi1->get_younger_sister(), bi2 = bi2->get_younger_sister())
        {
	vec  dr = bi1->get_pos() - bi2->get_pos();
	sum_sq_diff += dr*dr;
	
	if (!r_flag)
	    {
	    vec  dv = bi1->get_vel() - bi2->get_vel();
	    sum_sq_diff += dv*dv;
	    }
	}

    if (r_flag)
	return sqrt(sum_sq_diff / (3 * n1));
    else
	return sqrt(sum_sq_diff / (6 * n1));
}

main(int argc, char ** argv)
{
    bool  r_flag = FALSE;      // if TRUE: check only spatial differences

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "r";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c) {

	    case 'r': r_flag = TRUE;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            

    dyn * b1;
    dyn * b2;
    b1 = get_dyn();
    b2 = get_dyn();

    real  rmsdiff = dyndiff(b1, b2, r_flag);

    if (r_flag)
	cout << "  rms configuration space differences per coordinate:   "
	     << rmsdiff << endl;
    else
	cout << "  rms phasespace differences per coordinate:   "
	     << rmsdiff << endl;
}

#endif

// endof: dyndiff.C
