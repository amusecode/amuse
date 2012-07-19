
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Set up a 3-body system corresponding to the Pythagorean problem.
////
//// Usage:  makepyth [OPTIONS]
////
//// Options:
////             -c    add a comment to the output snapshot [false]
////             -C    output data in 'col' format [no]
////             -i    number the particles sequentially [don't number]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Nov 1994   Piet Hut               email: piet@iassns.bitnet
//                          Institute for Advanced Study, Princeton, NJ, USA
//.............................................................................
//     litt: Szebehely, V. and Peters, C.F., 1967, Astron. J. 72, 876.
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    int  i;
    bool c_flag = FALSE;
    bool C_flag = FALSE;
    bool i_flag = FALSE;
    char *comment;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:Ci";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    dyn b, b1, b2, b3;
    if (i_flag) {
        b1.set_label(1);
        b2.set_label(2);
        b3.set_label(3);
    }
    b.set_oldest_daughter(&b1);
    b1.set_parent(&b);
    b2.set_parent(&b);
    b3.set_parent(&b);

    b1.set_younger_sister(&b2);
    b2.set_elder_sister(&b1);
    b2.set_younger_sister(&b3);
    b3.set_elder_sister(&b2);

    if (c_flag)
        b.log_comment(comment);
    b.log_history(argc, argv);

    b1.set_mass(3);
    b2.set_mass(4);
    b3.set_mass(5);

    b1.set_pos(vec(1,3,0));
    b2.set_pos(vec(-2,-1,0));
    b3.set_pos(vec(1,-1,0));

    if (C_flag) b.set_col_output(true);

    b.to_com();
    put_dyn(&b);
}

#endif

