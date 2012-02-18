
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add a flag specifying that internal interactions should be suppressed
//// i.e. create a system of test particles.
////
//// Usage: set_ignore_internal [OPTIONS] < input > output
////
//// Options:
//// None.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 2001   Steve McMillan

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char *argv[])
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.6 $", _SRC_);

    dyn *b = get_dyn();
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);

    putrq(b->get_log_story(), "ignore_internal", 1);
    put_dyn(b);
}

#endif
