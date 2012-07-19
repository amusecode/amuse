
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print times of all snapshots in the input stream.
////
//// Usage: list_snap [OPTIONS] < input > output
////
//// Options:
//// None.
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.5 $", _SRC_);

    dyn *b = NULL;
    int count = 0;

    while (b = get_dyn()) {

	cerr << "snap " << ++count;
	real time = getrq(b->get_dyn_story(),"t");

	if (time > -VERY_LARGE_NUMBER)

	    cerr << "  time = " << time << endl;

	else

	    cerr << "  time unknown" << endl;

	rmtree(b);
    }
}

#endif
