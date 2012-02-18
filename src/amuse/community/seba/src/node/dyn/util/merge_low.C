
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Merge all low-level nodes.
////
//// Usage: merge_low [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -f    merge the specified fraction [1]
////              -o    specify how fraction f is chosen. [1]
////                    Option -1: random.
////                    Option  0: every 1/f-th binary.
////                    Option  n > 0: n-th block of size f.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  September 1999		Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void merge_low_level_nodes(dyn * b,	    // root node
			   real frac,	    // default = 1
			   int option)	    // default = 1
{
    if (frac <= 0 || frac > 1) return;

    int counter = 0;
    int nmerge = (int)(1/frac);		    // option = 0

    int start, stop;			    // option > 0
    if (option > 0) {
	int ncm = 0;
	for_all_daughters(dyn, b, bb)
	    if (bb->is_parent()) ncm++;
	nmerge = (int)(ncm*frac);
	start = (option-1)*nmerge;
	stop = start + nmerge;
    }
    int merged = 0;

    for_all_daughters(dyn, b, bb)
	if (bb->is_parent() && (
	    (option == -1 && randinter(0,1) <= frac) ||
	    (option ==  0 && ++counter == nmerge) ||
	    (option  >  0 && ++counter >= start && counter < stop)
	    )) {

	    dyn * od = bb->get_oldest_daughter();
	    dyn * yd = od->get_younger_sister();

	    // Assume that CM node is consistent, so all we need to
	    // do is delete the daughters.

	    rmtree(od);
	    rmtree(yd);

	    bb->set_oldest_daughter(NULL);
	    if (option == 0) counter = 0;
	    merged ++;
	}

    cerr << "merge_low: "; PRC(frac); PRC(option); PRL(merged);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  flatten_node() as a tool.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    real frac = 1;
    int option = 1;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:f:o:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.9 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'f': frac = atof(poptarg);
		      break;
	    case 'o': option = atoi(poptarg);
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

        merge_low_level_nodes(b, frac, option);
	put_dyn(b);
	rmtree(b);
    }
}

#endif
