
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Turn a text file into a snapshot.  Used in some demonstrations only.
//// File may be read top to bottom or left to right.
////
//// Usage:  makewrite [OPTIONS]
////
//// Options:
////           -c    add a comment to the output snapshot              [false]
////           -C    output data in 'col' format                          [no]
////           -d    toggle direction: top-bottom/left-right      [left-right]
////           -s    random seed                          [from system cloock]
////           -v    random velocity scale                                 [0]
////
//// Examples:
//// banner -w 40 STARLAB | makewrite -c "STARLAB"

#include "dyn.h"
#include <ctype.h>

#ifdef TOOLBOX

#define  SCALE_FACTOR   0.1

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = false;
    bool  C_flag = false;

    real v_rand = 0;

    int input_seed, actual_seed;
    bool s_flag = false;

    int dir = 0;

    check_help();
 
    extern char *poptarg;
    int c;
    const char *param_string = "c:Cdv:s:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.14 $", _SRC_)) != -1)
	switch(c) {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'd': dir = 1 - dir;
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'v': v_rand = atof(poptarg);
	    	      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    if (!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    dyn *b, *by, *bo;
    b = new dyn();
    b->set_root(b);

    bo = new dyn();
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    // Traditional banner writes from top to bottom.  Starlab banner
    // and file input writes left to right.  Want to translate into a
    // standard format.

    int i = 0;		// column number
    int j = 0;		// line number;

    int indx = 0;
    int n = 0;

    c = getchar();
    while(c != EOF) {
	while (c != '\n' && c != EOF) {
	    if (!isspace(c)) {
 		n++;
		real x, y;
		if (dir == 0) {

		    // Text reads left to right.

		    x = i * SCALE_FACTOR;
		    y = -j * SCALE_FACTOR;

		} else {

		    // Text reads top to bottom.

		    x = j * SCALE_FACTOR;
		    y = i * SCALE_FACTOR;

		}

		bo->set_pos(vec(x, y, 0));
		bo->set_vel(v_rand*vec(randinter(-1,1),
				       randinter(-1,1),
				       randinter(-1,1)));

		bo->set_label(++indx);

 		by = new dyn();
		bo->set_younger_sister(by);
		by->set_elder_sister(bo);
		bo = by;
	    }
	    i++;		// i increases with each character
	    c = getchar();
	}
	j++;			// j increases with each new line
	i = 0;
	if (c != EOF)
	    c = getchar();
    }

    bo = bo->get_elder_sister();
    delete by;
    bo->set_younger_sister(NULL);

    for_all_daughters(dyn, b, bb)
	bb->set_mass(1.0/n);

    if (c_flag == TRUE)
        b->log_comment(comment);
    b->log_history(argc, argv);

    if (C_flag) b->set_col_output(true);

    // Move to the center of mass frame (effectively centers the data,
    // so we don't have to worry too much about sizing and offsetting
    // the particles).

    b->to_com();
    put_dyn(b);
}

#endif

