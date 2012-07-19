
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Test/check Starlab node I/O functions.
////
//// Options:
//// None.
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "node.h"
#include "util_io.h"
#include "util_math.h"

#ifndef TOOLBOX

// Initialize all static node data here (?)

node * node::root           = NULL;

void node::print_static(ostream& s)		// default = cerr
{
    s << "root = " << root << endl;
}

bool node::check_and_correct_node(bool verbose)	// default = false
{
    // cerr << "node::check_and_correct_node: "; PRL(this);

    bool ok = true;

    if (oldest_daughter) {		// for node, just check that
					// the masses are consistent
	real m = 0;
	bool low = false;

	for_all_daughters(node, this, bb) {
	    if (bb->oldest_daughter)
		ok &= bb->check_and_correct_node(verbose);
	    m += bb->get_mass();
	}
	if (!ok) low = true;

	if (!twiddles(m, mass)) ok = false;
	mass = m;

	if (!ok && verbose && parent == NULL) {
	    cerr << "check_and_correct_node: applied ";
	    if (low)
		cerr << "low-level";
	    else
		cerr << "top-level";
	    cerr << " mass correction" << endl;
	}
    }

    return ok;
}

istream& node::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    // On entry, we have just read the "(Dynamics" [or "(D"] line
    // signifying the start of dyn_story input.  Keep reading and
    // storing until we encounter the corresponding closing
    // [END_DYNAMICS] line.

    while (get_line(s,input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (val) {

	    if (0) {   // trick to keep the else if() statements homogeneous

    	    }else if(!strcmp("m",keyword)){
		mass = strtod(val, NULL);
	    }else{
		add_story_line(dyn_story, input_line);
	    }
	}
    }
    return s;
}

ostream& node::print_dyn_story(ostream& s,
			       bool print_xreal,	// default = true
			       int short_output)	// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    put_real_number(s, "  m  =  ", mass);

    return s;
}

#else

main(int argc, char** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.13 $", _SRC_);

    node *b;
    while (b = get_node()) {
        put_node(b);
	pp2(b);
    }
    cerr << "Normal exit\n";
}
#endif
