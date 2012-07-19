
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Renumber named or unnumbered nodes to provide unique indices for all.
//// For use with applications requiring indices from input data having
//// only names (e.g. collisions).  For now, only handle leaves -- defer
//// the issue of determining unique indices for arbitrary CM nodes...
////
//// Usage: add_unique_index [OPTIONS]
////
//// Options:
////            -o      name-to-index renumbering option [not used]
////            -v      verbose output [quite operation]
////
//// Written by Steve McMillan
////
//// Report bugs to starlab@sns.ias.edu.

#include "node.h"

#ifdef TOOLBOX

local void reset_index(node *b, int fac, bool verbose)
{
    char *name = b->get_name();
    if (!name) return;

    // Make a list of separate integers found in the string.

    int ni = 0;
    int list[1024];
    char tmp[64];
    int nt = 0;
    bool plus = false;
    bool brace = false;

    for (int i = 0; i <= strlen(name); i++) {	// include trailing '\0'
	if (name[i] >= '0' && name[i] <= '9')
	    tmp[nt++] = name[i];
	else {
	    if (name[i] == '+')
		plus = true;
	    else if (name[i] == '<' || name[i] == '>')
		brace = true;
	    // PRC(name[i]); PRC(plus); PRC(brace); PRL(nt);
	    if (nt > 0) {
		tmp[nt] = '\0';
		list[ni++] = atoi(tmp);
		nt = 0;
	    }
	}
    }

    if (ni > 0) {

	// Index is fac times the sum of the last ni-1 elements, plus
	// the first.  Note that it makes more sense for the single
	// merger label a+b to be reinterpreted in the more general
	// form a<+1> for our purposes here.  Hence, if the string
	// contains "+" but no "<" or ">", then treat the last element
	// as 1.  Really should only do any of this if the label is of
	// the desired form, but we currently don't check carefully.

	// PRC(b->format_label()); PRC(plus); PRC(brace); PRL(ni);

	int sum = list[0];
	for (int i = 1; i < ni; i++) {
	    if (ni == 2 && plus && !brace)
		sum += fac;
	    else
		sum += fac*list[i];
	}
	b->set_index(sum);
	if (verbose)
	    cerr << "add_unique_index: index = " << b->get_index()
		 << " for leaf " << b->get_name() << endl;
    }
}

void add_unique_index(node * b, int option, bool verbose)
{
    // Add indices to unidentified nodes.  Start by finding a
    // convenient starting point for the extra numbers.

    int max_index = 0;
    for_all_nodes(node, b, bb)
	if (bb != b)
	    if (bb->get_index() > 0)
		max_index = Starlab::max(max_index, bb->get_index());

    if (verbose) {
	cerr << "add_unique_index: "; PRC(max_index);
    }

    // Start numbering at some well-separated value.

    int imax = (int) pow(10, (double)((int)log10((double)max_index) + 1));
    int index = imax+1;

    if (verbose) PRL(imax);

    // Add indices to otherwise unidentified leaves and create indices
    // for leaves with names.  Expectation is that the only current instance
    // of named leaves is collision products, which will have names of
    // the form a<+n>.

    for_all_leaves(node, b, bb)
	if (bb != b) {

	    if (bb->get_index() < 0 && bb->get_name() == NULL) {

		if (verbose)
		    cerr << "add_unique_index: attached index " << index
			 << " to unnamed leaf" << endl;
		bb->set_index(index++);

	    } else if (bb->get_name() != NULL)
		reset_index(bb, imax, verbose);

	}
}

int main(int argc, char ** argv)
{
    int option = 0;
    bool verbose = false;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "o:v";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.2 $", _SRC_)) != -1)
	switch(c) {

	    case 'o': option = atoi(poptarg);
		      break;
	    case 'v': verbose = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}


    node *b;
    while (b = get_node()) {
	b->log_history(argc, argv);
	add_unique_index(b, option, verbose);
	put_node(b);
	rmtree(b);
    }
}
#endif

