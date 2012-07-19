
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out the tree structure of the input snapshot(s) in normal form.
////
//// Options:
//// None.
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "node.h"

#ifndef TOOLBOX

local char* first_char(char* s)
{
    char* s1 = s;
    while (*s1 == '(' || *s1 == ')' || *s1 == ',') s1++;

    if (*s1 == '\0') err_exit("first_char: NULL label");

    return s1;
}

local bool precedes(char*s1, char* s2)
{
    // Return TRUE if s1 precedes s2 in our chosen ordering scheme.

    char* c1 = first_char(s1);
    char* c2 = first_char(s2);

    while (*c1 == *c2) {

	if (*c1 == '\0') err_exit("precedes: identical names!");

	c1++;
	c2++;
    }

    // Trick: short strings precede long ones because ")" and "," 
    // happen to precede alphanumeric characters in ASCII!

    if (*c1 > *c2) return FALSE;
    return TRUE;
}

local void sort_names(int n, char** list)	// Bubble sort!!
{
    for (int i = 0; i < n; i++) {
	for (int j = i + 1; j < n; j++) {
	    if (precedes(list[j], list[i])) {
		char* temp = list[i];
		list[i] = list[j];
		list[j] = temp;
	    }
	}
    }
}

void construct_node_name(node* b)
{
    if (!b->is_leaf()) {
	for_all_daughters(node, b, b1)
	    construct_node_name(b1);

	char** list = new char*[b->n_daughters()];

	int n = 0;
	for_all_daughters(node, b, b2) {
	    if (b2->get_name() == NULL) {

		// Create name if none exists.

		char temp2[128];
		sprintf(temp2, "#%d", b2->get_index());
		b2->set_label(temp2);
	    }
	    list[n++] = b2->get_name();
	}

	sort_names(n, list);

	// Don't worry about running out of space for now...

	char temp[1024];

	temp[0] = '(';
	temp[1] = '\0';
	int length = 2;

	for (int i = 0; i < n; i++) {

	    int sl = strlen(list[i]);
	    if (length + sl > 1024)
		err_exit("construct_node: buffer overflow.");

	    strcat(temp, list[i]);
	    if (i != n - 1) strcat(temp, ",");

	    length += sl + 1;
	}
	strcat(temp, ")");

	b->set_name(temp);
    }
}

void print_normal_form(node* b, ostream& s)
{
    construct_node_name(b);
    s << b->get_name() << endl;
}

char* get_normal_form(node* b)	// Unused...
{
    construct_node_name(b);
    return b->get_name();
}

/*===========================================================================*/

#else

/*-----------------------------------------------------------------------------
 *  main  --  driver to directly print out a tree structure
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
    {
    int i = 0;
    node *root;    // root node

    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.5 $", _SRC_);

    while (root = get_node())
	{
	// cerr << "snapshot #" << ++i << ":  n = "
	//      << root->n_leaves() << endl;
	print_normal_form(root, cerr);
	// rmtree(root);
	}
    }

#endif

/*===========================================================================*/

/* endof: print_normal.C */

