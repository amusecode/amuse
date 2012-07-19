
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Renumber stars in a specific order
////
//// Usage: renumber [OPTIONS]
////
//// Options:
////            -c c      add a comment to the output snapshot [false]
////            -I/i      start numbering number               [1]
////            -M        renumber the stars in order of mass
////                      (highest mass=I/lowest mass=i)       [false]
////            -N        name the stars                       [false]
////            -S/s      single fixed number for all stars    [false]
////                      except if a star was already numbered
////
//// Written by Piet Hut and Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Jan 1993   Piet Hut
//	 	 Feb 2001   Simon Portegies Zwart

#include "node.h"
#include "util_io.h"

#ifndef TOOLBOX

typedef  struct {
    node* str;
    real  mass;
} nm_pair, *nm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_mass  --  compare the masses of two particles
//-----------------------------------------------------------------------------

local int compare_mass(const void * pi, const void * pj)
{
    if (((nm_pair_ptr) pi)->mass < ((nm_pair_ptr) pj)->mass)
        return(1);
    else if (((nm_pair_ptr)pi)->mass > ((nm_pair_ptr)pj)->mass)
        return(-1);
    else
        return(0);
}

void renumber(node* b, int istart, bool mass_order,
	      bool name_nodes, bool single_number) {

    int i;
    if(!mass_order) {

      i = istart;
      for_all_leaves(node, b, bj) {
	bj->set_label(i);
	if(!single_number)
	  i++;
      }
    }
    else {

      // Renumber the stars in order of mass.
      // Highest mass gets smallest number (strange choise, but).

      int n = b->n_leaves();
      nm_pair_ptr nm_table = new nm_pair[n];
      if (nm_table == NULL) {
	cerr << "renumber: "
	     << "not enough memory left for nm_table\n";
	return;
      }

      i=0;
      for_all_daughters(node, b, bi) {
	nm_table[i].str = bi;
	nm_table[i].mass = bi->get_mass();
	i++;
      }

      qsort((void *)nm_table, (size_t)n, sizeof(nm_pair), compare_mass);

      for (i=0; i<n; i++) {
	nm_table[i].str->set_index(istart+i);
      }
      delete []nm_table;

    }

    if(name_nodes) {				// the new option that wasted
						// Steve's time...!
	char tmp[MAX_INPUT_LINE_LENGTH];
	for_all_leaves(node, b, bj) {
	    PRL(bj->get_index());
	    if (bj->get_index() >= 0) {
		sprintf(tmp, "%d", bj->get_index());
		bj->set_name(tmp);
	    }
	}
    }
}

#else

int main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    bool M_flag = false;
    bool N_flag = false;
    bool S_flag = false;
    int istart = 1;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "MmNI:i:sSc:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.12 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'i':
	    case 'I': istart = atoi(poptarg);
		      break;
	    case 's': 
	    case 'S': S_flag = true;
		      break;
	    case 'm':
	    case 'M': M_flag = true;
		      break;
	    case 'N': N_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}


    node *b;
    b = get_node();
    b->log_history(argc, argv);

    renumber(b, istart, M_flag, N_flag, S_flag);

    put_node(b);
    rmtree(b);
}
#endif

