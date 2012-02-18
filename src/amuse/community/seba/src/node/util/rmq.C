
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Remove a quantity from the dyn (or star or hydro) story of the
//// input snapshot(s) (FLAT TREE ONLY).
////
//// Usage rmq [OPTIONS]
////
//// Options:
////             -c     add a comment to snapshot [false]
////             -q     specify the quantity to remove [no default]
////             -h     remove from hydro story [false]
////             -s     remove from star story [false]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Jan 1993   Piet Hut

#include "node.h"

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  rm_all_dyn_q  --  only flat tree so far
 *-----------------------------------------------------------------------------
 */
void  rm_all_dyn_q(node * n, char * qname)
    {
    node * ni;

    for (ni=n->get_oldest_daughter(); ni != NULL; ni=ni->get_younger_sister())
	rmq(ni->get_dyn_story(), qname);
    }

/*-----------------------------------------------------------------------------
 *  rm_all_hydro_q  --  only flat tree so far
 *-----------------------------------------------------------------------------
 */
void  rm_all_hydro_q(node * n, char * qname)
    {
    node * ni;

    for (ni=n->get_oldest_daughter(); ni != NULL; ni=ni->get_younger_sister())
	rmq(ni->get_hydro_story(), qname);
    }

/*-----------------------------------------------------------------------------
 *  rm_all_star_q  --  only flat tree so far
 *-----------------------------------------------------------------------------
 */
void  rm_all_star_q(node * n, char * qname)
    {
    node * ni;

    for (ni=n->get_oldest_daughter(); ni != NULL; ni=ni->get_younger_sister())
	rmq(ni->get_star_story(), qname);
    }

/* endof: rmq.c */

/*-----------------------------------------------------------------------------
 *  main  --  removes a quantity with the name `name' from a dyn (or star or
 *            hydro) story, when invoked as  "rmq -q name".
 *     default: "rmq -q name" removes quantity from dyn story
 *   -h option: "rmq -h -q name" removes quantity from hydro story
 *   -s option: "rmq -s -q name" removes quantity from hydro story
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    bool  h_flag = FALSE;
    bool  s_flag = FALSE;
    bool  q_flag = FALSE;
    bool  c_flag = FALSE;
    char  *comment;
    char  *qname;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:q:hs";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'q': q_flag = TRUE;
		      qname = poptarg;
		      break;
	    case 'h': h_flag = TRUE;
		      break;
	    case 's': s_flag = TRUE;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            
    
    if (q_flag == FALSE)
	get_help();

    node *n;

    while (n = get_node()) {
        if (c_flag)
            n->log_comment(comment);
        n->log_history(argc, argv);

	if (h_flag)
	    rm_all_hydro_q(n, qname);
	else if (s_flag)
	    rm_all_star_q(n, qname);
	else
	    rm_all_dyn_q(n, qname);

	put_node(n);
	delete n;
    }
}

#endif
