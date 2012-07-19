
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Pass through only one out of every k snapshots from the input stream.
////
//// Usage snapprune [OPTIONS]
////
//// Options:
////            -c    add comment to snapshots [false]
////             -k    specify k [2]
////             -s    specify random seed [random from system clock]
////             -x    specify output snapshot after which to exit [at end]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Nov 1994   Piet Hut

#include "node.h"

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    int  k = 2;    // default: pass half of all snapshots
    int  s = 0;    // default: do not skip any snapshpt
    int  x = -1;   // default: exit after end of file
    bool  c_flag = FALSE;
    char  *comment;
    
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:k:s:x:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
    switch(c)
    {
	case 'c': c_flag = TRUE;
	          comment = poptarg;
	          break;
	case 'k': k = atoi(poptarg);
	          break;
	case 's': s = atoi(poptarg);
	          break;
	case 'x': x = atoi(poptarg);
	          break;
	case '?': params_to_usage(cerr, argv[0], param_string);
		  get_help();
		  exit(1);
	      }            

    node *b;
    for (int i = 0; i < s; i++){
	cerr << " skipping snapshot " << i << endl;
	if (!forget_node()) exit(1);
    }
    
    
    int node_read = 0;
    while (b = get_node())
    {
        node_read++;
        if (c_flag == TRUE)
	b->log_comment(comment);
        b->log_history(argc, argv);
	put_node(b);
	rmtree(b);
	
        for (int i = 1; i < k; i++)
	  if (!forget_node()) {
	    exit(1);
	  }

	if (x>0 && node_read>=x)
	  exit(1);
    }
}

/* endof: snapprune.c */
