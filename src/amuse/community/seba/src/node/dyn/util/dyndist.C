
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Get statistics over N 6-dimensional particle pair separations.  Compute
//// statistics for the separations between corresponding particles in two
//// N-body systems, using both configuration space and velocity data.  The
//// distribution of particle separations are reported as n-tiles.
////
//// Usage: dyndist [OPTIONS] < input > output
////
//// Options:   
////		  -n   specify number of n-tiles [4]
////
//// Written by the Starlab development group.
//// 
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  June 1990   Piet Hut           email: piet@iassns.bitnet
//                            Institute for Advanced Study, Princeton, NJ, USA
//                            Steve McMillan     email: steve@zonker.drexel.edu
//                            Physics Dept., Drexel Univ.,Philadelphia, PA, USA
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

//-----------------------------------------------------------------------------
//  print_sqrt_n_tiles  --  prints the values of the n_tiles .....
//-----------------------------------------------------------------------------

local void  print_sqrt_n_tiles(real * squares, int n_intervals, int n_squares)
    {
    int  i;
    real  factor = (n_squares - 1) / (real) n_intervals;

    for (i = 0; i <= n_intervals; i++)
        printf("  %.3g", sqrt(squares[(int) (factor*i + 0.5)]));
    printf("\n");
    }

//-----------------------------------------------------------------------------
//  realcompare  --  compares two number, used in qsort() above
//

local int  realcompare(const void * a, const void * b)
    {
    if (*((real *) a) > *((real *) b))
	return (1);
    else if (*((real *) a) < *((real *) b))
	return (-1);
    else
	return (0);
    }

//-----------------------------------------------------------------------------
//  dist_stats  --  ...
//                 note: the statistics are acquired on the first level in the 
//                       tree hierarchy, directly under the root node.
//                       If the positions and/or the velocities of the root 
//                       nodes are not equal, a warning message is printed.
//-----------------------------------------------------------------------------

void  dist_stats(dyn * b1, dyn * b2, int n_intervals)
    {
    real  rootrdist, rootvdist;
    vec  dr, dv;
    real * posdists;
    real * veldists;
    real * posdists_i;
    real * veldists_i;
    //    char *malloc();
    //    int  realcompare();


    // check the input data:

    int  n1, n2;
    dyn * bi1;
    dyn * bi2;
    for (n1 = 0, bi1 = b1->get_oldest_daughter(); bi1 != NULL;
         bi1 = bi1->get_younger_sister())
        n1++;
    for (n2 = 0, bi2 = b2->get_oldest_daughter(); bi2 != NULL;
         bi2 = bi2->get_younger_sister())
        n2++;
    if (n1 != n2)
        {
        cerr << "dyndiff: N1 = " << n1 << " != N2 = " << n2 << endl;
	exit(1);
	}

    // allocate memory:

    if ((posdists = (real *) malloc(n1 * sizeof(real))) == NULL)
	err_exit("dist_stats: not enough memory left for posdists[]");

    if ((veldists = (real *) malloc(n1 * sizeof(real))) == NULL)
	err_exit("dist_stats: not enough memory left for veldists[]");

    // compute the zero level distances:

    if ((rootrdist = abs(b1->get_pos() - b2->get_pos())) > 0)
        cerr << "dist_stats: the root nodes have a spatial separation of "
	     << rootrdist << endl;
    
    if ((rootvdist = abs(b1->get_vel() - b2->get_vel())) > 0)
        cerr << "dist_stats: the root nodes have a velocity separation of "
	     << rootvdist << endl;

    // compute the first level distances:

    for (bi1 = b1->get_oldest_daughter(), bi2 = b2->get_oldest_daughter(),
	 posdists_i = posdists, veldists_i = veldists;
	 bi1 != NULL;
         bi1 = bi1->get_younger_sister(), bi2 = bi2->get_younger_sister(),
	 posdists_i++, veldists_i++)
        {
	dr = bi1->get_pos() - bi2->get_pos();
	*posdists_i = dr*dr;
	dv = bi1->get_vel() - bi2->get_vel();
	*veldists_i = dv*dv;
	}

    qsort((void *)posdists, (size_t)n1, sizeof(real), realcompare);
    qsort((void *)veldists, (size_t)n1, sizeof(real), realcompare);

    print_sqrt_n_tiles(posdists, n_intervals, n1);
    print_sqrt_n_tiles(veldists, n_intervals, n1);
    }

main(int argc, char ** argv)
{
    int  n_intervals;
    bool  n_flag = FALSE;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "n:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'n': n_flag = TRUE;
		      n_intervals = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

    dyn * b1;
    dyn * b2;
    b1 = get_dyn();
    b2 = get_dyn();

    if (n_flag == FALSE)
	n_intervals = 4;                       // defaults: quartiles

    dist_stats(b1, b2, n_intervals);
}

#endif

// endof: dyndist.c
