
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Plot input N-body system(s) on a non-graphics screen.
////
//// Usage: starplot [OPTIONS] < input > output
////
//// Options:   
////		  -a    view axis [z]
////              -c    clear the screen before each new plot [scroll screen]
////              -l    specify plot limits [get from first snapshot]
////              -n    number of lines to use [100]
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
//.............................................................................
//  non-local functions: 
//    starplot
//.............................................................................
//  Plot positions of all bodies in an N-body system, in a primitive way:
//  the screen is divided in character positions, and the representation is:
//  every  x  stands for one body; every  *  stands for more than one body.
//  This is very primitive, but independent of terminal or plotting software
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

#define  HBINS  56		// number of horizontal bins
#define  VBINS  28		// number of vertical bins

//-----------------------------------------------------------------------------
//  starplot.c  --  projects the positions of all particles on the screen
//                  accepts:  pn: a pointer to a nbody system,
//		               k: the number of the coordinate axis along which
//		                  the projection is directed, with {k,x,y}
//		                  having a right-handend orientation,
//		            lmax: maximum lengths alond the remaining axes.
//-----------------------------------------------------------------------------

void  starplot(dyn * b, int k, real lmax, int nlines)
    {
    int  n;
    int  kx, ky;	 // the coordinates projected onto the x- and y-axes
    register int  i, j;
    real lvmax;
    char  screen[VBINS][HBINS];		             // pixels on the screen
    dyn *bi;
    char  *str;
    char  plot_symbol = 'o';

    switch (k)
	{
	case 1: 
	    kx = 2; ky = 1; break;
	case 2:
	    kx = 0; ky = 2; break;
	case 3:
	    kx = 1; ky = 0; break;
	default: 
	    cerr << "starplot: k = " << k
		 << " illegal value; choose from {1, 2, 3}\n";
	    exit(0);
	}

    //  quick fix to determine  n  for a flat tree and get true maximum:

    real temp, true_max = 0;
    for (n = 0, bi = b->get_oldest_daughter(); bi != NULL;
         bi = bi->get_younger_sister()) {
        n++;

	temp = bi->get_pos()[kx];
	if (temp < 0) temp = -temp;
	if (temp > true_max) true_max = temp;
	temp = bi->get_pos()[ky];
	if (temp < 0) temp = -temp;
	if (temp > true_max) true_max = temp;
    }
    if (lmax <= 0) {
	temp = 1.0/128;
	while (temp < true_max) temp *= 2;
	lmax = temp;
    }	

    if (nlines > VBINS)
        nlines = VBINS;

    for (i = 0; i < nlines; i++)
	for (j = 0; j < HBINS; j++)
	    screen[i][j] = ' ';		// no particles within this pixel

    lvmax = lmax * ((real)nlines / (real)VBINS);

    for (bi=b->get_oldest_daughter(); bi != NULL; bi=bi->get_younger_sister())
	{
	str = &plot_symbol;		// default

	i = (int)(nlines * 
		 (lvmax + bi->get_pos()[kx]) / (2.0*lvmax));
	j = (int)(HBINS *
		 (lmax + bi->get_pos()[ky]) / (2.0*lmax));

	if ( i >=0 && i < nlines && j >=0 && j < HBINS)
	    {
	    if (screen[i][j] == ' ')
		screen[i][j] = *str;
	    else if (screen[i][j] == *str || screen[i][j] == (*str + 'A'-'a'))
		screen[i][j] = *str + 'A'-'a';
	    else
		screen[i][j] = '*';
	    }
	}

    putchar(' ');				// top line
    putchar(' ');
    putchar('+');
    for (j = 0; j < HBINS; j++)
	putchar('-');
    putchar('+');
    putchar('\n');

    putchar('^');				// second line: N
    putchar(' ');
    putchar('|');
    for (j = 0; j < HBINS; j++)
        putchar(screen[nlines - 1][j]);
    putchar('|');
    //    cout << "   N = " << n;
    printf("   N = %d", n);

    printf("\n");

    putchar('|');				// third line
    putchar(' ');
    putchar('|');
    for (j = 0; j < HBINS; j++)
        putchar(screen[nlines - 2][j]);
    putchar('|');
    putchar('\n');

    putchar('|');				// fourth line: range
    putchar(' ');
    putchar('|');
    for (j = 0; j < HBINS; j++)
        putchar(screen[nlines - 3][j]);
    putchar('|');
//    cout << "   max = " << lmax;
    printf("   max = %.2f", lmax);
    printf("\n");

    switch (k)
	{
	case 1: putchar('z'); break;
	case 2: putchar('x'); break;
	case 3: putchar('y'); break;
	}
    putchar(' ');				// fifth line: y-axis label
    putchar('|');
    for (j = 0; j < HBINS; j++)
        putchar(screen[nlines - 4][j]);
    putchar('|');
    putchar('\n');

    for (i = nlines - 5; i >= 0; i--)
	{
	putchar(' ');
	putchar(' ');
	putchar('|');
	for (j = 0; j < HBINS; j++)
	    putchar(screen[i][j]);
	putchar('|');
	putchar('\n');
	}

    putchar(' ');
    putchar(' ');
    putchar('+');
    for (j = 0; j < HBINS; j++)
	putchar('-');
    putchar('+');
    putchar('\n');

    for (j = 0; j < HBINS - 3; j++)
	putchar(' ');
    switch (k)
	{
	case 1: putchar('y'); break;
	case 2: putchar('z'); break;
	case 3: putchar('x'); break;
	}
    printf(" --->\n");
    }

#define  BIG_NUMBER  100

//-----------------------------------------------------------------------------
//  main  --  driver to use  starplot()  as a tool. 
//               The argument -a is interpreted as the axis along which to
//            view the N-body system: the number of the coordinate axis along
//            which the projection is directed, with {k,x,y} having a
//            right-handend orientation.
//               If an argument -l is provided, it defines the maximum
//	      lengths alond the remaining axes.
//               If an argument -n is provided, it defines the maximum
//	      number of lines on the screen.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    int  k;
    int  nlines;
    real lmax;
    dyn *b;
    bool  a_flag = FALSE;       // if TRUE: a projection axis was prescribed
    bool  c_flag = FALSE;       // if TRUE: clear screen before each plot
    bool  l_flag = FALSE;       // if TRUE: a maximum length was prescribed
    bool  n_flag = FALSE;       // if TRUE: number of lines was prescribed

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "a:cl:n:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'a': a_flag = TRUE;
		      k = atoi(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      break;
	    case 'l': l_flag = TRUE;
		      lmax = atof(poptarg);
		      break;
	    case 'n': n_flag = TRUE;
		      nlines = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            

    if (a_flag == FALSE)
	k = 3;                             // default

    if (l_flag == FALSE)
        lmax = 0;                          // default
    else
        if (lmax < 0)
	   cerr << "starplot: a maximum length of " << lmax
	        << " < 0 not allowed\n";

    if (c_flag) cout << "\33[H\33[J\33[H";

    if (n_flag == FALSE)
        nlines = BIG_NUMBER;               // default

    while (b = get_dyn()) {
	if (c_flag) cout << "\33[H";
        starplot(b, k, lmax, nlines);
	rmtree(b);
    }
}

#endif

// endof: starplot.C

