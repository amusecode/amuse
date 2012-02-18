
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Plot Lagrangian radii for input N-body system(s), based on the
//// geometric center; takes output from "lagradplot -q" and smooths
//// over w output lines.  (Probably only marginally useful to use this
//// center -- better to use the density center or modified center of mass,
//// as in lagrad...)
//// Note: This really is a terrible kluge to smooth Lagrangian radii,
//// and in addition it is programmed in an opaque way, just to
//// get some quick and dirty results; it works is all I can say.
//// We really should clean this all up soon(ish) -- hope springs eternal.
////
//// Usage: lagradplotsmooth [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -w    number of output lines to smooth over
////
//// Examples:
////
//// mkplummer -n 128 -i -s 123 | ( kira -t 1000 -d 10 -D 0.2 -n 0 |
//// compute_density | to_cod | lagradplot -q > run128q.out ) >& run128q.log
////
//// lagradplotsmooth -w 30 < run128q.out
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.


//.............................................................................
//    version 1:  Dec 2000   Piet Hut
//.............................................................................
//  see also: lagradplot.C
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

#define  MAX_NUMBER_OF_INPUT_LINES     100
#define  MAX_NUMBER_OF_COLUMNS         180
#define  MAX_NUMBER_OF_OUTPUT_COLUMNS   79

//-----------------------------------------------------------------------------
//  use_square_window_smoothing  --  to smooth Lagrangian radii output
//-----------------------------------------------------------------------------

void  use_square_window_smoothing(
char in_screen[MAX_NUMBER_OF_INPUT_LINES][MAX_NUMBER_OF_COLUMNS], int w)
    {
    int c;
    int i, j;
    int k;
    int w_tot;
    w_tot = 0;
    real radius[12];
    for (i= 0; i < 12; i++)
        radius[i] = 0;

    bool input_flag;
    bool one_flag;
    bool two_flag;
    bool five_flag;
    for (i= 0; i < MAX_NUMBER_OF_INPUT_LINES; i++)
        {
        k = 2;
        input_flag = FALSE; 
        one_flag = FALSE; 
        two_flag = FALSE; 
        five_flag = FALSE; 
        for (j= 0; j < MAX_NUMBER_OF_COLUMNS; j++)
            {
            while (in_screen[i][j] == ' ' && j + 1 < MAX_NUMBER_OF_COLUMNS)
                j++;
            if (in_screen[i][j] == '`')
                {
                input_flag = TRUE;
                radius[0] += j;
                one_flag = TRUE;
                }
            if (in_screen[i][j] == '\"')
                {
                input_flag = TRUE;
                if (! one_flag)
                    {
                    radius[0] += j;
                    one_flag = TRUE;
                    }
                radius[1] += j;
                two_flag = TRUE;
                }
            if (in_screen[i][j] == ':')
                {
                input_flag = TRUE;
                if (! one_flag)
                    {
                    radius[0] += j;
                    one_flag = TRUE;
                    }
                if (! two_flag)
                    {
                    radius[1] += j;
                    two_flag = TRUE;
                    }
                radius[2] += j;
                five_flag = TRUE;
                }
            if (in_screen[i][j] == '|')
                {
                input_flag = TRUE;
                k++;
                if (! one_flag)
                    {
                    radius[0] += j;
                    one_flag = TRUE;
                    }
                if (! two_flag)
                    {
                    radius[1] += j;
                    two_flag = TRUE;
                    }
                if (! five_flag)
                    {
                    radius[2] += j;
                    five_flag = TRUE;
                    }
                radius[k] += j;
                }
            if (in_screen[i][j] == 'X')
                {
                input_flag = TRUE;
                k++;
                if (! one_flag)
                    {
                    radius[0] += j;
                    one_flag = TRUE;
                    }
                if (! two_flag)
                    {
                    radius[1] += j;
                    two_flag = TRUE;
                    }
                if (! five_flag)
                    {
                    radius[2] += j;
                    five_flag = TRUE;
                    }
                radius[k] += j;
                k++;
                radius[k] += j;
                }
            if (in_screen[i][j] == '@')
                {
                input_flag = TRUE;
                k++;
                if (! one_flag)
                    {
                    radius[0] += j;
                    one_flag = TRUE;
                    }
                if (! two_flag)
                    {
                    radius[1] += j;
                    two_flag = TRUE;
                    }
                if (! five_flag)
                    {
                    radius[2] += j;
                    five_flag = TRUE;
                    }
                radius[k] += j;
                k++;
                radius[k] += j;
                k++;
                radius[k] += j;  /* hoping that there are no more than 3     */
				 /* points involved; this is uncertain       */
                }
            }
        if (input_flag)
            w_tot++;
        }

//    cerr << "w_tot = " << w_tot << "\n";

    for (i= 0; i < 12; i++)
        radius[i] /= w_tot;

    j = 0;
    k = 3;
    while (j < MAX_NUMBER_OF_OUTPUT_COLUMNS)
        {
        if (radius[3] > j + 0.5)
            {
            if (radius[2] <= j + 0.5 && radius[2] > j - 0.5)
                printf(":");
            else if (radius[1] <= j + 0.5 && radius[1] > j - 0.5)
                printf("\"");
            else if (radius[0] <= j + 0.5 && radius[0] > j - 0.5)
                printf("`");
            else
		printf(" ");
            }
        else
            {
            if (radius[k] <= j + 0.5 && radius[k] > j - 0.5)
                {
                k++;
                if (radius[k] > j + 0.5)
		    printf("|");
                else if (radius[k + 1] > j + 0.5)
		    {
                    k++;
		    printf("X");
		    }
                else
		    {
                    k++;
                    k++;
		    printf("@");
		    }
                }
            else
		printf(" ");
            }
        j++;
        }
    printf("\n");
    }

//-----------------------------------------------------------------------------
//  main
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = FALSE;      /* if TRUE: a comment given on command line   */
    bool  w_flag = FALSE;      /* if TRUE: followed by the number of         */
                               /*          snapshots over which to           */
			       /*	   average in a square window        */
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:w:";

    int w;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'w': w_flag = TRUE;
		      w = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            

    if (w > MAX_NUMBER_OF_INPUT_LINES)
        {
        cerr << "lagradplotsmooth: w = " << w << " > " <<
                MAX_NUMBER_OF_INPUT_LINES << " (maximum value)\n";
	exit(1);
	}

    char  in_screen[MAX_NUMBER_OF_INPUT_LINES][MAX_NUMBER_OF_COLUMNS];
    int  i,j;
    for (i= 0; i < MAX_NUMBER_OF_INPUT_LINES; i++)
        for (j= 0; j < MAX_NUMBER_OF_COLUMNS; j++)
            in_screen[i][j] = ' ';        

    c = 1;

    int local_w;

    while (c != EOF)
        {
        local_w = w;
        i = 0;
        while (local_w && c != EOF)
            {
//cerr << "w = " << w << "\n";
            j = 0;
            c = getchar();
            while (c != '\n' && c != EOF)
                {
                in_screen[i][j++] = c;
                c = getchar();
                }
            i++;
            local_w--;
            }

        if (c != EOF)
            use_square_window_smoothing(in_screen, w);
        }
    }

#endif

// endof: lagradplot.C
