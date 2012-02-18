
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Cdetermine the mean mass of a power-law mass distribution.
////
//// Usage: m_bar [OPTIONS]
////
//// Options:
////            -f        multiply result by this factor [1]
////            -e/x      exponent [-2.35 (Salpeter)]
////            -l/L      lower mass limit [1]
////            -u/U      upper mass limit [1]
////             -v       set verbose mode [false]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "stdinc.h"

#ifdef TOOLBOX

main(int argc, char *argv[])
{
    real ml = -1.0, mu = -1.0, x = -2.35, factor = 1;
    bool verbose = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "e:f:l:L:u:U:vx:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c) {
	      case 'e':
	      case 'x':	x = atof(poptarg);
			break;
	      case 'f':	factor = atof(poptarg);
			break;
	      case 'l':
	      case 'L':	ml = atof(poptarg);
			break;
	      case 'u':
	      case 'U':	mu = atof(poptarg);
			break;
	      case 'v':	verbose = true;
			break;
	      case '?': params_to_usage(cerr, argv[0], param_string);
	                get_help();
	    	        exit(1);
	  }

    if (ml <= 0 || mu <= 0) exit(2);

    if (mu < ml) {
	real tmp = ml;
	ml = mu;
	mu = tmp;
    }

    real mbar = 1;

    if (ml < mu) {
	if (x != -1) {
	    mbar = ml * (1+x) * (pow(mu/ml, 2+x)-1)
				  / (pow(mu/ml, 1+x)-1) / (2+x);
	} else {
	    mbar = (mu - ml) / log(mu/ml);
	}
    }

    if (verbose) cout << "ml = " << ml << ", mu = " << mu << ", <m> = ";

    cout << mbar*factor << endl;
}

#endif
