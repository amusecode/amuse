#ifndef STDINC_H
#define STDINC_H

// Standard inclusion of various system files. 

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <mpi.h>

using namespace std;
typedef double real;
typedef real (*real2)[3];	// 2-d array: "real2 x" = "real (*x)[3]"

#include "vec.h"		// 3-d vector type

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef NN
    #define NN 2		// conditional nearest-neighbor computation
#endif				// on GPU: not fully tested

#define _INFINITY_ 1.0e300
#define _TINY_	   (pow(2.0, -52))	  // double limit on powers of 2

#define VERY_LARGE_INTEGER (1<<30)	  // assume we are working with
					  // standard 4-byte integers, but
					  // also want -(this) to be legal
#define LOW_PRECISION	 3
#define STD_PRECISION	 6
#define INT_PRECISION	10
#define HIGH_PRECISION	15

#define PRI(n) for (int __k = 0; __k < n; __k++) cout << " "
#define PR(x)  cout << #x << " = " << x << " "
#define PRC(x) cout << #x << " = " << x << ",  "
#define PRL(x) cout << #x << " = " << x << endl << flush
#define PRLL(x) cout << #x << " = " << x << endl << endl << flush

//#define PRC(r,x) cout << r << ": " << #x << " = " << x << ",  " << flush
//#define PRL(r,x) cout << r << ": " << #x << " = " << x << endl << flush

//-----------------------------------------------------------------------------
//  pos_angle  --  recasts an angular variable into the range [0, TWO_PI)
//  sym_angle  --  recasts an angular variable into the range [-PI, PI)
//-----------------------------------------------------------------------------

#define  pos_angle(phi)    ((phi) - 2*M_PI * floor((phi)/(2*M_PI)))
#define  sym_angle(phi)    ((phi) - 2*M_PI * floor(((phi)+M_PI)/(2*M_PI)))
#define  sign(x)	   ((x > 0 ? 1 : (x < 0 ? -1 : 0)))

// Functions in util.cc:

real get_elapsed_time();
void get_cpu_time(real& user_time, real& system_time);
void warning(const char *s);
void err_exit(const char *s);
real randinter(real a = 0, real b = 1);

#endif

