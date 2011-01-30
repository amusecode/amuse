#ifndef STDINC_H
#define STDINC_H

// Standard inclusion of various system files. 

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <unistd.h>
#include <mpi.h>

using namespace std;
typedef double real;
typedef real (*real2)[3];	// 2-d array: "real2 x" = "real (*x)[3]"

#include "vec.h"		// 3-d vectory type

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef NN
    #define NN 2		// conditional nearest-neighbor computation
#endif				// on GPU: not fully tested

const real tiny = 1.e-25;
const real huge = 1.e25;

#define PRI(n) for (int __k = 0; __k < n; __k++) cout << " "
#define PR(x)  cout << #x << " = " << x << " "
#define PRC(x) cout << #x << " = " << x << ",  "
#define PRL(x) cout << #x << " = " << x << endl << flush

//#define PRC(r,x) cout << r << ": " << #x << " = " << x << ",  " << flush
//#define PRL(r,x) cout << r << ": " << #x << " = " << x << endl << flush

real get_elapsed_time();
void get_cpu_time(real& user_time, real& system_time);

#endif

