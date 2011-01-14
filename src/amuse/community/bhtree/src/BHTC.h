
//#ifndef _MUSE_BHTC_
//#define _MUSE_BHTC_

#include<fstream>
#include <iostream>

#include <fstream>
#include <iostream>

#include <cstring>

#include "stdinc.h"
#include "vec.h"
#include "nbody_particle.h"

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef nbody_VF_ptr real_VF_ptr;
typedef nbody_RF_ptr real_RF_ptr;
typedef nbody_RRF_ptr real_RRF_ptr;

extern "C" double cpusec();
int  pgetopt(int argc, char ** argv,  char * optstr);
void pskipopt();

//#endif //_MUSE_BHTC_
