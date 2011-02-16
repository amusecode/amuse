#include "stdinc.h"
#include "jdata.h"

// Debugging functions:

void printq(int j, real2 q, const char *s);
void print_particle(int j, jdata &jd);
void print_list(int jlist[], int nlist, jdata &jd);
void print_list(jdata &jd);
bool twiddles(real q, real v, real tol = 1.e-3);
real total_energy(int jlist[], int n, jdata& jd);
