#include <iostream>
#include <cstdlib>
#include "nstab.h"
#include "f2c.h"
using namespace std;

// Read the command line and pass the first 7 arguments to nstab.

int main(int argc, char *argv[])
{
    if (argc > 7) {
	double sigma = atof(argv[1]);
	double ei0 = atof(argv[2]);
	double eo = atof(argv[3]);
	double relinc = atof(argv[4]);
	double m1 = atof(argv[5]);
	double m2 = atof(argv[6]);
	double m3 = atof(argv[7]);
	cout << nstab_(&sigma, &ei0, &eo, &relinc, &m1, &m2, &m3) << endl;
    }
    return 0;
}
