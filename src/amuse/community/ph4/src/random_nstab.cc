#include "stdinc.h"

#include <time.h>

int main(int argc, char *argv[])
{
    int n = 7;
    if (argc > 1) n = atoi(argv[1]);

    int seed = (int) (time(NULL) + 1000*getpid());
    srandom(seed);

    // Print n random doubles, chosen from appropriate ranges.

    cout << randinter(1,5) << " "		// sigma
	 << randinter(0,1) << " "		// ei0
	 << randinter(0,1) << " "		// eo
	 << randinter(-1,1) << " "		// relinc
	 << randinter(1,5) << " "		// m1
	 << randinter(1,5) << " "		// m2
	 << randinter(1,5) << " "		// m3
	 << endl;
}
