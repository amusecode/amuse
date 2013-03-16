
// General functions and random numbers.

#include "stdinc.h"

void warning(const char *s)  {cout << s << endl << flush;}
void err_exit(const char *s) {warning(s); exit(1);}

#include <sys/time.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#endif

// Note: the first call to a timer function sets its zero point.

static bool etset = false;
static timeval et0;
real get_elapsed_time()
{
    if (!etset) {
	gettimeofday(&et0, NULL);
	etset = true;
	cout << "initialized elapsed time counter" << endl << flush;
	return 0;
    } else {
	timeval et;			// et.tv_set = time(NULL), note
	gettimeofday(&et, NULL);
	return et.tv_sec - et0.tv_sec + 1.e-6*(et.tv_usec - et0.tv_usec);
    }
}

static bool ctset = false;
static real user0, sys0;
void get_cpu_time(real& user_time, real& system_time)
{
#ifdef _WIN32
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftUser, ftKernel;
    ULONGLONG user;
    GetProcessTimes (hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser); 
    user =  (((ULONGLONG) ftUser.dwHighDateTime) << 32) + ftUser.dwLowDateTime;
    user_time = user / 10000000U; 
    user =  (((ULONGLONG) ftKernel.dwHighDateTime) << 32) + ftKernel.dwLowDateTime;
    system_time = user / 10000000U; 
    
#else
    struct rusage tim;
    getrusage(RUSAGE_SELF, &tim);

    if (!ctset) {
	user0 = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec;
	sys0 = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec;
	user_time = system_time = 0;
	ctset = true;
	cout << "initialized CPU time counter" << endl << flush;
    } else {
	user_time = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec - user0;
	system_time = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec - sys0;
    }
#endif
}

//--------------------------------------------------------------------
// Self-contained simple random number generator, based on ran2 from
// Numerical Recipes.

// The internal random seed:

static int randx = 0;

// The initial seed used to start the current sequence:

static int initial_seed = 0;

// The number of times the current sequence has been called:

static int n_rand = 0;

// srandinter:  Accept a seed to start up randunit.
//              If seed = 0, then a positive number is chosen, which
//              is unique if no other call to srandinter is made within
//              2.0 seconds of the current call (the UNIX clock is used,
//              which returns the time in seconds).
//              The return value is the actual value of the seed (either
//              the specified non-zero value, or the clock time+PID).

int srandinter(int seed, int iter)
{
    // If no positive seed is provided, produce a random value, based
    // on the system clock (resolution 1 second) and the current
    // process ID.

    // cout << "srandinter 1: seed = " << seed << endl;
    if (seed == 0) seed = (int) (time(NULL) + 1000*getpid());

    initial_seed = abs(seed);
    n_rand = 0;

    randx = -initial_seed;	// randx < 0 ==> initialize ran1 or ran2
    // cout << "srandinter 2: randx = " << randx << endl;

    for (int i = 0; i < iter; i++) randinter(0,1);
    // cout << "srandinter 3: randx = " << randx << endl;

    return initial_seed;   	// to enable the user to verify the value
                           	// of the initialization
}

int getrandinter()
{
    return randx;
}

// From Numerical Recipes, verbatim:

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.e-14
#define RNMX (1.0-EPS)

real ran2(int& idum)
{
    int j;
    int k;
    static int idum2 = 123456789;
    static int iy = 0;
    static int iv[NTAB];
    real temp;

    if (idum <= 0) {		// idum < 0 ==> initialize
	if (-(idum) < 1)
	    idum = 1;
	else
	    idum = -(idum);
	idum2 = (idum);

	for (j = NTAB+7; j >= 0; j--) {
	    k = (idum)/IQ1;
	    idum = IA1*(idum-k*IQ1) - k*IR1;
	    if (idum < 0) idum += IM1;
	    if (j < NTAB) iv[j] = idum;
	}
	iy = iv[0];
    }
    k = (idum)/IQ1;
    idum = IA1*(idum-k*IQ1) - k*IR1;
    if (idum < 0) idum += IM1;

    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;

    j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;

    if ((temp = AM*iy) > RNMX)
	return RNMX;
    else
	return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

//------------------------------------------------------------------------

// randinter:  Return a random real number within an interval [a,b]
//	       by invoking a standard random number generator.

real randinter(real a, real b) {
    if (randx == 0) srandinter(0,0);	// srandinter was never called
    n_rand++;
    // cout << "randinter: randx = " << randx << endl;
    return a + (b-a)*ran2(randx);
}

// Old version used built-ins, but conflicts with openmpi...

static real scale = 1.0/(((unsigned long)1<<31)-1.0);
real oldrandinter(real a, real b)
{
#ifdef _WIN32
    unsigned long r = rand();
#else
    unsigned long r = random();
#endif
    real x = scale*r;
    return a + (b-a)*x;
}
