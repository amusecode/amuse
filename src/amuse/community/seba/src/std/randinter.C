
//  randinter.C: generates "random" numbers within an interval [a,b]
//.............................................................................
//    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
//    version 3:  Mar 1994   Steve McMillan, adapted Numerical Recipes routines
//.............................................................................
//  non-local functions: 
//    randinter, srandinter
//.............................................................................
//
//     The results are reproducible if a seed > 0 is specified to start up the
//  random number generator. Note that the generator can be started up with a
//  seed = 0, which is converted to a "random" seed (taken from the clock,
//  expressed in seconds).
//     The "randunit" version is copied after UNIX calls to  rand() , which
//  are unfortunately not very random. An improved version should be used for
//  serious applications.
//     The ran0, ran1 and ran2 routines are said to be better. In particular,
// ran2 is claimed to be "perfect"!
//
//     Note that ran1 and ran2 maintain a table of numbers in order to deal
//  with sequential correntations in their output. This makes it more difficult
//  to reproduce random numbers at an arbitrary point within the sequence,
//  so we must maintain some history information about calls to randinter.
//
//.............................................................................
//     Compiled with cc -DTOOLBOX , randinter # will model # throws of dice.
//.............................................................................


//// Check the Starlab random number generators.  Print out a specified
//// number (N) of random numbers uniformly distributed on [0,1).
////
//// Usage: randinter N seed
////
//// Options:
////     first argument = number of "throws of the dice," N  [no default]
////     second argument = random seed  [take from system clock]
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include <time.h>
#include "stdinc.h"

#ifndef TOOLBOX

//------------------------------------------------------------------------
//
// Initialization:
// --------------

// The internal random seed:

static int randx = 0;

// The initial seed used to start the current sequence:

static int initial_seed = 0;

// The number of times the current sequence has been called:

static int n_rand = 0;

// srandinter:  Accept a seed to start up randunit.
//              note:
//                  if seed = 0, then a positive number is chosen, which
//                  is unique if no other call to srandinter is made within
//                  2.0 seconds of the current call (the UNIX clock is used,
//                  which returns the time in seconds).
//                  The return value is the actual value of the seed (either
//                  the specified non-zero value, or the clock time).

int  srandinter(int seed, int iter) {

    if (seed == 0)               // no particular positive seed provided?
	seed = (int) time(NULL); // then give a random value,
				 // different every second.

    initial_seed = abs(seed);
    n_rand = 0;

    randx = -initial_seed;	 // randx < 0 ==> initialize ran1 or ran2

    for (int i = 0; i < iter; i++) randinter(0,1);

    return initial_seed;   // to enable the user to verify the actual value of
                           // the initialization for randx.
}

// get_initial_seed:  Return the initial seed for the current sequence

int get_initial_seed() {

    return initial_seed;
}

// get_rand_seed:  Return the current value of the random seed

int get_rand_seed() {

    return abs(randx);
}

// get_current_seed:  Lookalike for get_rand_seed

int get_current_seed() {

    return abs(randx);
}

// get_n_rand:  Return the number of invocations of the current sequence

int get_n_rand() {

    return n_rand;
}

//------------------------------------------------------------------------
// Random number generators:
// ------------------------

// randunit:  Return a random real number within the unit interval
//		  note: based on      @(#)rand.c   4.1 (Berkeley) 12/21/80,
//			but returning a positive number smaller than unity.

#define  MAXNUM	2147483647.0 	// the maximum value which rand() can return

real  randunit() {

    if (randx < 0) randx *= -1; // for compatibility with ran1, ran2
    randx &= 0x7fffffff;        // to mimic 32-bit int on a Cray

    return (real)((randx = randx * 1103515245 + 12345) & 0x7fffffff)/MAXNUM;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// nr_rand: Three congruential random-number generators from Numerical Recipes.
//	    Converted to C++ and 64-bit real precision by Steve McMillan

// Each routine takes a int random seed and returns a random number
// uniformly distributed between 0 and 1. The seed is modified on return,
// ready for the next call.

// Note: This package uses "int" to mean "long int". Starlab users on 16-bit
//       machines beware!!

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

real ran0(int& idum)
{
    int k;
    real ans;

    if (idum < 0) idum *= -1;	// for compatibility with ran1, ran2
    
    idum ^= MASK;
    k = (idum)/IQ;
    idum = IA*(idum-k*IQ) - IR*k;
    if (idum < 0) idum += IM;
    ans = AM*(idum);
    idum ^= MASK;
    
    return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.e-14
#define RNMX (1.0-EPS)

real ran1(int& idum)
{
    int j;
    int k;
    static int iy = 0;
    static int iv[NTAB];
    real temp;

    if (idum <= 0 || !iy) {	// idum < 0 ==> initialize
	if (-(idum) < 1)
	    idum = 1;
	else
	    idum = -(idum);
	for (j = NTAB+7; j >= 0; j--) {
	    k = (idum)/IQ;
	    idum = IA*(idum-k*IQ) - IR*k;
	    if (idum < 0) idum += IM;
	    if (j < NTAB) iv[j] = idum;
	}
	iy = iv[0];
    }
    
    k = (idum)/IQ;
    idum = IA*(idum-k*IQ) - IR*k;
    if (idum < 0) idum += IM;

    j = iy/NDIV;
    iy = iv[j];
    iv[j] = idum;
    
    if ((temp = AM*iy) > RNMX)
	return RNMX;
    else
	return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

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

    n_rand++;

//  Old Starlab:
//
//    return a + (b-a)*randunit();

//  New Starlab:
//
    return a + (b-a)*ran2(randx);

}

// gausrand: Return a random number distributed according to a Gaussian
//	     distribution with specified mean and standard deviation.

#define BIG 10

real gausrand(real mean, real sdev)
{
    // Select a number with zero mean and unit variance using the
    // rejection method (see Numerical Recipes).

    for(;;) {
	real x = randinter(-BIG, BIG);
	if (randinter(0,1) < exp(-0.5*x*x)) return mean + x*sdev;
    }
}

//------------------------------------------------------------------------

#else

// main:  Driver to test  randinter() . It also emulates N throws of dice,
//        when given the command  "randinter N" . If no second argument
//        is given, a random seed is used; otherwise the second argument
//        is used as the seed for the random number generator.

main(int argc, char **argv) {

    check_help();
    extern char *poptarg;
    int c;
    const char *param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.4 $", _SRC_)) != -1) {}

    if (argc < 2)
	err_exit("randinter N [S] (for N throws of dice, from seed S)");

    int seed;
    if (argc == 2)
	seed = srandinter(0);                   // "random" seed
    else
        seed = srandinter(atoi(argv[2]));

    cout.precision(8);
    int i;
    for (i = 1; i <= atoi(argv[1]); i++) {
	real x = randinter(0,1);
	cout << randinter(0, 1) << " ";
	if (i%5 == 0) cout << endl;
    }
    if ((i-1)%5 != 0) cout << endl;
}
#endif
