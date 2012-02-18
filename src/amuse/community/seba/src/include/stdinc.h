
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file stdinc.h  Standard Starlab include file

//----------------------------------------------------------------------

// Include Doxygen introductory documentation here.

/** \mainpage Starlab Classes and Library Functions

\section overview Overview

This is the (automatically generated) top-level page documenting
Starlab's internal structure.  It contains detailed listings of the
header files, major classes, and global functions found in the
package.  These pages are still under construction.  If there are
specific features you'd like to see covered in more detail, please
contact us at the address below.  (Note that "set" and "get" accessor
functions for member data are usually obvious, and are generally not
documented.) The goals of these pages are

<ul>
<li>to provide the application programmer with enough information to
write standalone programs using the Starlab libraries, and
<li>to document the internals well enough for programmers to
contribute new functions to the Starlab source tree.

\section structure Structure

The above links provide access to
<ul>
<li>Starlab namespace entries
<li>Starlab classes
<li>Starlab global functions
<li>Starlab header files
</ul>
Most of the listings are self-evident.  To see listings of global
functions, go to File Members and then select Functions.  For a list
of functions by file, select File List.

\section comments Comments

Send comments and suggestions to <a
href=mailto:steve@physics.drexel.edu>steve@physics.drexel.edu</a><br><br>
*/

//----------------------------------------------------------------------

//  version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
//
//     This file includes:
//  1) new naming conventions to add to or replace existing names in C
//  2) a string manipulation macro
//  3) mathematical constants
//  4) functions  abs()  min(,)  max(,)
//  5) macros to cast angular arguments in standard form

#ifndef  STARLAB_STDINC_H
#  define  STARLAB_STDINC_H

#include  <cstdlib>      // stdlib.h with namespace std
#include  <cmath>        // math.h   with namespace std
#include  <cstring>
#include  <string>
#include  <fstream>

#ifdef HAVE_CONFIG_H
#  include  "config.h"	// config.h is created by "configure" (see autoconf)
#endif

// New GNU, HP, and Sun all use stdiostream.h, which includes both
// the stdio.h (C) and the iostream.h (C++) header files.  SGI and
// old GNU apparently want iostream.h and stdio.h to be included
// separately...

#if defined(sgi) && !defined(__GNUC__) && !defined(SG_CPLUSPLUS)
# define SG_CPLUSPLUS 1
#endif

// Heck, for gcc3 this should do .... forget about the rest for now.
#include <iostream>

// May not need include <cstdio> here??
#include <cstdio>

#if 0
#if defined GNU_CPLUSPLUS
#  ifdef GNU_POST_29
#    include  <iostream.h>
#  else
#    ifdef GNU_PRE_26
#      include  <iostream.h>
#      include  <stdio.h>
#    else
#      include  <stdiostream.h>
#    endif
#  endif
#elif defined SUN_CPLUSPLUS
#  include  <stdiostream.h>
#elif defined HP_CPLUSPLUS
#  include  <stdiostream.h>
#elif defined DEC_CPLUSPLUS
#  include  <stdiostream.h>
#elif defined SG_CPLUSPLUS
#  include  <stdio.h>
#  include  <iostream.h>
#endif
#endif

// So we get the good old standard cerr, cout, cin, ....

using namespace std;

//=============================================================================
//  New naming conventions to add to or replace existing names in C :
//=============================================================================

//-----------------------------------------------------------------------------
//  real  --  a more general name for the standard floating-point data type
//-----------------------------------------------------------------------------

typedef double  real;

//-----------------------------------------------------------------------------
//  bool  --  another name for int, to indicate use in logical operations
//-----------------------------------------------------------------------------

// g++ 2.6 (and later versions) already has a "bool" data type...

#if !defined(SG_CPLUSPLUS) && !defined(__GNUG__)
    typedef int bool;
#   define  false  0
#   define  true   1
#endif

// Convenient definitions (the g++ >2.6 "bool" type has enumerated values
// "false" and "true," which may be used interchangeably with these):

#  define  FALSE  0
#  define  TRUE   1

//-----------------------------------------------------------------------------
//  local  --  a more descriptive name for variables or functions which
//             are invisible outside the file in which they are defined.
//-----------------------------------------------------------------------------

#define  local      static

//=============================================================================
//  A  string manipulation macro :
//=============================================================================

//-----------------------------------------------------------------------------
//  streq  --  a macro which returns 1 if two strings are equal, 0 otherwise
//-----------------------------------------------------------------------------

#define  streq(x,y)  (strcmp((x), (y)) == 0)

//=============================================================================
// Simple output (#param is ANSI C, but presently only works with g"++
//=============================================================================

#define PRI(x) {for (int __pri__ = 0; __pri__ < x; __pri__++) cerr << " ";}

//#ifdef GNU_CPLUSPLUS

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << endl

//#else

//#define PR(x)  cerr << "x" << " = " << x << " "
//#define PRC(x) cerr << "x" << " = " << x << ",  "
//#define PRL(x) cerr << "x" << " = " << x << endl

//#endif

//=============================================================================
//  Mathematical constants : 
//=============================================================================

//-----------------------------------------------------------------------------
//  pi, etc.  --  mathematical constants, as well as `infinity'
//-----------------------------------------------------------------------------

#ifndef PI
#  define   PI         3.14159265358979323846
#endif
#ifndef M_PI
#  define   M_PI	PI
#endif
#define   TWO_PI     (2 * (PI))
#define   HALF_PI    (0.5 * (PI))
#define   ONE_THIRD  0.33333333333333333333
#define   ONE_SIXTH  0.16666666666666666667

#define VERY_LARGE_NUMBER 1.0e300
#define VERY_SMALL_NUMBER (pow(2.0, -52)) // dynamic limit on powers of
					  // 2 for double precision

#define VERY_LARGE_INTEGER (1<<30)	  // assume we are working with
					  // standard 4-byte integers, but
					  // we also want the quantity
					  // -VERY_LARGE_INTEGER to be legal

#define SETBIT(i, n) ((i) |= (1 << (n)))  // set bit n of i
#define GETBIT(i, n) ((i) & (1 << (n)))   // get bit n of i
#define CLRBIT(i, n) ((i) &= ~(1 << (n))) // clear bit n of i
void printbits(unsigned int i);		  // print significant nonzero bits

//=============================================================================
//  Other constants : 
//=============================================================================

#define LOW_PRECISION	 3
#define STD_PRECISION	 6
#define INT_PRECISION	10
#define HIGH_PRECISION	15

//=============================================================================
//  Functions  abs()  min(,)  max(,)  square()
//=============================================================================
//
//-----------------------------------------------------------------------------
//  abs  --  returns the absolute value of its argument
//  max  --  returns the argument with the highest value
//  min  --  returns the argument with the lowest value
//  square - returns the square of its argument (cf. vec version)
//-----------------------------------------------------------------------------

namespace Starlab {
#ifndef HP_CPLUSPLUS
  inline real      abs(real x)                	{return (x < 0) ? -x : x;}
  inline long int  abs(long int x)            	{return (x < 0) ? -x : x;}
#   ifndef GNU_CPLUSPLUS
#   ifndef SG_CPLUSPLUS
      inline int       abs(int x)             	{return (x < 0) ? -x : x;}
#   endif
#   endif
#endif
}

#ifdef _WIN32
  extern double erf(double);
  extern double erfc(double);

  inline real      rint(real x)			{return floor((x) + 0.5);}
#endif

#undef max
#undef min

namespace Starlab {

  inline real		max(real x, real y)	{return (x > y) ? x : y;}
  inline int		max(int x, int y)	{return (x > y) ? x : y;}
  inline long int	max(long int x, long int y)  {return (x > y) ? x : y;}

  inline real		min(real x, real y)	{return (x < y) ? x : y;}
  inline int		min(int x, int y)	{return (x < y) ? x : y;}
  inline long int	min(long int x, long int y)  {return (x < y) ? x : y;}

// Force mixed "max" and "min" calculations to be done as real, to
// avoid "ambiguous" errors from g++

#ifdef GNU_PRE_26
  inline real	max(int x, real y)         	{return (x > y) ? x : y;}
  inline real	max(real x, int y)         	{return (x > y) ? x : y;}

  inline real	min(real x, int y)         	{return (x < y) ? x : y;}
  inline real	min(int x, real y)         	{return (x < y) ? x : y;}
#endif
}

inline real	square(real x)			{return x*x;}

//=============================================================================
//  Macros to cast angular arguments in standard form :
//=============================================================================

//-----------------------------------------------------------------------------
//  pos_angle  --  recasts an angular variable into the range [0, TWO_PI)
//  sym_angle  --  recasts an angular variable into the range [-PI, PI)
//                   (recasting: transforming modulo 2 pi)
//         example:
//                 to map an angular variable 'phi' into the smallest positive
//                 value, use
//
//                     phi = pos_angle(phi);
//
//                 to map an angular variable 'phi' into the smallest value,
//                 positive or negative, use
//
//                     phi = sym_angle(phi);
//
//-----------------------------------------------------------------------------

#define  pos_angle(phi)    ((phi) - TWO_PI * floor((phi)/TWO_PI ))
#define  sym_angle(phi)    ((phi) - TWO_PI * floor(((phi)+PI)/TWO_PI ))

//=============================================================================
// xreal -- extended precision real, for use with time and system_time
//=============================================================================

// Long long is an ongoing irritant...  See std/xreal.C
// *** Options should now be properly set up by configure. (Steve, 9/04) ***

#if 0			// old (pre-/04)

#  ifdef HAVE_STRTOLL

//    e.g. Linux.

#     define STRTOL strtoll
#     define STRTOUL strtoull

#  else

//    e.g. Dec UNIX

#     define STRTOL strtol
#     define STRTOUL strtoul

#  endif

#else			// new

// Configure should have done all the work, but just in case the size of
// a long long isn't 8 bytes, silently disable USE_XREAL.  (If this has
// been handled properly, configure should already have warned the user.)

#  ifdef USE_XREAL

//    Make sure...

#     ifdef HAVE_LONG_LONG
#        ifdef SIZEOF_LONG_LONG
#           if (SIZEOF_LONG_LONG != 8)
//#              warning "Undefining USE_XREAL because long long size is not 8"
#              undef USE_XREAL
#           endif
#        endif
#     else
//#        warning "Undefining USE_XREAL because long long is not defined"
#        undef USE_XREAL
#     endif
#
#   endif

#endif

#ifdef USE_XREAL
#  include "xreal.h"
#else
   typedef real xreal;
#endif

#ifndef HAVE_STPCPY
  char *stpcpy(char *restrict0, const char *restrict1);
#endif

// Intended for xreal, but referenced in the real version too.

void xprint(xreal x,
	    ostream & s = cerr,
	    bool newline = true);
real fmod2(xreal x, real y);
xreal get_xreal(char *str);
void identify_xreal(ostream& s = cerr);

//=============================================================================
//  Various declarations
//=============================================================================

void print_message(const char*);
void warning(const char*);
void err_exit(const char*);

int pgetopt(int argc, char** argv, const char *optstr,
	    const char *cvs_id = NULL, const char *source = NULL);

void pskipopt();
void params_to_usage(ostream&, char*, const char*);

int  srandinter(int, int n = 0);
int  get_initial_seed();
int  get_rand_seed();
int  get_current_seed();
int  get_n_rand();
real randinter(real, real);
real gausrand(real, real);

void cpu_init();
real cpu_time();

void starlab_wait(int iwait);

int set_starlab_precision(ostream&);
int adjust_starlab_precision(int p);
int  get_starlab_precision();

char * gethist(int, char **);

// Convenient invocation of run-time help function.

void check_runtime_help(int argc, char** argv,
			const char* source_file,
			const char* date, const char *time);
void get_runtime_help(const char* source_file,
		      const char* date, const char *time, int level = 1);

// Note: assuming here that the macros __DATE__ and __TIME__ are standard...
// Macro _SRC_ will be provided by configure/make.

#define check_help() check_runtime_help(argc, argv, _SRC_, __DATE__, __TIME__);
#define get_help() get_runtime_help(_SRC_, __DATE__, __TIME__);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/stdinc.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
